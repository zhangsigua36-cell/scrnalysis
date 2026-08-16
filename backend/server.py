"""Local Scrnalysis analysis service.

This small standard-library server accepts a 10x Matrix Market zip or the
three matrix files, runs the checked-in R workflow, and serves its outputs.
It is intentionally local-first: the browser and this service should run on
the same computer, while the browser may be opened from another device on
the same network.
"""

from __future__ import annotations

import gzip
import json
import mimetypes
import os
import secrets
import shutil
import subprocess
import threading
import time
import zipfile
from email.parser import BytesParser
from email.policy import default
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse


ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = ROOT / "runs"
ANALYSIS_SCRIPT = ROOT / "test_runs" / "GSE279818_test_20260816" / "scripts" / "run_test_analysis.R"
RSCRIPT = Path(os.environ.get("SCRNALYSIS_RSCRIPT", r"D:\R\R-4.5.2\bin\Rscript.exe"))
MAX_UPLOAD_BYTES = 1024 * 1024 * 1024
JOBS: dict[str, dict] = {}
JOBS_LOCK = threading.Lock()
ANALYSIS_SLOT = threading.Semaphore(1)


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def safe_filename(value: str) -> str:
    name = Path(value).name
    return "".join(char for char in name if char.isalnum() or char in ".-_()") or "upload.bin"


def is_inside(parent: Path, child: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def parse_multipart(handler: BaseHTTPRequestHandler) -> tuple[list[tuple[str, str, bytes]], dict[str, str]]:
    content_type = handler.headers.get("Content-Type", "")
    if not content_type.startswith("multipart/form-data"):
        raise ValueError("请上传文件表单")
    length = int(handler.headers.get("Content-Length", "0"))
    if length <= 0 or length > MAX_UPLOAD_BYTES:
        raise ValueError("上传文件为空，或超过 1 GB 限制")
    body = handler.rfile.read(length)
    raw = (f"Content-Type: {content_type}\r\nMIME-Version: 1.0\r\n\r\n").encode() + body
    message = BytesParser(policy=default).parsebytes(raw)
    files: list[tuple[str, str, bytes]] = []
    fields: dict[str, str] = {}
    for part in message.iter_parts():
        field_name = part.get_param("name", header="content-disposition") or ""
        filename = part.get_filename()
        payload = part.get_payload(decode=True) or b""
        if filename:
            files.append((field_name, safe_filename(filename), payload))
        else:
            fields[field_name] = payload.decode("utf-8", errors="replace")
    if not files:
        raise ValueError("没有收到数据文件")
    return files, fields


def extract_uploads(files: list[tuple[str, str, bytes]], input_dir: Path) -> None:
    input_dir.mkdir(parents=True, exist_ok=True)
    for _, filename, payload in files:
        destination = input_dir / filename
        if filename.lower().endswith(".zip"):
            zip_path = input_dir / filename
            zip_path.write_bytes(payload)
            with zipfile.ZipFile(zip_path) as archive:
                for member in archive.infolist():
                    target = (input_dir / member.filename).resolve()
                    if member.is_dir():
                        continue
                    if not is_inside(input_dir, target):
                        raise ValueError("压缩包包含不安全的路径")
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(member) as source, target.open("wb") as target_file:
                        shutil.copyfileobj(source, target_file)
            zip_path.unlink(missing_ok=True)
        else:
            destination.write_bytes(payload)


def find_named(input_dir: Path, names: tuple[str, ...]) -> Path | None:
    candidates = [path for path in input_dir.rglob("*") if path.is_file() and path.name.lower() in names]
    return candidates[0] if candidates else None


def normalize_10x_input(input_dir: Path) -> None:
    expected = {
        "matrix.mtx.gz": ("matrix.mtx.gz", "matrix.mtx"),
        "barcodes.tsv.gz": ("barcodes.tsv.gz", "barcodes.tsv"),
        "features.tsv.gz": ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"),
    }
    for target_name, source_names in expected.items():
        source = find_named(input_dir, tuple(name.lower() for name in source_names))
        if source is None:
            raise ValueError("需要 10x 三个文件：matrix.mtx.gz、barcodes.tsv.gz、features.tsv.gz")
        target = input_dir / target_name
        if source.resolve() == target.resolve():
            continue
        if source.name.lower().endswith(".gz"):
            shutil.copyfile(source, target)
        else:
            with source.open("rb") as source_file, gzip.open(target, "wb") as target_file:
                shutil.copyfileobj(source_file, target_file)


def create_contract(run_dir: Path, config: dict) -> None:
    metadata_dir = run_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    write_json(metadata_dir / "user_config.json", config)
    yaml_lines = [
        f"project_id: {config['job_id']}",
        "mode: uploaded_local_run",
        "raw_data_read_only: true",
        "user_fixed:",
        f"  background: {json.dumps(config.get('analysisBackground', ''), ensure_ascii=False)}",
        f"  purpose: {json.dumps(config.get('researchPurpose', ''), ensure_ascii=False)}",
        f"  comparison_groups: {json.dumps(config.get('comparisonGroups', ''), ensure_ascii=False)}",
        "parameters:",
    ]
    for key, value in config.items():
        if key in {"job_id", "analysisBackground", "researchPurpose", "comparisonGroups"}:
            continue
        yaml_lines.append(f"  {key}: {json.dumps(value, ensure_ascii=False)}")
    (metadata_dir / "run_contract.yaml").write_text("\n".join(yaml_lines) + "\n", encoding="utf-8")


def update_job(job_id: str, **changes: object) -> None:
    with JOBS_LOCK:
        JOBS[job_id].update(changes)


def public_job(job_id: str) -> dict:
    with JOBS_LOCK:
        job = dict(JOBS[job_id])
    job.pop("run_dir", None)
    job.pop("process", None)
    if job.get("status") == "complete":
        base = f"/api/jobs/{job_id}/files/"
        job["artifacts"] = {
            "umap": base + "results/figures/umap_by_cluster.png",
            "tsne": base + "results/figures/tsne_by_cluster.png",
            "qc": base + "results/figures/qc_after_filtering.png",
            "report": base + "reports/final_analysis_report.md",
            "umapCoordinates": base + "results/tables/umap_coordinates.csv",
            "markers": base + "results/tables/cluster_markers.csv",
            "qcSummary": base + "results/tables/qc_summary.csv",
        }
    return job


def run_analysis(job_id: str, run_dir: Path, input_dir: Path, config: dict) -> None:
    acquired = ANALYSIS_SLOT.acquire(blocking=True)
    try:
        update_job(job_id, status="running", progress=5, message="正在启动单细胞分析引擎", started_at=now())
        log_path = run_dir / "logs" / "backend.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        args = [
            str(RSCRIPT), str(ANALYSIS_SCRIPT),
            "--run-dir", str(run_dir), "--raw-dir", str(input_dir), "--dataset-id", job_id,
            "--species", config["species"], "--tissue", config["tissue"],
            "--normalization", config["normalization"], "--batch-correction", config["batchCorrection"],
            "--gene-filter-mode", config["geneFilterMode"], "--gene-lower-quantile", str(config["geneLowerQuantile"]),
            "--gene-upper-quantile", str(config["geneUpperQuantile"]), "--min-genes", str(config["minGenes"]), "--max-genes", str(config["maxGenes"]),
            "--umi-filter-mode", config["umiFilterMode"], "--umi-lower-quantile", str(config["umiLowerQuantile"]),
            "--umi-upper-quantile", str(config["umiUpperQuantile"]), "--min-counts", str(config["minCounts"]), "--max-counts", str(config["maxCounts"]),
            "--mito-filter-mode", config["mitoFilterMode"], "--max-mito", str(config["maxMito"]),
            "--ribo-filter-mode", config["riboFilterMode"], "--max-ribo", str(config["maxRibo"]),
            "--hb-filter-mode", config["hbFilterMode"], "--max-hb", str(config["maxHb"]),
            "--min-cells-per-gene-mode", config["minCellsPerGeneMode"], "--min-cells-per-gene", str(config["minCellsPerGene"]),
            "--dims", str(config["dims"]), "--resolution", str(config["resolution"]),
            "--neighbors", str(config["neighbors"]), "--umap-neighbors", str(config["umapNeighbors"]),
            "--umap-min-dist", str(config["umapMinDist"]), "--tsne-perplexity", str(config["tsnePerplexity"]),
            "--doublet-method", config["doubletMethod"], "--ambient-method", config["ambientMethod"],
        ]
        with log_path.open("w", encoding="utf-8") as log:
            log.write(f"{now()} | command started\n")
            process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
            update_job(job_id, process=process.pid)
            for line in process.stdout or []:
                log.write(line)
                log.flush()
                if "STEP 1/9" in line:
                    update_job(job_id, progress=15, message="正在读取 10x 矩阵")
                elif "STEP 2/9" in line or "STEP 3/9" in line:
                    update_job(job_id, progress=28, message="正在进行 QC 质控")
                elif "STEP 4/9" in line:
                    update_job(job_id, progress=40, message="正在检查 doublet")
                elif "STEP 6/9" in line:
                    update_job(job_id, progress=58, message="正在标准化并计算 PCA")
                elif "STEP 7/9" in line:
                    update_job(job_id, progress=74, message="正在生成聚类和 UMAP")
                elif "STEP 8/9" in line:
                    update_job(job_id, progress=88, message="正在导出 t-SNE 和 marker")
            return_code = process.wait()
        if return_code != 0:
            update_job(job_id, status="error", progress=100, message="分析未完成，请查看错误信息", error=f"Rscript exit code {return_code}", finished_at=now())
            return
        if not (run_dir / "results" / "figures" / "umap_by_cluster.png").exists():
            update_job(job_id, status="error", progress=100, message="分析结束但没有找到 UMAP 图", error="UMAP output missing", finished_at=now())
            return
        update_job(job_id, status="complete", progress=100, message="分析完成，已生成 UMAP 和结果文件", finished_at=now())
    except Exception as exc:  # keep the browser-facing service alive for later jobs
        update_job(job_id, status="error", progress=100, message="分析服务遇到问题", error=str(exc), finished_at=now())
    finally:
        if acquired:
            ANALYSIS_SLOT.release()


class Handler(BaseHTTPRequestHandler):
    server_version = "ScrnalysisLocal/0.1"

    def end_headers(self) -> None:
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        super().end_headers()

    def do_OPTIONS(self) -> None:
        self.send_response(HTTPStatus.NO_CONTENT)
        self.end_headers()

    def send_json(self, payload: object, status: int = 200) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        path = unquote(urlparse(self.path).path)
        if path == "/api/health":
            self.send_json({"ok": True, "rscript": str(RSCRIPT), "rscript_exists": RSCRIPT.exists()})
            return
        prefix = "/api/jobs/"
        if path.startswith(prefix):
            rest = path[len(prefix):]
            job_id, _, file_path = rest.partition("/files/")
            with JOBS_LOCK:
                exists = job_id in JOBS
            if not exists:
                self.send_json({"error": "任务不存在"}, 404)
                return
            if not file_path:
                self.send_json(public_job(job_id))
                return
            run_dir = Path(JOBS[job_id]["run_dir"])
            target = (run_dir / file_path).resolve()
            if not is_inside(run_dir, target) or not target.is_file():
                self.send_json({"error": "结果文件不存在"}, 404)
                return
            body = target.read_bytes()
            content_type = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Content-Disposition", f'inline; filename="{target.name}"')
            self.end_headers()
            self.wfile.write(body)
            return
        self.send_json({"error": "not found"}, 404)

    def do_POST(self) -> None:
        if self.path.split("?", 1)[0] != "/api/jobs":
            self.send_json({"error": "not found"}, 404)
            return
        try:
            files, fields = parse_multipart(self)
            config = json.loads(fields.get("config", "{}"))
            job_id = time.strftime("%Y%m%d_%H%M%S") + "_" + secrets.token_hex(3)
            run_dir = RUNS_DIR / job_id
            input_dir = run_dir / "input"
            run_dir.mkdir(parents=True, exist_ok=True)
            extract_uploads(files, input_dir)
            normalize_10x_input(input_dir)
            config = {
                "job_id": job_id,
                "analysisBackground": str(config.get("analysisBackground", "")),
                "researchPurpose": str(config.get("researchPurpose", "")),
                "comparisonGroups": str(config.get("comparisonGroups", "")),
                "species": str(config.get("species", "human")),
                "tissue": str(config.get("tissue", "synovium")),
                "normalization": str(config.get("normalization", "log")),
                "batchCorrection": str(config.get("batchCorrection", "harmony")),
                "geneFilterMode": str(config.get("geneFilterMode", "quantile")),
                "geneLowerQuantile": float(config.get("geneLowerQuantile", 0.025)),
                "geneUpperQuantile": float(config.get("geneUpperQuantile", 0.025)),
                "minGenes": int(config.get("minGenes", 200)),
                "maxGenes": int(config.get("maxGenes", 6000)),
                "umiFilterMode": str(config.get("umiFilterMode", "quantile")),
                "umiLowerQuantile": float(config.get("umiLowerQuantile", 0.025)),
                "umiUpperQuantile": float(config.get("umiUpperQuantile", 0.025)),
                "minCounts": int(config.get("minCounts", 0)),
                "maxCounts": int(config.get("maxCounts", 25000)),
                "mitoFilterMode": str(config.get("mitoFilterMode", "fixed")),
                "maxMito": float(config.get("maxMito", 20)),
                "riboFilterMode": str(config.get("riboFilterMode", "none")),
                "maxRibo": float(config.get("maxRibo", 30)),
                "hbFilterMode": str(config.get("hbFilterMode", "none")),
                "maxHb": float(config.get("maxHb", 20)),
                "minCellsPerGeneMode": str(config.get("minCellsPerGeneMode", "fixed")),
                "minCellsPerGene": int(config.get("minCellsPerGene", 3)),
                "dims": int(config.get("dims", 30)),
                "resolution": float(config.get("resolution", 0.6)),
                "neighbors": int(config.get("neighbors", 30)),
                "umapNeighbors": int(config.get("umapNeighbors", 30)),
                "umapMinDist": float(config.get("umapMinDist", 0.3)),
                "tsnePerplexity": int(config.get("tsnePerplexity", 30)),
                "doubletMethod": str(config.get("doubletMethod", "scDblFinder")),
                "ambientMethod": str(config.get("ambientMethod", "diagnostic_only")),
            }
            create_contract(run_dir, config)
            with JOBS_LOCK:
                JOBS[job_id] = {"job_id": job_id, "status": "queued", "progress": 2, "message": "文件已上传，等待分析", "created_at": now(), "run_dir": str(run_dir)}
            threading.Thread(target=run_analysis, args=(job_id, run_dir, input_dir, config), daemon=True).start()
            self.send_json(public_job(job_id), 202)
        except (ValueError, json.JSONDecodeError, zipfile.BadZipFile) as exc:
            self.send_json({"error": str(exc)}, 400)
        except Exception as exc:
            self.send_json({"error": f"服务器无法创建任务：{exc}"}, 500)


def main() -> None:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    if not RSCRIPT.exists():
        print(f"Warning: Rscript not found at {RSCRIPT}")
    server = ThreadingHTTPServer(("0.0.0.0", int(os.environ.get("SCRNALYSIS_PORT", "8000"))), Handler)
    print(f"Scrnalysis analysis service listening on http://0.0.0.0:{server.server_address[1]}")
    server.serve_forever()


if __name__ == "__main__":
    main()
