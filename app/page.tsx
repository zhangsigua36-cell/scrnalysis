"use client";

import { ChangeEvent, DragEvent, useEffect, useMemo, useRef, useState } from "react";

type RunStatus = "ready" | "running" | "complete" | "error";
type ViewMode = "celltype" | "sample" | "cluster";
type EmbeddingMode = "umap" | "tsne";
type SidePanel = "workspace" | "runs" | "markers" | "files";
type OutputFigure = "umap" | "tsne" | "dotplot" | "featureplot" | "qc_violin";

type AnnotationRow = {
  cluster: string;
  celltype_major: string;
  markers: string;
  score: number;
  confidence: string;
};

type MarkerDisplayRow = {
  cluster: string;
  label: string;
  markers: string;
  score: string;
  confidence: string;
};

type RunRecord = {
  job_id: string;
  created_at: string;
  status: RunStatus;
  message: string;
  files: string[];
};

type AnalysisResult = {
  job_id: string;
  status: RunStatus;
  progress: number;
  message: string;
  error?: string;
  artifacts?: Record<string, string>;
  annotationRows?: AnnotationRow[];
  insight?: { summary: string; confidence: number };
};

type ValidationErrors = Partial<Record<"files" | "background" | "purpose" | "groups" | "species" | "tissue" | "normalization" | "batch" | "annotation" | "outputs", string>>;

const API_BASE = process.env.NEXT_PUBLIC_SCRNALYSIS_API || (typeof window !== "undefined" ? `http://${window.location.hostname}:8000` : "http://localhost:8000");

const pipeline = [
  ["01", "上传数据", "读取矩阵与样本信息"],
  ["02", "背景与目的", "明确研究问题与比较分组"],
  ["03", "分析参数", "QC / 降维 / 整合参数"],
  ["04", "细胞注释", "marker 证据与 AI 建议"],
  ["05", "输出文件", "选择需要保存的图形"],
] as const;

const cellTypes = [
  { name: "T cells", count: "2,418", color: "#e88663" },
  { name: "Myeloid", count: "1,967", color: "#5b8def" },
  { name: "B cells", count: "1,104", color: "#58b88e" },
  { name: "Stromal", count: "736", color: "#d6a14c" },
  { name: "Epithelial", count: "428", color: "#a577e8" },
];

const markerRows = [
  ["T cells", "IL7R, LTB, TRBC2", "0.94", "高"],
  ["Myeloid", "LYZ, FCER1G, CTSS", "0.91", "高"],
  ["B cells", "CD79A, MS4A1, CD37", "0.88", "高"],
  ["Stromal", "COL1A1, DCN, COL3A1", "0.82", "中"],
];

const markerLibrary = [
  ["T cells", "CD3D, CD3E, TRBC1, IL7R, LTB"],
  ["NK cells", "NKG7, GNLY, FCGR3A"],
  ["Myeloid", "LYZ, LST1, FCER1G, CTSS"],
  ["B cells", "CD79A, MS4A1, CD74, CD37"],
  ["Plasma cells", "MZB1, JCHAIN, SDC1"],
  ["Fibroblast / stromal", "COL1A1, COL3A1, DCN, LUM"],
  ["Endothelial", "PECAM1, VWF, KDR"],
  ["Epithelial", "EPCAM, KRT8, KRT18, KRT19"],
];

function createPoints() {
  const groups = [
    { color: "#e88663", cx: 27, cy: 31, spread: 13, count: 28 },
    { color: "#5b8def", cx: 68, cy: 34, spread: 15, count: 25 },
    { color: "#58b88e", cx: 44, cy: 69, spread: 12, count: 20 },
    { color: "#d6a14c", cx: 76, cy: 72, spread: 9, count: 13 },
    { color: "#a577e8", cx: 18, cy: 74, spread: 8, count: 9 },
  ];

  return groups.flatMap((group, groupIndex) =>
    Array.from({ length: group.count }, (_, index) => {
      const angle = (index * 2.399 + groupIndex) % (Math.PI * 2);
      const radius = group.spread * (0.28 + ((index * 17 + groupIndex * 11) % 100) / 150);
      return {
        left: group.cx + Math.cos(angle) * radius,
        top: group.cy + Math.sin(angle) * radius * 0.72,
        color: group.color,
        size: 5 + ((index + groupIndex) % 3),
      };
    }),
  );
}

export default function Home() {
  const fileInput = useRef<HTMLInputElement>(null);
  const [selectedFiles, setSelectedFiles] = useState<File[]>([]);
  const [runStatus, setRunStatus] = useState<RunStatus>("ready");
  const [analysisResult, setAnalysisResult] = useState<AnalysisResult | null>(null);
  const [errorMessage, setErrorMessage] = useState("");
  const [activeStep, setActiveStep] = useState(0);
  const [viewMode, setViewMode] = useState<ViewMode>("celltype");
  const [embeddingMode, setEmbeddingMode] = useState<EmbeddingMode>("umap");
  const [sidePanel, setSidePanel] = useState<SidePanel>("workspace");
  const [projects, setProjects] = useState(["RA_synovium_demo", "SLE_PBMC_01", "scRNA_exploration"]);
  const [projectName, setProjectName] = useState("RA_synovium_demo");
  const [runHistory, setRunHistory] = useState<RunRecord[]>([]);
  const [savedFileNames, setSavedFileNames] = useState<string[]>([]);
  const [markerQuery, setMarkerQuery] = useState("");
  const [aiAnnotationEnabled, setAiAnnotationEnabled] = useState(true);
  const [manualEditsEnabled, setManualEditsEnabled] = useState(true);
  const [manualLabels, setManualLabels] = useState<Record<string, string>>({});
  const [selectedOutputs, setSelectedOutputs] = useState<OutputFigure[]>(["umap", "tsne", "dotplot", "featureplot", "qc_violin"]);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [maxMito, setMaxMito] = useState(20);
  const [minGenes, setMinGenes] = useState(200);
  const [maxGenes, setMaxGenes] = useState(6000);
  const [geneFilterMode, setGeneFilterMode] = useState("quantile");
  const [geneLowerQuantile, setGeneLowerQuantile] = useState(0.025);
  const [geneUpperQuantile, setGeneUpperQuantile] = useState(0.025);
  const [dims, setDims] = useState(30);
  const [resolution, setResolution] = useState(0.6);
  const [minCounts, setMinCounts] = useState(0);
  const [maxCounts, setMaxCounts] = useState(25000);
  const [umiFilterMode, setUmiFilterMode] = useState("quantile");
  const [umiLowerQuantile, setUmiLowerQuantile] = useState(0.025);
  const [umiUpperQuantile, setUmiUpperQuantile] = useState(0.025);
  const [mitoFilterMode, setMitoFilterMode] = useState("fixed");
  const [riboFilterMode, setRiboFilterMode] = useState("none");
  const [hbFilterMode, setHbFilterMode] = useState("none");
  const [maxHb, setMaxHb] = useState(20);
  const [minCellsPerGeneMode, setMinCellsPerGeneMode] = useState("fixed");
  const [minCellsPerGene, setMinCellsPerGene] = useState(3);
  const [maxRibo, setMaxRibo] = useState(30);
  const [doubletMethod, setDoubletMethod] = useState("scDblFinder");
  const [ambientMethod, setAmbientMethod] = useState("decontX");
  const [neighbors, setNeighbors] = useState(30);
  const [umapNeighbors, setUmapNeighbors] = useState(30);
  const [umapMinDist, setUmapMinDist] = useState(0.3);
  const [tsnePerplexity, setTsnePerplexity] = useState(30);
  const [species, setSpecies] = useState("human");
  const [tissue, setTissue] = useState("synovium");
  const [normalization, setNormalization] = useState("log");
  const [batchCorrection, setBatchCorrection] = useState("harmony");
  const [analysisBackground, setAnalysisBackground] = useState("");
  const [researchPurpose, setResearchPurpose] = useState("");
  const [comparisonGroups, setComparisonGroups] = useState("");
  const [annotationScope, setAnnotationScope] = useState("major");
  const [validationErrors, setValidationErrors] = useState<ValidationErrors>({});
  const points = useMemo(() => createPoints(), []);

  useEffect(() => {
    try {
      const savedRuns = window.localStorage.getItem("scrnalysis.runHistory");
      const savedProjects = window.localStorage.getItem("scrnalysis.projects");
      const savedFiles = window.localStorage.getItem("scrnalysis.fileNames");
      const savedLabels = window.localStorage.getItem("scrnalysis.manualLabels");
      if (savedRuns) setRunHistory(JSON.parse(savedRuns) as RunRecord[]);
      if (savedProjects) setProjects(JSON.parse(savedProjects) as string[]);
      if (savedFiles) setSavedFileNames(JSON.parse(savedFiles) as string[]);
      if (savedLabels) setManualLabels(JSON.parse(savedLabels) as Record<string, string>);
    } catch {
      // Local browser history is optional and should never block analysis.
    }
  }, []);

  const clearValidationError = (field: keyof ValidationErrors) => {
    setValidationErrors((current) => {
      if (!current[field]) return current;
      const next = { ...current };
      delete next[field];
      return next;
    });
  };

  const chooseFiles = (files?: FileList | File[]) => {
    if (!files?.length) return;
    setSelectedFiles(Array.from(files));
    setRunStatus("ready");
    setAnalysisResult(null);
    setErrorMessage("");
    setActiveStep(0);
    const names = Array.from(files).map((file) => file.name);
    setSavedFileNames((current) => {
      const next = Array.from(new Set([...names, ...current])).slice(0, 20);
      window.localStorage.setItem("scrnalysis.fileNames", JSON.stringify(next));
      return next;
    });
    clearValidationError("files");
  };

  const onFileChange = (event: ChangeEvent<HTMLInputElement>) => chooseFiles(event.target.files ?? undefined);
  const onDrop = (event: DragEvent<HTMLDivElement>) => {
    event.preventDefault();
    chooseFiles(event.dataTransfer.files);
  };

  const pollJob = async (jobId: string) => {
    const response = await fetch(`${API_BASE}/api/jobs/${jobId}`);
    const job = (await response.json()) as AnalysisResult;
    setAnalysisResult(job);
    setActiveStep(Math.min(pipeline.length - 1, Math.floor((job.progress || 0) / 20)));
    if (job.status === "complete") {
      setRunStatus("complete");
      setActiveStep(pipeline.length - 1);
      const record: RunRecord = { job_id: job.job_id, created_at: new Date().toISOString(), status: job.status, message: job.message, files: selectedFiles.map((file) => file.name) };
      setRunHistory((current) => {
        const next = [record, ...current.filter((item) => item.job_id !== record.job_id)].slice(0, 20);
        window.localStorage.setItem("scrnalysis.runHistory", JSON.stringify(next));
        return next;
      });
      return;
    }
    if (job.status === "error") {
      setRunStatus("error");
      setErrorMessage(job.error || job.message || "分析服务遇到问题");
      return;
    }
    window.setTimeout(() => void pollJob(jobId), 1200);
  };

  const startAnalysis = async () => {
    const nextErrors: ValidationErrors = {};
    if (!selectedFiles.length) nextErrors.files = "请先上传 10x 三个文件，或上传包含这三个文件的 zip 压缩包。";
    if (!analysisBackground.trim()) nextErrors.background = "请填写研究背景。";
    if (!researchPurpose.trim()) nextErrors.purpose = "请填写研究目的。";
    if (!comparisonGroups.trim()) nextErrors.groups = "请填写比较分组；如果没有对照组，请填写“单组探索”。";
    if (!species) nextErrors.species = "请选择物种。";
    if (!tissue) nextErrors.tissue = "请选择组织。";
    if (!normalization) nextErrors.normalization = "请选择标准化方法。";
    if (!batchCorrection) nextErrors.batch = "请选择批次校正方式。";
    if (!annotationScope) nextErrors.annotation = "请选择注释粒度。";
    if (!selectedOutputs.length) nextErrors.outputs = "请至少选择一种输出图形。";
    setValidationErrors(nextErrors);
    if (Object.keys(nextErrors).length) {
      setRunStatus("ready");
      setErrorMessage("请先补充红色星号标记的必填信息。");
      return;
    }
    setActiveStep(0);
    setRunStatus("running");
    setErrorMessage("");
    setAnalysisResult(null);
    const formData = new FormData();
    selectedFiles.forEach((file) => formData.append("files", file, file.name));
    formData.append("config", JSON.stringify({ species, tissue, normalization, batchCorrection, analysisBackground, researchPurpose, comparisonGroups, annotationScope, aiAnnotationEnabled, manualEditsEnabled, outputFigures: selectedOutputs, mitoFilterMode, maxMito, geneFilterMode, geneLowerQuantile, geneUpperQuantile, minGenes, maxGenes, umiFilterMode, umiLowerQuantile, umiUpperQuantile, minCounts, maxCounts, riboFilterMode, maxRibo, hbFilterMode, maxHb, minCellsPerGeneMode, minCellsPerGene, doubletMethod, ambientMethod, dims, resolution, neighbors, umapNeighbors, umapMinDist, tsnePerplexity }));
    try {
      const response = await fetch(`${API_BASE}/api/jobs`, { method: "POST", body: formData });
      const body = await response.json();
      if (!response.ok) throw new Error(body.error || "上传失败");
      await pollJob(body.job_id as string);
    } catch (error) {
      setRunStatus("error");
      setErrorMessage(error instanceof Error ? error.message : "无法连接分析服务");
    }
  };

  const plotKey = embeddingMode === "tsne"
    ? viewMode === "celltype" ? "tsneCelltype" : viewMode === "sample" ? "tsneSample" : "tsne"
    : viewMode === "celltype" ? "umapCelltype" : viewMode === "sample" ? "umapSample" : "umap";
  const plotEnabled = selectedOutputs.includes(embeddingMode);
  const plotUrl = plotEnabled ? analysisResult?.artifacts?.[plotKey] : undefined;
  const displayedMarkerRows: MarkerDisplayRow[] = analysisResult?.annotationRows?.length
    ? analysisResult.annotationRows.map((row) => ({ cluster: row.cluster, label: manualLabels[row.cluster] || row.celltype_major, markers: row.markers, score: row.score.toFixed(2), confidence: row.confidence }))
    : markerRows.map(([label, markers, score, confidence]) => ({ cluster: "", label, markers, score, confidence }));
  const filteredMarkerLibrary = markerLibrary.filter(([name, genes]) => `${name} ${genes}`.toLowerCase().includes(markerQuery.toLowerCase()));
  const toggleOutput = (output: OutputFigure) => setSelectedOutputs((current) => current.includes(output) ? current.filter((item) => item !== output) : [...current, output]);
  const updateManualLabel = (cluster: string, label: string) => setManualLabels((current) => {
    const next = { ...current, [cluster]: label };
    window.localStorage.setItem("scrnalysis.manualLabels", JSON.stringify(next));
    return next;
  });
  const selectPanel = (panel: SidePanel) => setSidePanel(panel);
  const addProject = () => {
    const name = window.prompt("请输入新项目名称", "新单细胞项目");
    if (!name?.trim()) return;
    const next = Array.from(new Set([name.trim(), ...projects]));
    setProjects(next);
    setProjectName(name.trim());
    window.localStorage.setItem("scrnalysis.projects", JSON.stringify(next));
  };

  return (
    <main className="app-shell">
      <aside className="sidebar">
        <div className="brand">
          <div className="brand-mark"><span /><span /><span /></div>
          <div>
            <div className="brand-name">Scrnalysis</div>
            <div className="brand-caption">single-cell workspace</div>
          </div>
        </div>

        <div className="workspace-label">我的工作区</div>
        <nav className="side-nav" aria-label="主导航">
          <button className={`nav-item ${sidePanel === "workspace" ? "active" : ""}`} onClick={() => selectPanel("workspace")}><span className="nav-icon">⌘</span>分析工作台</button>
          <button className={`nav-item ${sidePanel === "runs" ? "active" : ""}`} onClick={() => selectPanel("runs")}><span className="nav-icon">◷</span>运行记录 <span className="nav-count">{runHistory.length}</span></button>
          <button className={`nav-item ${sidePanel === "markers" ? "active" : ""}`} onClick={() => selectPanel("markers")}><span className="nav-icon">⌁</span>Marker 库</button>
          <button className={`nav-item ${sidePanel === "files" ? "active" : ""}`} onClick={() => selectPanel("files")}><span className="nav-icon">▣</span>数据文件</button>
        </nav>

        <div className="sidebar-divider" />
        <div className="workspace-label project-label">最近项目 <button className="tiny-add" onClick={addProject} aria-label="新建项目">+</button></div>
        <div className="project-list">
          {projects.map((project, index) => <button className={`project-item ${projectName === project ? "selected" : ""}`} key={project} onClick={() => setProjectName(project)}><span className={`project-dot ${index % 3 === 0 ? "green" : index % 3 === 1 ? "blue" : "orange"}`} /><span><strong>{project}</strong><small>{index === 0 ? "刚刚更新" : index === 1 ? "昨天" : "最近使用"}</small></span></button>)}
        </div>

        <div className="sidebar-bottom">
          <div className="tip-card"><span className="tip-spark">✦</span><div><strong>AI 分析助手</strong><p>从数据摘要开始，逐步给出可解释的建议。</p></div></div>
          <button className="profile"><span className="avatar">L</span><span><strong>我的账户</strong><small>本地工作区</small></span><span className="profile-more">···</span></button>
        </div>
      </aside>

      <section className="main-content">
        <header className="topbar">
          <div className="breadcrumbs"><span>工作区</span><b>/</b><strong>{projectName}</strong><span className="status-pill"><i /> {runStatus === "complete" ? "已完成" : "草稿"}</span></div>
          <div className="top-actions"><button className="icon-button" aria-label="帮助" onClick={() => window.alert("先上传 10x 三文件或 zip，再填写研究背景、研究目的和比较分组，然后点击开始分析。")}>?</button><button className="icon-button" aria-label="通知" onClick={() => selectPanel("runs")}>♧</button><button className="outline-button" onClick={() => window.localStorage.setItem("scrnalysis.projectDraft", JSON.stringify({ projectName, selectedFiles: selectedFiles.map((file) => file.name) }))}>保存项目</button><button className="primary-button compact" onClick={startAnalysis}>运行分析 <span>→</span></button></div>
        </header>

        {sidePanel !== "workspace" && <aside className="workspace-drawer" aria-label="工作区面板">
          <div className="drawer-head"><div><span className="section-kicker">WORKSPACE</span><h2>{sidePanel === "runs" ? "运行记录" : sidePanel === "markers" ? "Marker 库" : "数据文件"}</h2></div><button className="icon-button" onClick={() => selectPanel("workspace")} aria-label="关闭">×</button></div>
          {sidePanel === "runs" && <div className="drawer-list">{runHistory.length ? runHistory.map((run) => <div className="drawer-item" key={run.job_id}><strong>{run.job_id}</strong><span>{run.status === "complete" ? "已完成" : run.status}</span><small>{run.files.join("、") || "未记录文件名"}</small>{run.status === "complete" && analysisResult?.job_id === run.job_id && <a href={`${API_BASE}${analysisResult.artifacts?.report || ""}`} target="_blank" rel="noreferrer">查看本次报告 →</a>}</div>) : <div className="drawer-empty">还没有完成的分析。完成一次真实分析后，记录会自动保存在这台浏览器中。</div>}</div>}
          {sidePanel === "markers" && <div className="drawer-list"><input className="drawer-search" value={markerQuery} onChange={(event) => setMarkerQuery(event.target.value)} placeholder="搜索细胞类型或 marker" />{filteredMarkerLibrary.map(([name, genes]) => <div className="drawer-item" key={name}><strong>{name}</strong><small>{genes}</small></div>)}</div>}
          {sidePanel === "files" && <div className="drawer-list">{selectedFiles.length > 0 && <div className="drawer-note">本次已选择 {selectedFiles.length} 个文件。</div>}{savedFileNames.length ? savedFileNames.map((name) => <div className="drawer-item" key={name}><strong>{name}</strong><small>浏览器记录的文件名；原始文件不会上传到云端列表。</small></div>) : <div className="drawer-empty">选择数据文件后，这里会显示最近使用过的文件名。</div>}</div>}
        </aside>}

        <div className="content-wrap">
          <div className="page-heading">
            <div><div className="eyebrow">ANALYSIS WORKSPACE <span>·</span> RUN 004</div><h1>让数据自己说话。</h1><p>上传你的单细胞数据，Scrnalysis 会把复杂流程拆成清晰、可检查的步骤。</p></div>
            <div className="heading-note"><span className="ai-orb">✦</span><span><strong>AI 建议已开启</strong><small>只会建议参数，不会替你隐藏证据</small></span></div>
          </div>

          <div className="run-progress">
            {pipeline.map(([number, name, detail], index) => <div className={`progress-step ${index < activeStep || runStatus === "complete" ? "done" : ""} ${index === activeStep && runStatus !== "complete" ? "current" : ""}`} key={number}><span className="step-number">{index < activeStep || runStatus === "complete" ? "✓" : number}</span><span><strong>{name}</strong><small>{detail}</small></span></div>)}
          </div>

          <div className="workspace-grid">
            <div className="left-column">
              <section className="card upload-card">
                <div className="card-title-row"><div><span className="section-kicker">STEP 01</span><h2>上传数据 <em className="required-mark" aria-hidden="true">*</em></h2></div><span className="soft-tag">支持 10x Matrix Market</span></div>
                <div className={`dropzone ${selectedFiles.length ? "has-file" : ""} ${validationErrors.files ? "invalid" : ""}`} onClick={() => fileInput.current?.click()} onDragOver={(event) => event.preventDefault()} onDrop={onDrop} role="button" tabIndex={0} aria-invalid={Boolean(validationErrors.files)} onKeyDown={(event) => event.key === "Enter" && fileInput.current?.click()}>
                  <input ref={fileInput} type="file" hidden multiple accept=".zip,.gz,.mtx,.tsv" onChange={onFileChange} />
                  <div className="upload-icon">↥</div>
                  {selectedFiles.length ? <><strong>{selectedFiles.length === 1 ? selectedFiles[0].name : `已选择 ${selectedFiles.length} 个文件`}</strong><span>{(selectedFiles.reduce((total, file) => total + file.size, 0) / 1024 / 1024).toFixed(2)} MB · 已准备分析</span><button className="text-button" onClick={(event) => { event.stopPropagation(); setSelectedFiles([]); setAnalysisResult(null); setRunStatus("ready"); }}>移除文件</button></> : <><strong>拖放 10x 文件到这里，或点击上传</strong><span>支持 zip 压缩包，或同时选择 matrix / barcodes / features 三个文件</span><button className="upload-button" onClick={(event) => { event.stopPropagation(); fileInput.current?.click(); }}>选择文件</button></>}
                </div>
                {validationErrors.files && <div className="field-error upload-error" role="alert">{validationErrors.files}</div>}
                <div className="upload-footnote"><span className="lock">⌑</span> 文件会发送到运行 Scrnalysis 的这台电脑，只用于本次分析。</div>
              </section>

              <section className="card context-card">
                <div className="card-title-row"><div><span className="section-kicker">STEP 02 · PROJECT CONTEXT</span><h2>分析背景与目的</h2></div><span className="soft-tag">帮助 AI 理解方向</span></div>
                <label className={`context-field ${validationErrors.background ? "has-error" : ""}`}><span>研究背景 <em className="required-mark" aria-hidden="true">*</em></span><textarea value={analysisBackground} aria-invalid={Boolean(validationErrors.background)} aria-describedby={validationErrors.background ? "background-error" : undefined} onChange={(event) => { setAnalysisBackground(event.target.value); clearValidationError("background"); }} placeholder="例如：本项目来自类风湿关节炎滑膜组织，关注免疫细胞与基质细胞的组成变化。" rows={3} />{validationErrors.background && <small id="background-error" className="field-error" role="alert">{validationErrors.background}</small>}</label>
                <label className={`context-field ${validationErrors.purpose ? "has-error" : ""}`}><span>研究目的 <em className="required-mark" aria-hidden="true">*</em></span><textarea value={researchPurpose} aria-invalid={Boolean(validationErrors.purpose)} aria-describedby={validationErrors.purpose ? "purpose-error" : undefined} onChange={(event) => { setResearchPurpose(event.target.value); clearValidationError("purpose"); }} placeholder="例如：比较疾病组与对照组的细胞比例、亚群状态和差异表达。" rows={3} />{validationErrors.purpose && <small id="purpose-error" className="field-error" role="alert">{validationErrors.purpose}</small>}</label>
                <div className="context-inline"><label className={`context-field ${validationErrors.groups ? "has-error" : ""}`}><span>比较分组 <em className="required-mark" aria-hidden="true">*</em></span><input value={comparisonGroups} aria-invalid={Boolean(validationErrors.groups)} aria-describedby={validationErrors.groups ? "groups-error" : undefined} onChange={(event) => { setComparisonGroups(event.target.value); clearValidationError("groups"); }} placeholder="例如：RA vs Control；无对照可填：单组探索" />{validationErrors.groups && <small id="groups-error" className="field-error" role="alert">{validationErrors.groups}</small>}</label><label className={`context-field ${validationErrors.annotation ? "has-error" : ""}`}><span>注释粒度 <em className="required-mark" aria-hidden="true">*</em></span><select value={annotationScope} aria-invalid={Boolean(validationErrors.annotation)} onChange={(event) => { setAnnotationScope(event.target.value); clearValidationError("annotation"); }}><option value="major_and_detail">主要类型 + 细分亚群</option><option value="major">仅主要细胞类型</option><option value="detail">尽可能细的亚群</option></select>{validationErrors.annotation && <small className="field-error" role="alert">{validationErrors.annotation}</small>}</label></div>
                <div className="context-note"><span>✦</span> 这些信息会写入本次分析的运行记录，帮助 AI 选择 marker、解释异常群体，并避免脱离研究问题自动改参数。</div>
              </section>

              <section className="card parameters-card">
                <div className="card-title-row"><div><span className="section-kicker">STEP 03</span><h2>分析参数</h2></div><button className={`advanced-toggle ${showAdvanced ? "open" : ""}`} onClick={() => setShowAdvanced(!showAdvanced)}>高级参数 <span>⌄</span></button></div>
                <div className="parameter-grid">
                  <label className={validationErrors.species ? "has-error" : ""}><span>物种 <em className="required-mark" aria-hidden="true">*</em></span><select value={species} onChange={(event) => { setSpecies(event.target.value); clearValidationError("species"); }}><option value="human">Human / 人</option><option value="mouse">Mouse / 小鼠</option></select>{validationErrors.species && <small className="field-error" role="alert">{validationErrors.species}</small>}</label>
                  <label className={validationErrors.tissue ? "has-error" : ""}><span>组织 <em className="required-mark" aria-hidden="true">*</em></span><select value={tissue} onChange={(event) => { setTissue(event.target.value); clearValidationError("tissue"); }}><option value="synovium">Synovium / 滑膜</option><option value="pbmc">PBMC / 外周血</option><option value="blood">Whole blood / 全血</option><option value="lung">Lung / 肺</option><option value="skin">Skin / 皮肤</option><option value="kidney">Kidney / 肾脏</option><option value="muscle">Muscle / 肌肉</option><option value="salivary_gland">Salivary gland / 唾液腺</option><option value="bone_marrow">Bone marrow / 骨髓</option><option value="lymph_node">Lymph node / 淋巴结</option><option value="bal">BAL / 支气管肺泡灌洗液</option></select>{validationErrors.tissue && <small className="field-error" role="alert">{validationErrors.tissue}</small>}</label>
                  <label className={validationErrors.normalization ? "has-error" : ""}><span>标准化方法 <em className="required-mark" aria-hidden="true">*</em></span><select value={normalization} onChange={(event) => { setNormalization(event.target.value); clearValidationError("normalization"); }}><option value="log">LogNormalize</option><option value="sct">SCTransform</option></select>{validationErrors.normalization && <small className="field-error" role="alert">{validationErrors.normalization}</small>}</label>
                  <label className={validationErrors.batch ? "has-error" : ""}><span>批次校正 <em className="required-mark" aria-hidden="true">*</em></span><select value={batchCorrection} onChange={(event) => { setBatchCorrection(event.target.value); clearValidationError("batch"); }}><option value="harmony">Harmony</option><option value="anchors">Seurat anchors</option><option value="none">不进行整合</option></select>{validationErrors.batch && <small className="field-error" role="alert">{validationErrors.batch}</small>}</label>
                </div>
                {showAdvanced && <div className="advanced-panel">
                  <div className="advanced-section full-span qc-section-new"><div className="advanced-heading"><strong>QC 与质量控制</strong><span>可选固定阈值、上下分位数，或对单项不设置</span></div><div className="qc-rule-grid">
                    <div className="qc-rule-card"><div className="qc-rule-head"><strong>检测基因数</strong><span>默认去除上下各 2.5%</span></div><label className="qc-mode-field"><span>筛选方式</span><select value={geneFilterMode} onChange={(event) => setGeneFilterMode(event.target.value)}><option value="quantile">按分位数</option><option value="fixed">固定阈值</option><option value="none">不设置</option></select></label>{geneFilterMode === "fixed" && <div className="qc-input-grid"><label><span>最少基因数</span><div className="number-row"><input type="number" min="0" step="50" value={minGenes} onChange={(event) => setMinGenes(Number(event.target.value))} /><small>genes</small></div></label><label><span>最多基因数</span><div className="number-row"><input type="number" min="0" step="50" value={maxGenes} onChange={(event) => setMaxGenes(Number(event.target.value))} /><small>genes</small></div></label></div>}{geneFilterMode === "quantile" && <div className="qc-input-grid"><label><span>去除下方</span><div className="number-row"><input type="number" min="0" max="49" step="0.5" value={geneLowerQuantile * 100} onChange={(event) => setGeneLowerQuantile(Number(event.target.value) / 100)} /><small>%</small></div></label><label><span>去除上方</span><div className="number-row"><input type="number" min="0" max="49" step="0.5" value={geneUpperQuantile * 100} onChange={(event) => setGeneUpperQuantile(Number(event.target.value) / 100)} /><small>%</small></div></label></div>}{geneFilterMode === "none" && <div className="qc-disabled-note">本项不限制基因数。</div>}</div>
                    <div className="qc-rule-card"><div className="qc-rule-head"><strong>UMI / nCount</strong><span>默认去除上下各 2.5%</span></div><label className="qc-mode-field"><span>筛选方式</span><select value={umiFilterMode} onChange={(event) => setUmiFilterMode(event.target.value)}><option value="quantile">按分位数</option><option value="fixed">固定阈值</option><option value="none">不设置</option></select></label>{umiFilterMode === "fixed" && <div className="qc-input-grid"><label><span>最少 UMI</span><div className="number-row"><input type="number" min="0" step="500" value={minCounts} onChange={(event) => setMinCounts(Number(event.target.value))} /><small>UMIs</small></div></label><label><span>最多 UMI</span><div className="number-row"><input type="number" min="0" step="500" value={maxCounts} onChange={(event) => setMaxCounts(Number(event.target.value))} /><small>UMIs</small></div></label></div>}{umiFilterMode === "quantile" && <div className="qc-input-grid"><label><span>去除下方</span><div className="number-row"><input type="number" min="0" max="49" step="0.5" value={umiLowerQuantile * 100} onChange={(event) => setUmiLowerQuantile(Number(event.target.value) / 100)} /><small>%</small></div></label><label><span>去除上方</span><div className="number-row"><input type="number" min="0" max="49" step="0.5" value={umiUpperQuantile * 100} onChange={(event) => setUmiUpperQuantile(Number(event.target.value) / 100)} /><small>%</small></div></label></div>}{umiFilterMode === "none" && <div className="qc-disabled-note">本项不限制 UMI 数。</div>}</div>
                    <div className="qc-rule-card"><div className="qc-rule-head"><strong>线粒体比例</strong><span>常用上限 20%</span></div><label className="qc-mode-field"><span>筛选方式</span><select value={mitoFilterMode} onChange={(event) => setMitoFilterMode(event.target.value)}><option value="fixed">固定上限</option><option value="none">不设置</option></select></label>{mitoFilterMode === "fixed" ? <label><span>线粒体比例上限</span><div className="number-row"><input type="number" min="0" max="100" step="1" value={maxMito} onChange={(event) => setMaxMito(Number(event.target.value))} /><small>%</small></div></label> : <div className="qc-disabled-note">本项不限制线粒体比例。</div>}</div>
                    <div className="qc-rule-card"><div className="qc-rule-head"><strong>核糖体比例</strong><span>默认不限制</span></div><label className="qc-mode-field"><span>筛选方式</span><select value={riboFilterMode} onChange={(event) => setRiboFilterMode(event.target.value)}><option value="none">不设置</option><option value="fixed">固定上限</option></select></label>{riboFilterMode === "fixed" ? <label><span>核糖体比例上限</span><div className="number-row"><input type="number" min="0" max="100" step="1" value={maxRibo} onChange={(event) => setMaxRibo(Number(event.target.value))} /><small>%</small></div></label> : <div className="qc-disabled-note">本项不限制核糖体比例。</div>}</div>
                    <div className="qc-rule-card"><div className="qc-rule-head"><strong>血红蛋白 / Hb</strong><span>血液样本可按需限制</span></div><label className="qc-mode-field"><span>筛选方式</span><select value={hbFilterMode} onChange={(event) => setHbFilterMode(event.target.value)}><option value="none">不设置</option><option value="fixed">固定上限</option></select></label>{hbFilterMode === "fixed" ? <label><span>Hb 比例上限</span><div className="number-row"><input type="number" min="0" max="100" step="1" value={maxHb} onChange={(event) => setMaxHb(Number(event.target.value))} /><small>%</small></div></label> : <div className="qc-disabled-note">本项不限制 Hb 比例。</div>}</div>
                    <div className="qc-rule-card"><div className="qc-rule-head"><strong>基因最低出现细胞数</strong><span>默认至少 3 个细胞</span></div><label className="qc-mode-field"><span>筛选方式</span><select value={minCellsPerGeneMode} onChange={(event) => setMinCellsPerGeneMode(event.target.value)}><option value="fixed">固定阈值</option><option value="none">不设置</option></select></label>{minCellsPerGeneMode === "fixed" ? <label><span>最低出现细胞数</span><div className="number-row"><input type="number" min="0" max="100" step="1" value={minCellsPerGene} onChange={(event) => setMinCellsPerGene(Number(event.target.value))} /><small>cells</small></div></label> : <div className="qc-disabled-note">保留所有至少出现过 1 次的基因。</div>}</div>
                  </div></div>
                  <div className="advanced-section full-span"><div className="advanced-heading"><strong>QC 与质量控制</strong><span>按样本记录阈值与过滤前后细胞数</span></div><div className="advanced-controls"><label><span>线粒体比例上限</span><div className="number-row"><input type="number" min="5" max="40" step="1" value={maxMito} onChange={(event) => setMaxMito(Number(event.target.value))} /><small>%</small></div></label><label><span>最少基因数</span><div className="number-row"><input type="number" min="50" max="1000" step="50" value={minGenes} onChange={(event) => setMinGenes(Number(event.target.value))} /><small>genes</small></div></label><label><span>最多 UMI 数</span><div className="number-row"><input type="number" min="5000" max="100000" step="500" value={maxCounts} onChange={(event) => setMaxCounts(Number(event.target.value))} /><small>UMIs</small></div></label><label><span>核糖体比例上限</span><div className="number-row"><input type="number" min="5" max="60" step="1" value={maxRibo} onChange={(event) => setMaxRibo(Number(event.target.value))} /><small>%</small></div></label><label><span>基因最低出现细胞数</span><div className="number-row"><input type="number" min="1" max="10" step="1" value={minCellsPerGene} onChange={(event) => setMinCellsPerGene(Number(event.target.value))} /><small>cells</small></div></label></div></div>
                  <div className="advanced-section full-span"><div className="advanced-heading"><strong>双细胞与环境 RNA</strong><span>默认保留原始对象，并输出处理前后对照</span></div><div className="advanced-selects"><label><span>双细胞检测</span><select value={doubletMethod} onChange={(event) => setDoubletMethod(event.target.value)}><option value="scDblFinder">scDblFinder</option><option value="doubletFinder">DoubletFinder</option><option value="none">不检测</option></select></label><label><span>环境 RNA 诊断</span><select value={ambientMethod} onChange={(event) => setAmbientMethod(event.target.value)}><option value="decontX">decontX</option><option value="soupX">SoupX</option><option value="diagnostic_only">仅诊断，不校正</option></select></label></div></div>
                  <div className="advanced-section full-span"><div className="advanced-heading"><strong>降维、聚类与图形</strong><span>这些参数会进入可复现运行记录</span></div><div className="advanced-controls"><label><span>PCA 维度</span><div className="number-row"><input type="number" min="10" max="50" step="1" value={dims} onChange={(event) => setDims(Number(event.target.value))} /><small>dims</small></div></label><label><span>聚类分辨率</span><div className="number-row"><input type="number" min="0.1" max="1.2" step="0.1" value={resolution} onChange={(event) => setResolution(Number(event.target.value))} /><small>resolution</small></div></label><label><span>邻居数 k.param</span><div className="number-row"><input type="number" min="10" max="100" step="1" value={neighbors} onChange={(event) => setNeighbors(Number(event.target.value))} /><small>neighbors</small></div></label><label><span>UMAP 邻居数</span><div className="number-row"><input type="number" min="5" max="100" step="1" value={umapNeighbors} onChange={(event) => setUmapNeighbors(Number(event.target.value))} /><small>neighbors</small></div></label><label><span>UMAP min.dist</span><div className="number-row"><input type="number" min="0.05" max="0.9" step="0.05" value={umapMinDist} onChange={(event) => setUmapMinDist(Number(event.target.value))} /><small>min.dist</small></div></label><label><span>t-SNE perplexity</span><div className="number-row"><input type="number" min="5" max="100" step="1" value={tsnePerplexity} onChange={(event) => setTsnePerplexity(Number(event.target.value))} /><small>perplexity</small></div></label></div></div>
                </div>}
                <div className="ai-suggestion"><span className="suggestion-icon">✦</span><span><strong>AI 建议</strong> 当前参数适合滑膜免疫细胞数据。运行前会根据每个样本的 QC 分布再次检查阈值。</span><button>查看依据 →</button></div>
              </section>

              <section className="card annotation-card">
                <div className="card-title-row"><div><span className="section-kicker">STEP 04</span><h2>注释策略</h2></div><span className="ready-tag"><i /> 可开始</span></div>
                <div className="annotation-choice"><div className="choice-icon">✦</div><div><strong>AI 辅助注释</strong><p>运行时根据 cluster markers、经典 marker 和组织背景生成初步细胞类型建议。</p></div><button className={`switch ${aiAnnotationEnabled ? "on" : ""}`} aria-pressed={aiAnnotationEnabled} onClick={() => setAiAnnotationEnabled((value) => !value)}><i /></button></div>
                <div className="annotation-choice"><div className="choice-icon manual">⌁</div><div><strong>保留手动修改</strong><p>保留 cluster 到细胞类型的修改入口，并在结果证据中记录原始建议。</p></div><button className={`switch ${manualEditsEnabled ? "on" : ""}`} aria-pressed={manualEditsEnabled} onClick={() => setManualEditsEnabled((value) => !value)}><i /></button></div>
                <div className="annotation-status"><span>运行后会输出 cluster 注释建议、marker 证据和置信度。</span><strong>{aiAnnotationEnabled ? "AI 已开启" : "仅输出未注释 cluster"}</strong></div>
              </section>

              <section className={`card output-card ${validationErrors.outputs ? "has-output-error" : ""}`}>
                <div className="card-title-row"><div><span className="section-kicker">STEP 05</span><h2>输出文件</h2></div><span className="soft-tag">可在运行前选择</span></div>
                <p className="output-intro">勾选后，分析脚本会生成你选择的基础图形；坐标表、分析报告和 QC 汇总表会始终保留，便于复核。</p>
                <div className="output-grid">{([ ["umap", "UMAP", "按细胞类型、样本和 cluster 输出"], ["tsne", "t-SNE", "输出 t-SNE 坐标和分组图"], ["dotplot", "DotPlot", "主要 marker 的表达比例和平均表达"], ["featureplot", "FeaturePlot", "经典 marker 的空间表达"], ["qc_violin", "QC violin plot", "nFeature、nCount、线粒体比例"], ] as const).map(([value, label, detail]) => <label className="output-option" key={value}><input type="checkbox" checked={selectedOutputs.includes(value)} onChange={() => toggleOutput(value)} /><span><strong>{label}</strong><small>{detail}</small></span></label>)}</div>
                {validationErrors.outputs && <small className="field-error" role="alert">{validationErrors.outputs}</small>}
              </section>
            </div>

            <div className="right-column">
              <section className="card results-card">
                <div className="card-title-row"><div><span className="section-kicker">LIVE PREVIEW</span><h2>结果预览</h2></div><span className={`run-tag ${runStatus}`}>{runStatus === "running" ? "分析中 · " + (activeStep + 1) + "/5" : runStatus === "complete" ? "已完成" : runStatus === "error" ? "需要处理" : "等待上传"}</span></div>
                <div className="plot-toolbar"><div className="segmented"><button className={embeddingMode === "umap" ? "selected" : ""} onClick={() => setEmbeddingMode("umap")}>UMAP</button><button className={embeddingMode === "tsne" ? "selected" : ""} onClick={() => setEmbeddingMode("tsne")}>t-SNE</button></div><div className="segmented"><button className={viewMode === "celltype" ? "selected" : ""} onClick={() => setViewMode("celltype")}>细胞类型</button><button className={viewMode === "sample" ? "selected" : ""} onClick={() => setViewMode("sample")}>样本</button><button className={viewMode === "cluster" ? "selected" : ""} onClick={() => setViewMode("cluster")}>Cluster</button></div><button className="plot-action" onClick={() => plotUrl && window.open(`${API_BASE}${plotUrl}`, "_blank", "noopener,noreferrer")}>⤢ 放大</button></div>
                <div className={`demo-note ${analysisResult?.status === "complete" ? "real-note" : ""}`}><span>{analysisResult?.status === "complete" ? "✓" : "ⓘ"}</span> {analysisResult?.status === "complete" ? `这是本次上传数据生成的真实 ${embeddingMode.toUpperCase()}，当前按${viewMode === "celltype" ? "细胞类型" : viewMode === "sample" ? "样本" : "cluster"}着色。` : "上传数据后，这里会显示真实的 UMAP 或 t-SNE；未上传时显示界面示例。"}</div>
                {plotUrl ? <div className="real-plot-wrap"><img className="real-plot" src={`${API_BASE}${plotUrl}`} alt={`本次单细胞分析生成的 ${embeddingMode} 图`} /></div> : <div className="plot-area"><div className="axis-label y">{embeddingMode === "umap" ? "UMAP_2" : "t-SNE_2"}</div><div className="axis-label x">{embeddingMode === "umap" ? "UMAP_1" : "t-SNE_1"}</div>{points.map((point, index) => <span key={index} className="plot-point" style={{ left: `${point.left}%`, top: `${point.top}%`, background: point.color, width: point.size, height: point.size }} />)}</div>}
                <div className="plot-footer"><div className="legend">{analysisResult?.status === "complete" ? <span><i style={{ background: "#157c72" }} />真实 {embeddingMode.toUpperCase()}</span> : cellTypes.map((type) => <span key={type.name}><i style={{ background: type.color }} />{type.name}</span>)}</div><span className="plot-meta">{analysisResult?.status === "complete" ? <a href={`${API_BASE}${analysisResult.artifacts?.[embeddingMode === "umap" ? "umapCoordinates" : "tsneCoordinates"] || ""}`} target="_blank" rel="noreferrer">下载 {embeddingMode.toUpperCase()} 坐标</a> : "等待真实数据"}</span></div>
                {analysisResult?.status === "complete" && <div className="result-downloads"><span>本次已生成：</span>{selectedOutputs.map((output) => { const artifact = output === "umap" ? "umap" : output === "tsne" ? "tsne" : output === "qc_violin" ? "qcViolinAfter" : output; const label = output === "qc_violin" ? "QC violin" : output === "featureplot" ? "FeaturePlot" : output === "dotplot" ? "DotPlot" : output.toUpperCase(); return <a key={output} href={`${API_BASE}${analysisResult.artifacts?.[artifact] || ""}`} target="_blank" rel="noreferrer">{label}</a>; })}<a href={`${API_BASE}${analysisResult.artifacts?.report || ""}`} target="_blank" rel="noreferrer">分析报告</a></div>}
                {errorMessage && <div className="service-error">{errorMessage}</div>}
              </section>

              <section className="card insight-card"><div className="insight-head"><span className="ai-orb small">✦</span><div><span className="section-kicker">AI INSIGHT</span><h2>初步观察</h2></div><button className="more-button" onClick={() => selectPanel("runs")}>···</button></div><p>{analysisResult?.insight?.summary || "上传并完成分析后，这里会根据保留细胞数、cluster marker 和注释置信度生成本次运行的初步观察。"}</p><div className="insight-bottom"><span><i className="confidence-dot" /> 置信度 {analysisResult?.insight?.confidence ? `${Math.round(analysisResult.insight.confidence * 100)}%` : "待分析"}</span><button onClick={() => selectPanel("markers")}>查看 marker 证据 →</button></div></section>

              <section className="card markers-card"><div className="card-title-row"><div><span className="section-kicker">ANNOTATION EVIDENCE</span><h2>注释证据</h2></div><button className="text-link" onClick={() => selectPanel("markers")}>全部查看 →</button></div><div className="marker-table"><div className="marker-row header"><span>细胞类型</span><span>代表性 markers</span><span>评分</span><span /></div>{displayedMarkerRows.map((row, index) => <div className="marker-row" key={`${row.cluster || row.label}-${index}`}><span className="marker-name"><i style={{ background: cellTypes[index % cellTypes.length].color }} />{row.label}</span><span className="marker-genes">{row.markers}</span><span className="score">{row.score}</span><span className={`confidence ${row.confidence === "高" ? "high" : "medium"}`}>{row.confidence}</span></div>)}</div>{manualEditsEnabled && analysisResult?.annotationRows?.length ? <div className="manual-labels"><strong>手动修改 cluster 标签</strong><span>只修改网页显示标签；原始 AI 建议仍保留在证据表和分析报告中。</span>{analysisResult.annotationRows.map((row) => <label key={row.cluster}><span>Cluster {row.cluster}</span><input value={manualLabels[row.cluster] || row.celltype_major} onChange={(event) => updateManualLabel(row.cluster, event.target.value)} /></label>)}</div> : null}</section>
            </div>
          </div>

          <section className="run-banner"><div className="run-banner-icon">{runStatus === "running" ? "◌" : runStatus === "complete" ? "✓" : runStatus === "error" ? "!" : "✦"}</div><div><strong>{runStatus === "running" ? `正在执行：${analysisResult?.message || pipeline[activeStep][1]}` : runStatus === "complete" ? "本次真实分析已完成" : runStatus === "error" ? "分析没有完成" : "准备好开始你的第一次分析了吗？"}</strong><p>{runStatus === "running" ? `${analysisResult?.progress || 0}% · ${analysisResult?.message || pipeline[activeStep][2]}` : runStatus === "complete" ? "你可以查看真实 UMAP、QC、marker 和结果报告。" : runStatus === "error" ? errorMessage : "请先上传 10x 三文件或 zip 压缩包，再开始分析。"}</p></div><button className="primary-button" onClick={() => void startAnalysis()} disabled={runStatus === "running"}>{runStatus === "running" ? "分析进行中…" : runStatus === "complete" ? "重新运行" : "开始分析"}<span>→</span></button></section>
          <footer className="footer"><span>Scrnalysis · 为自己的数据保留完整证据链</span><span>v0.2 <i /> 本地分析服务已连接</span></footer>
        </div>
      </section>
    </main>
  );
}
