"use client";

import { ChangeEvent, DragEvent, useMemo, useRef, useState } from "react";

type RunStatus = "ready" | "running" | "complete" | "error";
type ViewMode = "celltype" | "sample" | "cluster";

type AnalysisResult = {
  job_id: string;
  status: RunStatus;
  progress: number;
  message: string;
  error?: string;
  artifacts?: Record<string, string>;
};

type ValidationErrors = Partial<Record<"files" | "background" | "purpose" | "groups" | "species" | "tissue" | "normalization" | "batch" | "annotation", string>>;

const API_BASE = process.env.NEXT_PUBLIC_SCRNALYSIS_API || (typeof window !== "undefined" ? `http://${window.location.hostname}:8000` : "http://localhost:8000");

const pipeline = [
  ["01", "上传数据", "读取矩阵与样本信息"],
  ["02", "背景与目的", "明确研究问题与比较分组"],
  ["03", "分析参数", "QC / 降维 / 整合参数"],
  ["04", "细胞注释", "marker 证据与 AI 建议"],
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
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [maxMito, setMaxMito] = useState(20);
  const [minGenes, setMinGenes] = useState(200);
  const [dims, setDims] = useState(30);
  const [resolution, setResolution] = useState(0.6);
  const [maxCounts, setMaxCounts] = useState(25000);
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
  const [annotationScope, setAnnotationScope] = useState("major_and_detail");
  const [validationErrors, setValidationErrors] = useState<ValidationErrors>({});
  const points = useMemo(() => createPoints(), []);

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
    setActiveStep(Math.min(pipeline.length - 1, Math.floor((job.progress || 0) / 25)));
    if (job.status === "complete") {
      setRunStatus("complete");
      setActiveStep(pipeline.length - 1);
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
    formData.append("config", JSON.stringify({ species, tissue, normalization, batchCorrection, analysisBackground, researchPurpose, comparisonGroups, annotationScope, maxMito, minGenes, maxCounts, minCellsPerGene, maxRibo, doubletMethod, ambientMethod, dims, resolution, neighbors, umapNeighbors, umapMinDist, tsnePerplexity }));
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
          <button className="nav-item active"><span className="nav-icon">⌘</span>分析工作台</button>
          <button className="nav-item"><span className="nav-icon">◷</span>运行记录 <span className="nav-count">3</span></button>
          <button className="nav-item"><span className="nav-icon">⌁</span>Marker 库</button>
          <button className="nav-item"><span className="nav-icon">▣</span>数据文件</button>
        </nav>

        <div className="sidebar-divider" />
        <div className="workspace-label project-label">最近项目 <button className="tiny-add">+</button></div>
        <div className="project-list">
          <button className="project-item selected"><span className="project-dot green" /><span><strong>RA_synovium_demo</strong><small>刚刚更新</small></span></button>
          <button className="project-item"><span className="project-dot blue" /><span><strong>SLE_PBMC_01</strong><small>昨天</small></span></button>
          <button className="project-item"><span className="project-dot orange" /><span><strong>scRNA_exploration</strong><small>8 月 12 日</small></span></button>
        </div>

        <div className="sidebar-bottom">
          <div className="tip-card"><span className="tip-spark">✦</span><div><strong>AI 分析助手</strong><p>从数据摘要开始，逐步给出可解释的建议。</p></div></div>
          <button className="profile"><span className="avatar">L</span><span><strong>我的账户</strong><small>本地工作区</small></span><span className="profile-more">···</span></button>
        </div>
      </aside>

      <section className="main-content">
        <header className="topbar">
          <div className="breadcrumbs"><span>工作区</span><b>/</b><strong>RA_synovium_demo</strong><span className="status-pill"><i /> 草稿</span></div>
          <div className="top-actions"><button className="icon-button" aria-label="帮助">?</button><button className="icon-button" aria-label="通知">♧</button><button className="outline-button">保存项目</button><button className="primary-button compact" onClick={startAnalysis}>运行分析 <span>→</span></button></div>
        </header>

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
                  <label className={validationErrors.tissue ? "has-error" : ""}><span>组织 <em className="required-mark" aria-hidden="true">*</em></span><select value={tissue} onChange={(event) => { setTissue(event.target.value); clearValidationError("tissue"); }}><option value="synovium">Synovium / 滑膜</option><option value="pbmc">PBMC</option><option value="lung">Lung / 肺</option></select>{validationErrors.tissue && <small className="field-error" role="alert">{validationErrors.tissue}</small>}</label>
                  <label className={validationErrors.normalization ? "has-error" : ""}><span>标准化方法 <em className="required-mark" aria-hidden="true">*</em></span><select value={normalization} onChange={(event) => { setNormalization(event.target.value); clearValidationError("normalization"); }}><option value="log">LogNormalize</option><option value="sct">SCTransform</option></select>{validationErrors.normalization && <small className="field-error" role="alert">{validationErrors.normalization}</small>}</label>
                  <label className={validationErrors.batch ? "has-error" : ""}><span>批次校正 <em className="required-mark" aria-hidden="true">*</em></span><select value={batchCorrection} onChange={(event) => { setBatchCorrection(event.target.value); clearValidationError("batch"); }}><option value="harmony">Harmony</option><option value="anchors">Seurat anchors</option><option value="none">不进行整合</option></select>{validationErrors.batch && <small className="field-error" role="alert">{validationErrors.batch}</small>}</label>
                </div>
                {showAdvanced && <div className="advanced-panel">
                  <div className="advanced-section full-span"><div className="advanced-heading"><strong>QC 与质量控制</strong><span>按样本记录阈值与过滤前后细胞数</span></div><div className="advanced-controls"><label><span>线粒体比例上限</span><div className="number-row"><input type="number" min="5" max="40" step="1" value={maxMito} onChange={(event) => setMaxMito(Number(event.target.value))} /><small>%</small></div></label><label><span>最少基因数</span><div className="number-row"><input type="number" min="50" max="1000" step="50" value={minGenes} onChange={(event) => setMinGenes(Number(event.target.value))} /><small>genes</small></div></label><label><span>最多 UMI 数</span><div className="number-row"><input type="number" min="5000" max="100000" step="500" value={maxCounts} onChange={(event) => setMaxCounts(Number(event.target.value))} /><small>UMIs</small></div></label><label><span>核糖体比例上限</span><div className="number-row"><input type="number" min="5" max="60" step="1" value={maxRibo} onChange={(event) => setMaxRibo(Number(event.target.value))} /><small>%</small></div></label><label><span>基因最低出现细胞数</span><div className="number-row"><input type="number" min="1" max="10" step="1" value={minCellsPerGene} onChange={(event) => setMinCellsPerGene(Number(event.target.value))} /><small>cells</small></div></label></div></div>
                  <div className="advanced-section full-span"><div className="advanced-heading"><strong>双细胞与环境 RNA</strong><span>默认保留原始对象，并输出处理前后对照</span></div><div className="advanced-selects"><label><span>双细胞检测</span><select value={doubletMethod} onChange={(event) => setDoubletMethod(event.target.value)}><option value="scDblFinder">scDblFinder</option><option value="doubletFinder">DoubletFinder</option><option value="none">不检测</option></select></label><label><span>环境 RNA 诊断</span><select value={ambientMethod} onChange={(event) => setAmbientMethod(event.target.value)}><option value="decontX">decontX</option><option value="soupX">SoupX</option><option value="diagnostic_only">仅诊断，不校正</option></select></label></div></div>
                  <div className="advanced-section full-span"><div className="advanced-heading"><strong>降维、聚类与图形</strong><span>这些参数会进入可复现运行记录</span></div><div className="advanced-controls"><label><span>PCA 维度</span><div className="number-row"><input type="number" min="10" max="50" step="1" value={dims} onChange={(event) => setDims(Number(event.target.value))} /><small>dims</small></div></label><label><span>聚类分辨率</span><div className="number-row"><input type="number" min="0.1" max="1.2" step="0.1" value={resolution} onChange={(event) => setResolution(Number(event.target.value))} /><small>resolution</small></div></label><label><span>邻居数 k.param</span><div className="number-row"><input type="number" min="10" max="100" step="1" value={neighbors} onChange={(event) => setNeighbors(Number(event.target.value))} /><small>neighbors</small></div></label><label><span>UMAP 邻居数</span><div className="number-row"><input type="number" min="5" max="100" step="1" value={umapNeighbors} onChange={(event) => setUmapNeighbors(Number(event.target.value))} /><small>neighbors</small></div></label><label><span>UMAP min.dist</span><div className="number-row"><input type="number" min="0.05" max="0.9" step="0.05" value={umapMinDist} onChange={(event) => setUmapMinDist(Number(event.target.value))} /><small>min.dist</small></div></label><label><span>t-SNE perplexity</span><div className="number-row"><input type="number" min="5" max="100" step="1" value={tsnePerplexity} onChange={(event) => setTsnePerplexity(Number(event.target.value))} /><small>perplexity</small></div></label></div></div>
                </div>}
                <div className="ai-suggestion"><span className="suggestion-icon">✦</span><span><strong>AI 建议</strong> 当前参数适合滑膜免疫细胞数据。运行前会根据每个样本的 QC 分布再次检查阈值。</span><button>查看依据 →</button></div>
              </section>

              <section className="card annotation-card">
                <div className="card-title-row"><div><span className="section-kicker">STEP 04</span><h2>注释策略</h2></div><span className="ready-tag"><i /> 可开始</span></div>
                <div className="annotation-choice"><div className="choice-icon">✦</div><div><strong>AI 辅助注释</strong><p>结合 cluster markers、经典 marker 和组织背景生成建议。</p></div><span className="switch on"><i /></span></div>
                <div className="annotation-choice"><div className="choice-icon manual">⌁</div><div><strong>保留手动修改</strong><p>你可以在 UMAP 结果页修改任意 cluster 的最终标签。</p></div><span className="switch on"><i /></span></div>
              </section>
            </div>

            <div className="right-column">
              <section className="card results-card">
                <div className="card-title-row"><div><span className="section-kicker">LIVE PREVIEW</span><h2>结果预览</h2></div><span className={`run-tag ${runStatus}`}>{runStatus === "running" ? "分析中 · " + (activeStep + 1) + "/4" : runStatus === "complete" ? "已完成" : runStatus === "error" ? "需要处理" : "等待上传"}</span></div>
                <div className="plot-toolbar"><div className="segmented"><button className={viewMode === "celltype" ? "selected" : ""} onClick={() => setViewMode("celltype")}>细胞类型</button><button className={viewMode === "sample" ? "selected" : ""} onClick={() => setViewMode("sample")}>样本</button><button className={viewMode === "cluster" ? "selected" : ""} onClick={() => setViewMode("cluster")}>Cluster</button></div><button className="plot-action">⤢ 放大</button></div>
                <div className={`demo-note ${analysisResult?.status === "complete" ? "real-note" : ""}`}><span>{analysisResult?.status === "complete" ? "✓" : "ⓘ"}</span> {analysisResult?.status === "complete" ? "这是本次上传数据生成的真实 UMAP。Cluster 编号还不是细胞类型名称。" : "上传数据后，这里会显示真实的 UMAP；未上传时显示界面示例。"}</div>
                {analysisResult?.artifacts?.umap ? <div className="real-plot-wrap"><img className="real-plot" src={`${API_BASE}${analysisResult.artifacts.umap}`} alt="本次单细胞分析生成的 UMAP 图" /></div> : <div className="plot-area"><div className="axis-label y">UMAP_2</div><div className="axis-label x">UMAP_1</div>{points.map((point, index) => <span key={index} className="plot-point" style={{ left: `${point.left}%`, top: `${point.top}%`, background: point.color, width: point.size, height: point.size }} />)}</div>}
                <div className="plot-footer"><div className="legend">{analysisResult?.status === "complete" ? <span><i style={{ background: "#157c72" }} />真实 Cluster UMAP</span> : cellTypes.map((type) => <span key={type.name}><i style={{ background: type.color }} />{type.name}</span>)}</div><span className="plot-meta">{analysisResult?.status === "complete" ? <a href={`${API_BASE}${analysisResult.artifacts?.umapCoordinates || ""}`} target="_blank" rel="noreferrer">下载 UMAP 坐标</a> : "等待真实数据"}</span></div>
                {errorMessage && <div className="service-error">{errorMessage}</div>}
              </section>

              <section className="card insight-card"><div className="insight-head"><span className="ai-orb small">✦</span><div><span className="section-kicker">AI INSIGHT</span><h2>初步观察</h2></div><button className="more-button">···</button></div><p>当前预览显示主要细胞群之间分离良好。建议下一步重点检查 <strong>Myeloid</strong> 亚群中的 LYZ / FCER1G 表达，以及是否存在双细胞。</p><div className="insight-bottom"><span><i className="confidence-dot" /> 置信度 86%</span><button>查看 marker 证据 →</button></div></section>

              <section className="card markers-card"><div className="card-title-row"><div><span className="section-kicker">ANNOTATION EVIDENCE</span><h2>注释证据</h2></div><button className="text-link">全部查看 →</button></div><div className="marker-table"><div className="marker-row header"><span>细胞类型</span><span>代表性 markers</span><span>评分</span><span /></div>{markerRows.map((row) => <div className="marker-row" key={row[0]}><span className="marker-name"><i style={{ background: cellTypes.find((type) => type.name === row[0])?.color }} />{row[0]}</span><span className="marker-genes">{row[1]}</span><span className="score">{row[2]}</span><span className={`confidence ${row[3] === "高" ? "high" : "medium"}`}>{row[3]}</span></div>)}</div></section>
            </div>
          </div>

          <section className="run-banner"><div className="run-banner-icon">{runStatus === "running" ? "◌" : runStatus === "complete" ? "✓" : runStatus === "error" ? "!" : "✦"}</div><div><strong>{runStatus === "running" ? `正在执行：${analysisResult?.message || pipeline[activeStep][1]}` : runStatus === "complete" ? "本次真实分析已完成" : runStatus === "error" ? "分析没有完成" : "准备好开始你的第一次分析了吗？"}</strong><p>{runStatus === "running" ? `${analysisResult?.progress || 0}% · ${analysisResult?.message || pipeline[activeStep][2]}` : runStatus === "complete" ? "你可以查看真实 UMAP、QC、marker 和结果报告。" : runStatus === "error" ? errorMessage : "请先上传 10x 三文件或 zip 压缩包，再开始分析。"}</p></div><button className="primary-button" onClick={() => void startAnalysis()} disabled={runStatus === "running"}>{runStatus === "running" ? "分析进行中…" : runStatus === "complete" ? "重新运行" : "开始分析"}<span>→</span></button></section>
          <footer className="footer"><span>Scrnalysis · 为自己的数据保留完整证据链</span><span>v0.2 <i /> 本地分析服务已连接</span></footer>
        </div>
      </section>
    </main>
  );
}

