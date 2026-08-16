options(stringsAsFactors = FALSE)
set.seed(42)

argv <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  idx <- match(name, argv)
  if (!is.na(idx) && idx < length(argv)) argv[[idx + 1]] else default
}
num_arg <- function(name, default) {
  value <- suppressWarnings(as.numeric(arg_value(name, as.character(default))))
  if (is.na(value)) default else value
}

run_dir <- normalizePath(arg_value("--run-dir", "D:/ResearchPro/scrnalysis/test_runs/GSE279818_test_20260816_actual"), mustWork = FALSE)
raw_dir <- normalizePath(arg_value("--raw-dir", "D:/ResearchPro/scrnalysis/testdata/GSE279818_RAW"), mustWork = TRUE)
dataset_id <- arg_value("--dataset-id", "GSE279818_RAW")
min_features <- num_arg("--min-genes", 200)
max_features <- num_arg("--max-genes", 6000)
max_counts <- num_arg("--max-counts", Inf)
max_mito <- num_arg("--max-mito", 20)
max_ribo <- num_arg("--max-ribo", Inf)
min_cells_per_gene <- num_arg("--min-cells-per-gene", 3)
dims_requested <- num_arg("--dims", 20)
resolution_requested <- num_arg("--resolution", 0.6)
neighbors_requested <- num_arg("--neighbors", 30)
umap_neighbors_requested <- num_arg("--umap-neighbors", 30)
umap_min_dist_requested <- num_arg("--umap-min-dist", 0.3)
tsne_perplexity_requested <- num_arg("--tsne-perplexity", 30)
doublet_method <- arg_value("--doublet-method", "scDblFinder")
ambient_method <- arg_value("--ambient-method", "diagnostic_only")
fig_dir <- file.path(run_dir, "results", "figures")
table_dir <- file.path(run_dir, "results", "tables")
object_dir <- file.path(run_dir, "results", "objects")
log_dir <- file.path(run_dir, "logs")
report_dir <- file.path(run_dir, "reports")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "analysis.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

say <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""), "\n", sep = "")
  flush.console()
}
write_csv_base <- function(x, path, row.names = FALSE) {
  utils::write.csv(x, path, row.names = row.names, na = "")
}
save_plot <- function(plot_obj, path, width = 9, height = 7, dpi = 160) {
  ggplot2::ggsave(path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

say("START Scrnalysis real test analysis")
say("Raw input is read-only for this run: ", raw_dir)
say("Dataset: ", dataset_id)
say("Run contract: ", file.path(run_dir, "metadata", "run_contract.yaml"))
say("R version: ", paste(R.version$major, R.version$minor, sep = "."))

required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "readr", "tibble")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

matrix_path <- file.path(raw_dir, "matrix.mtx.gz")
barcodes_path <- file.path(raw_dir, "barcodes.tsv.gz")
features_path <- file.path(raw_dir, "features.tsv.gz")
if (!all(file.exists(c(matrix_path, barcodes_path, features_path)))) {
  stop("The testdata folder must contain matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz")
}

say("STEP 1/9 Read 10x Matrix Market files")
counts <- Seurat::ReadMtx(
  mtx = matrix_path,
  cells = barcodes_path,
  features = features_path,
  feature.column = 2,
  unique.features = TRUE
)
say("Loaded features=", nrow(counts), ", cells=", ncol(counts), ", nonzero=", length(counts@x))

obj <- Seurat::CreateSeuratObject(
  counts = counts,
  project = dataset_id,
  min.cells = min_cells_per_gene,
  min.features = 0
)
obj$sample_id <- "GSE279818_RAW"
obj$species <- "human_inferred"
obj$tissue <- "unknown"
obj$condition <- "unknown"

say("STEP 2/9 Calculate technical QC metrics")
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^(RPS|RPL)")
obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^HB")
meta_before <- obj[[]]
meta_before$cell_id <- rownames(meta_before)
meta_before$pass_qc <- with(meta_before, nFeature_RNA >= min_features & nFeature_RNA <= max_features & nCount_RNA <= max_counts & percent.mt <= max_mito & percent.ribo <= max_ribo)
qc_summary <- tibble::tibble(
  metric = c("cells_before_qc", "cells_after_qc", "genes_in_matrix", "median_genes_before", "median_umi_before", "median_mito_percent_before", "median_genes_after", "median_umi_after", "median_mito_percent_after", "qc_min_features", "qc_max_features", "qc_max_mito_percent"),
  value = c(
    nrow(meta_before), sum(meta_before$pass_qc), nrow(obj),
    stats::median(meta_before$nFeature_RNA), stats::median(meta_before$nCount_RNA), stats::median(meta_before$percent.mt),
    stats::median(meta_before$nFeature_RNA[meta_before$pass_qc]), stats::median(meta_before$nCount_RNA[meta_before$pass_qc]), stats::median(meta_before$percent.mt[meta_before$pass_qc]),
    min_features, max_features, max_mito
  )
)
write_csv_base(qc_summary, file.path(table_dir, "qc_summary.csv"))
write_csv_base(meta_before, file.path(table_dir, "cell_metadata_before_qc.csv"))

qc_plot <- ggplot2::ggplot(meta_before, ggplot2::aes(x = nFeature_RNA, y = nCount_RNA, color = percent.mt)) +
  ggplot2::geom_point(size = 0.55, alpha = 0.65) +
  ggplot2::scale_color_viridis_c(option = "C") +
  ggplot2::geom_vline(xintercept = c(min_features, max_features), linetype = "dashed", color = "#C2410C") +
  ggplot2::labs(title = "QC before filtering", x = "Detected genes per cell", y = "UMI counts per cell", color = "Mitochondrial %") +
  ggplot2::theme_minimal(base_size = 12)
save_plot(qc_plot, file.path(fig_dir, "qc_before_filtering.png"))

say("STEP 3/9 Apply fixed exploratory QC thresholds")
obj <- subset(obj, subset = nFeature_RNA >= min_features & nFeature_RNA <= max_features & nCount_RNA <= max_counts & percent.mt <= max_mito & percent.ribo <= max_ribo)
if (ncol(obj) < 50) stop("Fewer than 50 cells remained after QC; stop to avoid misleading downstream results")
obj$pass_qc <- TRUE
meta_after <- obj[[]]
meta_after$cell_id <- rownames(meta_after)
write_csv_base(meta_after, file.path(table_dir, "cell_metadata_after_qc.csv"))

qc_after_plot <- ggplot2::ggplot(meta_after, ggplot2::aes(x = nFeature_RNA, y = nCount_RNA, color = percent.mt)) +
  ggplot2::geom_point(size = 0.65, alpha = 0.7) +
  ggplot2::scale_color_viridis_c(option = "C") +
  ggplot2::labs(title = "QC after filtering", x = "Detected genes per cell", y = "UMI counts per cell", color = "Mitochondrial %") +
  ggplot2::theme_minimal(base_size = 12)
save_plot(qc_after_plot, file.path(fig_dir, "qc_after_filtering.png"))

say("STEP 4/9 Run doublet check when scDblFinder is available")
doublet_status <- "not_run"
if (tolower(doublet_method) == "none") {
  doublet_status <- "disabled_by_user"
  doublet_result <- FALSE
} else if (requireNamespace("scDblFinder", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  doublet_result <- tryCatch({
    sce <- Seurat::as.SingleCellExperiment(obj)
    sce$sample_id <- obj$sample_id[colnames(sce)]
    sce <- scDblFinder::scDblFinder(sce, samples = "sample_id", verbose = FALSE)
    pred <- as.data.frame(SummarizedExperiment::colData(sce))
    pred <- pred[colnames(obj), , drop = FALSE]
    obj$doublet_class <- as.character(pred$scDblFinder.class)
    obj$doublet_score <- as.numeric(pred$scDblFinder.score)
    doublet_status <<- paste0("completed; singlet=", sum(obj$doublet_class == "singlet", na.rm = TRUE), ", doublet=", sum(obj$doublet_class == "doublet", na.rm = TRUE))
    TRUE
  }, error = function(e) {
    say("Doublet check skipped after error: ", conditionMessage(e))
    FALSE
  })
} else {
  say("scDblFinder or SingleCellExperiment is not installed; doublet check skipped")
  doublet_status <- "unavailable"
  doublet_result <- FALSE
}
if (!doublet_result) {
  obj$doublet_class <- "not_evaluated"
  obj$doublet_score <- NA_real_
}
doublet_removed <- 0L
if (doublet_result) {
  doublet_removed <- sum(obj$doublet_class == "doublet", na.rm = TRUE)
  obj <- subset(obj, subset = doublet_class == "singlet")
  say("Removed suspected doublets before downstream embedding: ", doublet_removed, "; cells remaining=", ncol(obj))
}
meta_after_doublet <- obj[[]]
meta_after_doublet$cell_id <- rownames(meta_after_doublet)
write_csv_base(meta_after_doublet, file.path(table_dir, "cell_metadata_after_doublet.csv"))
qc_summary <- dplyr::bind_rows(qc_summary, tibble::tibble(metric = "cells_after_doublet", value = ncol(obj)))
write_csv_base(qc_summary, file.path(table_dir, "qc_summary.csv"))
write_csv_base(obj[[]] |> tibble::rownames_to_column("cell_id"), file.path(table_dir, "cell_metadata.csv"))

say("STEP 5/9 Ambient RNA diagnostic")
ambient_status <- paste0("diagnostic_only_no_correction; requested=", ambient_method)
ambient_table <- tibble::tibble(
  diagnostic = c("raw_droplet_matrix_available", "correction_applied", "reason"),
  value = c("no", "no", "Only filtered 10x matrix was supplied; raw droplet diagnostics are not possible")
)
write_csv_base(ambient_table, file.path(table_dir, "ambient_rna_diagnostic.csv"))

say("STEP 6/9 Normalize, find variable genes and run PCA")
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000, nrow(obj)), verbose = FALSE)
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), verbose = FALSE)
max_pca_dims <- min(50, ncol(obj) - 1)
npcs_to_run <- min(max(30, dims_requested), max_pca_dims)
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), npcs = npcs_to_run, verbose = FALSE)
elbow <- Seurat::ElbowPlot(obj, ndims = npcs_to_run) + ggplot2::labs(title = "PCA elbow plot")
save_plot(elbow, file.path(fig_dir, "pca_elbow.png"), width = 8, height = 6)

say("STEP 7/9 Neighbors, clusters and UMAP")
dims_use <- min(dims_requested, npcs_to_run, ncol(obj) - 1)
k_param_use <- min(neighbors_requested, ncol(obj) - 1)
umap_neighbors_use <- min(umap_neighbors_requested, ncol(obj) - 1)
obj <- Seurat::FindNeighbors(obj, dims = seq_len(dims_use), k.param = k_param_use, verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = resolution_requested, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, dims = seq_len(dims_use), n.neighbors = umap_neighbors_use, min.dist = umap_min_dist_requested, seed.use = 42, verbose = FALSE)
umap_df <- as.data.frame(Seurat::Embeddings(obj, "umap"))
umap_df$cell_id <- rownames(umap_df)
umap_df$cluster <- as.character(obj$seurat_clusters[umap_df$cell_id])
umap_df$sample_id <- as.character(obj$sample_id[umap_df$cell_id])
write_csv_base(umap_df, file.path(table_dir, "umap_coordinates.csv"))
umap_plot <- Seurat::DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggplot2::labs(title = "UMAP by unsupervised cluster", subtitle = "Exploratory run; clusters are not biological annotations")
save_plot(umap_plot, file.path(fig_dir, "umap_by_cluster.png"), width = 9, height = 7)

say("STEP 8/9 Try t-SNE and export cluster markers")
tsne_status <- "not_run"
try({
  obj <- Seurat::RunTSNE(obj, dims = seq_len(dims_use), perplexity = min(tsne_perplexity_requested, floor((ncol(obj) - 1) / 3)), seed.use = 42, verbose = FALSE)
  tsne_df <- as.data.frame(Seurat::Embeddings(obj, "tsne"))
  tsne_df$cell_id <- rownames(tsne_df)
  tsne_df$cluster <- as.character(obj$seurat_clusters[tsne_df$cell_id])
  write_csv_base(tsne_df, file.path(table_dir, "tsne_coordinates.csv"))
  tsne_plot <- Seurat::DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggplot2::labs(title = "t-SNE by unsupervised cluster")
  save_plot(tsne_plot, file.path(fig_dir, "tsne_by_cluster.png"), width = 9, height = 7)
  tsne_status <- "completed"
}, silent = TRUE)
if (tsne_status != "completed") say("t-SNE was not produced; UMAP remains the primary embedding")

markers <- tryCatch({
  Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
}, error = function(e) {
  say("Marker export warning: ", conditionMessage(e))
  data.frame()
})
if (nrow(markers) > 0) write_csv_base(markers, file.path(table_dir, "cluster_markers.csv")) else write_csv_base(data.frame(), file.path(table_dir, "cluster_markers.csv"))

say("STEP 9/9 Save object and reports")
saveRDS(obj, file.path(object_dir, "analysis_object.rds"), compress = TRUE)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
cluster_counts <- sort(table(obj$seurat_clusters), decreasing = TRUE)
marker_note <- if (nrow(markers) > 0) paste0(nrow(markers), " marker rows exported") else "marker export was empty"
report_lines <- c(
  "# Scrnalysis 测试分析报告",
  "",
  "## 本次运行结论",
  "",
  paste0("这是一份真实的 10x 单细胞矩阵测试。输入包含 ", nrow(counts), " 个基因和 ", ncol(counts), " 个细胞；QC 后保留 ", sum(meta_before$pass_qc), " 个细胞，去除疑似 doublet 后最终用于降维和聚类的细胞数为 ", ncol(obj), "。"),
  "本次运行成功完成了技术质控、LogNormalize、PCA、邻居图、无监督聚类和 UMAP。",
  "由于 testdata 没有组织、疾病/处理分组、供体信息或研究目的，本报告不把 cluster 编号直接解释为具体细胞类型。",
  "",
  "## 关键参数",
  "",
  paste0("- QC：检测基因数 ", min_features, "–", max_features, "; UMI ≤", max_counts, "; 线粒体比例 ≤", max_mito, "%; 核糖体比例 ≤", max_ribo, "%"),
  "- 归一化：LogNormalize，scale factor 10000",
  paste0("- PCA：运行 ", npcs_to_run, " 个主成分；下游使用前 ", dims_use, " 个主成分"),
  paste0("- 聚类：resolution ", resolution_requested, "; k.param ", k_param_use),
  paste0("- UMAP：n.neighbors ", umap_neighbors_use, "; min.dist ", umap_min_dist_requested, "; seed 42"),
  paste0("- Doublet（", doublet_method, "）：", doublet_status, "; removed=", doublet_removed),
  paste0("- Ambient RNA：", ambient_status),
  paste0("- t-SNE：", tsne_status),
  "",
  "## QC 与聚类概览",
  "",
  paste0("- QC 前细胞数：", nrow(meta_before)),
  paste0("- QC 后细胞数：", sum(meta_before$pass_qc)),
  paste0("- 去除疑似 doublet 后细胞数：", ncol(obj)),
  paste0("- 聚类数量：", length(cluster_counts)),
  paste0("- 各 cluster 细胞数：", paste(names(cluster_counts), as.integer(cluster_counts), sep = ":", collapse = "; ")),
  paste0("- Marker：", marker_note),
  "",
  "## 输出文件",
  "",
  "- `results/figures/umap_by_cluster.png`：UMAP 图",
  "- `results/tables/umap_coordinates.csv`：每个细胞的 UMAP 坐标和 cluster",
  "- `results/tables/qc_summary.csv`：QC 汇总",
  "- `results/tables/cell_metadata_after_doublet.csv`：去除疑似 doublet 后的细胞信息",
  "- `results/tables/cluster_markers.csv`：每个 cluster 的 marker 候选",
  "- `results/objects/analysis_object.rds`：可继续分析的 Seurat 对象",
  "",
  "## 限制",
  "",
  "这不是最终生物学结论。下一步需要补充样本背景、比较分组和细胞注释依据，再进行 marker 审核、注释和差异分析。"
)
writeLines(report_lines, file.path(report_dir, "final_analysis_report.md"), useBytes = TRUE)
say("DONE UMAP and analysis outputs written to ", file.path(run_dir, "results"))

