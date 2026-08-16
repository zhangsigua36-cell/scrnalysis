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
species_requested <- arg_value("--species", "human")
tissue_requested <- arg_value("--tissue", "synovium")
normalization_requested <- arg_value("--normalization", "log")
batch_correction_requested <- arg_value("--batch-correction", "harmony")
annotation_scope <- arg_value("--annotation-scope", "major")
ai_annotation_enabled <- tolower(arg_value("--ai-annotation", "true")) == "true"
manual_edits_enabled <- tolower(arg_value("--manual-edits", "true")) == "true"
output_figures <- unique(strsplit(arg_value("--output-figures", "umap,tsne,dotplot,featureplot,qc_violin"), ",", fixed = TRUE)[[1]])
gene_filter_mode <- tolower(arg_value("--gene-filter-mode", "quantile"))
gene_lower_quantile <- num_arg("--gene-lower-quantile", 0.025)
gene_upper_quantile <- num_arg("--gene-upper-quantile", 0.025)
min_features_fixed <- num_arg("--min-genes", 200)
max_features_fixed <- num_arg("--max-genes", 6000)
umi_filter_mode <- tolower(arg_value("--umi-filter-mode", "quantile"))
umi_lower_quantile <- num_arg("--umi-lower-quantile", 0.025)
umi_upper_quantile <- num_arg("--umi-upper-quantile", 0.025)
min_counts_fixed <- num_arg("--min-counts", 0)
max_counts_fixed <- num_arg("--max-counts", Inf)
mito_filter_mode <- tolower(arg_value("--mito-filter-mode", "fixed"))
max_mito <- num_arg("--max-mito", 20)
ribo_filter_mode <- tolower(arg_value("--ribo-filter-mode", "none"))
max_ribo <- num_arg("--max-ribo", 30)
hb_filter_mode <- tolower(arg_value("--hb-filter-mode", "none"))
max_hb <- num_arg("--max-hb", 20)
min_cells_per_gene_mode <- tolower(arg_value("--min-cells-per-gene-mode", "fixed"))
min_cells_per_gene <- num_arg("--min-cells-per-gene", 3)
dims_requested <- num_arg("--dims", 20)
resolution_requested <- num_arg("--resolution", 0.6)
neighbors_requested <- num_arg("--neighbors", 30)
umap_neighbors_requested <- num_arg("--umap-neighbors", 30)
umap_min_dist_requested <- num_arg("--umap-min-dist", 0.3)
tsne_perplexity_requested <- num_arg("--tsne-perplexity", 30)
doublet_method <- arg_value("--doublet-method", "scDblFinder")
ambient_method <- arg_value("--ambient-method", "diagnostic_only")

tail_quantile <- function(values, tail_probability, upper = FALSE) {
  values <- values[is.finite(values)]
  if (!length(values) || !is.finite(tail_probability)) return(if (upper) Inf else -Inf)
  tail_probability <- max(0, min(0.49, tail_probability))
  probability <- if (upper) 1 - tail_probability else tail_probability
  as.numeric(stats::quantile(values, probs = probability, na.rm = TRUE, names = FALSE, type = 7))
}
resolve_lower_bound <- function(mode, fixed, values, tail_probability) {
  if (mode == "none") -Inf else if (mode == "quantile") tail_quantile(values, tail_probability, upper = FALSE) else fixed
}
resolve_upper_bound <- function(mode, fixed, values, tail_probability) {
  if (mode == "none") Inf else if (mode == "quantile") tail_quantile(values, tail_probability, upper = TRUE) else fixed
}
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
figure_enabled <- function(name) name %in% output_figures

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
  min.cells = if (min_cells_per_gene_mode == "none") 0 else min_cells_per_gene,
  min.features = 0
)
obj$sample_id <- "GSE279818_RAW"
obj$species <- species_requested
obj$tissue <- tissue_requested
obj$condition <- "unknown"

say("STEP 2/9 Calculate technical QC metrics")
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^(RPS|RPL)")
obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^HB")
meta_before <- obj[[]]
meta_before$cell_id <- rownames(meta_before)
resolved_min_features <- resolve_lower_bound(gene_filter_mode, min_features_fixed, meta_before$nFeature_RNA, gene_lower_quantile)
resolved_max_features <- resolve_upper_bound(gene_filter_mode, max_features_fixed, meta_before$nFeature_RNA, gene_upper_quantile)
resolved_min_counts <- resolve_lower_bound(umi_filter_mode, min_counts_fixed, meta_before$nCount_RNA, umi_lower_quantile)
resolved_max_counts <- resolve_upper_bound(umi_filter_mode, max_counts_fixed, meta_before$nCount_RNA, umi_upper_quantile)
resolved_max_mito <- if (mito_filter_mode == "none") Inf else max_mito
resolved_max_ribo <- if (ribo_filter_mode == "none") Inf else max_ribo
resolved_max_hb <- if (hb_filter_mode == "none") Inf else max_hb
meta_before$pass_qc <- meta_before$nFeature_RNA >= resolved_min_features & meta_before$nFeature_RNA <= resolved_max_features & meta_before$nCount_RNA >= resolved_min_counts & meta_before$nCount_RNA <= resolved_max_counts & meta_before$percent.mt <= resolved_max_mito & meta_before$percent.ribo <= resolved_max_ribo & meta_before$percent.hb <= resolved_max_hb
meta_before$pass_qc[is.na(meta_before$pass_qc)] <- FALSE
qc_summary <- tibble::tibble(
  metric = c("cells_before_qc", "cells_after_qc", "genes_in_matrix", "median_genes_before", "median_umi_before", "median_mito_percent_before", "median_genes_after", "median_umi_after", "median_mito_percent_after", "qc_gene_filter_mode", "qc_min_features", "qc_max_features", "qc_umi_filter_mode", "qc_min_umi", "qc_max_umi", "qc_mito_filter_mode", "qc_max_mito_percent", "qc_ribo_filter_mode", "qc_max_ribo_percent", "qc_hb_filter_mode", "qc_max_hb_percent", "qc_min_cells_per_gene_mode", "qc_min_cells_per_gene"),
  value = as.character(c(
    nrow(meta_before), sum(meta_before$pass_qc), nrow(obj),
    stats::median(meta_before$nFeature_RNA), stats::median(meta_before$nCount_RNA), stats::median(meta_before$percent.mt),
    stats::median(meta_before$nFeature_RNA[meta_before$pass_qc]), stats::median(meta_before$nCount_RNA[meta_before$pass_qc]), stats::median(meta_before$percent.mt[meta_before$pass_qc]),
    gene_filter_mode, resolved_min_features, resolved_max_features, umi_filter_mode, resolved_min_counts, resolved_max_counts, mito_filter_mode, resolved_max_mito, ribo_filter_mode, resolved_max_ribo, hb_filter_mode, resolved_max_hb, min_cells_per_gene_mode, min_cells_per_gene
  ))
)
write_csv_base(qc_summary, file.path(table_dir, "qc_summary.csv"))
write_csv_base(meta_before, file.path(table_dir, "cell_metadata_before_qc.csv"))

qc_plot <- ggplot2::ggplot(meta_before, ggplot2::aes(x = nFeature_RNA, y = nCount_RNA, color = percent.mt)) +
  ggplot2::geom_point(size = 0.55, alpha = 0.65) +
  ggplot2::scale_color_viridis_c(option = "C") +
  ggplot2::geom_vline(xintercept = c(resolved_min_features, resolved_max_features)[is.finite(c(resolved_min_features, resolved_max_features))], linetype = "dashed", color = "#C2410C") +
  ggplot2::labs(title = "QC before filtering", x = "Detected genes per cell", y = "UMI counts per cell", color = "Mitochondrial %") +
  ggplot2::theme_minimal(base_size = 12)
if (figure_enabled("qc_violin")) {
  qc_before_long <- rbind(
    data.frame(value = meta_before$nFeature_RNA, metric = "nFeature_RNA"),
    data.frame(value = meta_before$nCount_RNA, metric = "nCount_RNA"),
    data.frame(value = meta_before$percent.mt, metric = "percent.mt")
  )
  qc_violin_before <- ggplot2::ggplot(qc_before_long, ggplot2::aes(x = metric, y = value, fill = metric)) + ggplot2::geom_violin(scale = "width", trim = TRUE) + ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.25) + ggplot2::labs(title = "QC violin plots before filtering", x = NULL, y = "Value") + ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "none")
  save_plot(qc_violin_before, file.path(fig_dir, "qc_violin_before_filtering.png"), width = 9, height = 6)
}

say("STEP 3/9 Apply configured QC thresholds")
obj <- subset(obj, subset = nFeature_RNA >= resolved_min_features & nFeature_RNA <= resolved_max_features & nCount_RNA >= resolved_min_counts & nCount_RNA <= resolved_max_counts & percent.mt <= resolved_max_mito & percent.ribo <= resolved_max_ribo & percent.hb <= resolved_max_hb)
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
if (figure_enabled("qc_violin")) {
  qc_after_long <- rbind(
    data.frame(value = meta_after$nFeature_RNA, metric = "nFeature_RNA"),
    data.frame(value = meta_after$nCount_RNA, metric = "nCount_RNA"),
    data.frame(value = meta_after$percent.mt, metric = "percent.mt")
  )
  qc_violin_after <- ggplot2::ggplot(qc_after_long, ggplot2::aes(x = metric, y = value, fill = metric)) + ggplot2::geom_violin(scale = "width", trim = TRUE) + ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.25) + ggplot2::labs(title = "QC violin plots after filtering", x = NULL, y = "Value") + ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "none")
  save_plot(qc_violin_after, file.path(fig_dir, "qc_violin_after_filtering.png"), width = 9, height = 6)
}

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
qc_summary <- dplyr::bind_rows(qc_summary, tibble::tibble(metric = "cells_after_doublet", value = as.character(ncol(obj))))
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
normalization_used <- "LogNormalize"
if (tolower(normalization_requested) == "sct" && requireNamespace("sctransform", quietly = TRUE)) {
  obj <- Seurat::SCTransform(obj, verbose = FALSE)
  normalization_used <- "SCTransform"
} else {
  if (tolower(normalization_requested) == "sct") say("SCTransform requested but sctransform is unavailable; falling back to LogNormalize")
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
}
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
if (figure_enabled("umap")) save_plot(umap_plot, file.path(fig_dir, "umap_by_cluster.png"), width = 9, height = 7)

say("STEP 8/9 Try t-SNE and export cluster markers")
tsne_status <- "not_run"
try({
  obj <- Seurat::RunTSNE(obj, dims = seq_len(dims_use), perplexity = min(tsne_perplexity_requested, floor((ncol(obj) - 1) / 3)), seed.use = 42, verbose = FALSE)
  tsne_df <- as.data.frame(Seurat::Embeddings(obj, "tsne"))
  tsne_df$cell_id <- rownames(tsne_df)
  tsne_df$cluster <- as.character(obj$seurat_clusters[tsne_df$cell_id])
  write_csv_base(tsne_df, file.path(table_dir, "tsne_coordinates.csv"))
  tsne_plot <- Seurat::DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggplot2::labs(title = "t-SNE by unsupervised cluster")
  if (figure_enabled("tsne")) save_plot(tsne_plot, file.path(fig_dir, "tsne_by_cluster.png"), width = 9, height = 7)
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

say("STEP 8/9 Build marker-evidence annotation and requested plots")
marker_sets <- list(
  "T cells" = c("CD3D", "CD3E", "TRBC1", "IL7R", "LTB"),
  "NK cells" = c("NKG7", "GNLY", "FCGR3A"),
  "Myeloid" = c("LYZ", "LST1", "FCER1G", "CTSS"),
  "B cells" = c("CD79A", "MS4A1", "CD74", "CD37"),
  "Plasma cells" = c("MZB1", "JCHAIN", "SDC1"),
  "Fibroblast / stromal" = c("COL1A1", "COL3A1", "DCN", "LUM"),
  "Endothelial" = c("PECAM1", "VWF", "KDR"),
  "Epithelial" = c("EPCAM", "KRT8", "KRT18", "KRT19")
)
annotation_evidence <- data.frame(cluster = character(), celltype_major = character(), markers = character(), score = numeric(), confidence = character(), stringsAsFactors = FALSE)
cluster_labels <- setNames(rep("Unassigned", length(levels(obj$seurat_clusters))), levels(obj$seurat_clusters))
if (ai_annotation_enabled) {
  avg_expr <- tryCatch(Seurat::AverageExpression(obj, assays = DefaultAssay(obj), group.by = "seurat_clusters", slot = "data", verbose = FALSE)[[DefaultAssay(obj)]], error = function(e) NULL)
  if (!is.null(avg_expr)) {
    for (cluster in names(cluster_labels)) {
      column <- which(colnames(avg_expr) %in% cluster)
      if (!length(column)) column <- grep(paste0("(^|_)", cluster, "$"), colnames(avg_expr))
      score_values <- vapply(marker_sets, function(genes) {
        available <- intersect(genes, rownames(avg_expr))
        if (!length(available) || !length(column)) 0 else mean(as.numeric(avg_expr[available, column, drop = TRUE]), na.rm = TRUE)
      }, numeric(1))
      best_name <- names(which.max(score_values))
      best_score <- max(score_values, na.rm = TRUE)
      second_score <- if (length(score_values) > 1) sort(score_values, decreasing = TRUE)[2] else 0
      confidence_score <- min(0.99, max(0.05, 0.5 + (best_score - second_score) / (abs(best_score) + 1)))
      confidence_label <- if (confidence_score >= 0.75) "高" else if (confidence_score >= 0.55) "中" else "低"
      if (best_score > 0) cluster_labels[[cluster]] <- best_name
      marker_rows_for_cluster <- if (nrow(markers) > 0 && "cluster" %in% colnames(markers)) markers[as.character(markers$cluster) == cluster, , drop = FALSE] else data.frame()
      if (nrow(marker_rows_for_cluster) > 0 && "avg_log2FC" %in% colnames(marker_rows_for_cluster)) marker_rows_for_cluster <- marker_rows_for_cluster[order(marker_rows_for_cluster$avg_log2FC, decreasing = TRUE), , drop = FALSE]
      evidence_genes <- if (nrow(marker_rows_for_cluster) > 0 && "gene" %in% colnames(marker_rows_for_cluster)) head(as.character(marker_rows_for_cluster$gene), 5) else character()
      evidence_genes <- unique(c(evidence_genes, intersect(marker_sets[[best_name]], rownames(avg_expr))))
      annotation_evidence <- rbind(annotation_evidence, data.frame(cluster = cluster, celltype_major = cluster_labels[[cluster]], markers = paste(head(evidence_genes, 8), collapse = ", "), score = round(confidence_score, 3), confidence = confidence_label, stringsAsFactors = FALSE))
    }
  }
}
if (!nrow(annotation_evidence)) annotation_evidence <- data.frame(cluster = names(cluster_labels), celltype_major = "Unassigned", markers = "", score = 0, confidence = "低", stringsAsFactors = FALSE)
obj$celltype_major <- unname(cluster_labels[as.character(obj$seurat_clusters)])
obj$celltype_detail <- if (annotation_scope == "major") obj$celltype_major else obj$celltype_major
write_csv_base(annotation_evidence, file.path(table_dir, "annotation_evidence.csv"))
write_csv_base(obj[[]] |> tibble::rownames_to_column("cell_id"), file.path(table_dir, "cell_metadata.csv"))

if (figure_enabled("umap")) {
  umap_celltype <- Seurat::DimPlot(obj, reduction = "umap", group.by = "celltype_major", label = TRUE, repel = TRUE) + ggplot2::labs(title = "UMAP by marker-evidence cell type")
  umap_sample <- Seurat::DimPlot(obj, reduction = "umap", group.by = "sample_id") + ggplot2::labs(title = "UMAP by sample")
  save_plot(umap_celltype, file.path(fig_dir, "umap_by_celltype.png"), width = 9, height = 7)
  save_plot(umap_sample, file.path(fig_dir, "umap_by_sample.png"), width = 9, height = 7)
}
if (figure_enabled("tsne") && tsne_status == "completed") {
  tsne_celltype <- Seurat::DimPlot(obj, reduction = "tsne", group.by = "celltype_major", label = TRUE, repel = TRUE) + ggplot2::labs(title = "t-SNE by marker-evidence cell type")
  tsne_sample <- Seurat::DimPlot(obj, reduction = "tsne", group.by = "sample_id") + ggplot2::labs(title = "t-SNE by sample")
  save_plot(tsne_celltype, file.path(fig_dir, "tsne_by_celltype.png"), width = 9, height = 7)
  save_plot(tsne_sample, file.path(fig_dir, "tsne_by_sample.png"), width = 9, height = 7)
}
available_features <- rownames(obj)
canonical_features <- intersect(c("CD3D", "IL7R", "NKG7", "LYZ", "MS4A1", "COL1A1", "EPCAM", "HLA-DRA"), available_features)
if (figure_enabled("dotplot") && length(canonical_features) >= 2) {
  dotplot <- Seurat::DotPlot(obj, features = canonical_features, group.by = "celltype_major") + Seurat::RotatedAxis() + ggplot2::labs(title = "Canonical marker DotPlot")
  save_plot(dotplot, file.path(fig_dir, "dotplot_markers.png"), width = 10, height = 7)
}
if (figure_enabled("featureplot") && length(canonical_features) >= 1) {
  featureplot <- Seurat::FeaturePlot(obj, features = head(canonical_features, 4), reduction = "umap", combine = TRUE) + ggplot2::labs(title = "Canonical marker FeaturePlot")
  save_plot(featureplot, file.path(fig_dir, "featureplot_markers.png"), width = 11, height = 8)
}

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
  paste0("本次运行完成了技术质控、", normalization_used, "、PCA、邻居图、无监督聚类和 UMAP。组织：", tissue_requested, "；物种：", species_requested, "；批次校正设置：", batch_correction_requested, "。"),
  "由于 testdata 没有组织、疾病/处理分组、供体信息或研究目的，本报告不把 cluster 编号直接解释为具体细胞类型。",
  "",
  "## 关键参数",
  "",
  paste0("- QC：基因数模式 ", gene_filter_mode, "（", resolved_min_features, "–", resolved_max_features, "）；UMI 模式 ", umi_filter_mode, "（", resolved_min_counts, "–", resolved_max_counts, "）；线粒体模式 ", mito_filter_mode, "（上限 ", resolved_max_mito, "）；核糖体模式 ", ribo_filter_mode, "（上限 ", resolved_max_ribo, "）；Hb 模式 ", hb_filter_mode, "（上限 ", resolved_max_hb, "）。"),
  paste0("- 归一化：", normalization_used, "；scale factor 10000（如适用）"),
  paste0("- PCA：运行 ", npcs_to_run, " 个主成分；下游使用前 ", dims_use, " 个主成分"),
  paste0("- 聚类：resolution ", resolution_requested, "; k.param ", k_param_use),
  paste0("- UMAP：n.neighbors ", umap_neighbors_use, "; min.dist ", umap_min_dist_requested, "; seed 42"),
  paste0("- Doublet（", doublet_method, "）：", doublet_status, "; removed=", doublet_removed),
  paste0("- Ambient RNA：", ambient_status),
  paste0("- t-SNE：", tsne_status),
  paste0("- 注释：AI 注释 ", ifelse(ai_annotation_enabled, "开启", "关闭"), "；手动修改入口 ", ifelse(manual_edits_enabled, "保留", "关闭"), "；注释粒度 ", annotation_scope),
  paste0("- 输出图形：", paste(output_figures, collapse = ", ")), 
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
