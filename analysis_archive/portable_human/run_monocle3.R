suppressPackageStartupMessages({
  library(monocle3)
  library(Matrix)
  library(hdf5r)
})

set.seed(20260802)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: run_monocle3.R <prepared_dir> <output_dir>")
prepared_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_sparse_h5 <- function(path) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  sparseMatrix(
    i = as.integer(h5[["indices"]][]) + 1L,
    p = as.integer(h5[["indptr"]][]),
    x = as.numeric(h5[["data"]][]),
    dims = as.integer(h5[["shape"]][]),
    repr = "C",
    dimnames = list(
      as.character(h5[["genes"]][]),
      as.character(h5[["cells"]][])
    )
  )
}

counts <- read_sparse_h5(file.path(
  prepared_dir, "monocle3_epithelial_counts.h5"
))
metadata <- read.delim(
  file.path(prepared_dir, "monocle3_epithelial_metadata.tsv"),
  check.names = FALSE
)
umap <- read.delim(
  file.path(prepared_dir, "monocle3_epithelial_umap.tsv"),
  check.names = FALSE
)
stopifnot(
  identical(colnames(counts), metadata$cell),
  identical(metadata$cell, umap$cell)
)
rownames(metadata) <- metadata$cell
gene_metadata <- data.frame(
  gene_short_name = rownames(counts), row.names = rownames(counts)
)
cds <- new_cell_data_set(
  counts, cell_metadata = metadata, gene_metadata = gene_metadata
)
reducedDims(cds)$UMAP <- as.matrix(umap[, c("UMAP_1", "UMAP_2")])
rownames(reducedDims(cds)$UMAP) <- umap$cell

cds <- cluster_cells(
  cds, reduction_method = "UMAP", k = 20L, resolution = 1e-3,
  random_seed = 20260802, verbose = FALSE
)
cds <- learn_graph(
  cds, use_partition = TRUE, close_loop = FALSE,
  learn_graph_control = list(
    minimal_branch_len = 10, geodesic_distance_ratio = 0.5
  ),
  verbose = FALSE
)
root_pool <- rownames(metadata)[metadata$sub_labels == "KC_Basal"]
if (!length(root_pool)) stop("KC_Basal root cells missing.")
root_umap <- reducedDims(cds)$UMAP[root_pool, , drop = FALSE]
root_centroid <- colMeans(root_umap)
distance <- rowSums(
  (root_umap - matrix(
    root_centroid, nrow = nrow(root_umap), ncol = 2L, byrow = TRUE
  ))^2
)
root_cells <- names(sort(distance))[seq_len(min(50L, length(distance)))]
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
pseudotime_values <- pseudotime(cds)
trajectory <- data.frame(
  cell = colnames(cds),
  patient = colData(cds)$patient,
  group = colData(cds)$group,
  bulk_labels = colData(cds)$bulk_labels,
  sub_labels = colData(cds)$sub_labels,
  pseudotime = pseudotime_values,
  stringsAsFactors = FALSE
)
trajectory_connection <- gzfile(
  file.path(output_dir, "M09_epithelial_pseudotime.tsv.gz"), "wt"
)
write.table(
  trajectory, trajectory_connection,
  sep = "\t", quote = FALSE, row.names = FALSE
)
close(trajectory_connection)
root_table <- data.frame(cell = root_cells, root_definition = "KC_Basal_UMAP_centroid")
write.table(
  root_table, file.path(output_dir, "M09_root_cells.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

finite <- is.finite(trajectory$pseudotime)
summary_table <- aggregate(
  pseudotime ~ patient + group + sub_labels,
  data = trajectory[finite, ],
  FUN = function(x) c(
    n = length(x), median = median(x), q1 = unname(quantile(x, 0.25)),
    q3 = unname(quantile(x, 0.75))
  )
)
summary_expanded <- data.frame(
  summary_table[, c("patient", "group", "sub_labels")],
  do.call(rbind, summary_table$pseudotime)
)
write.table(
  summary_expanded, file.path(output_dir, "M09_patient_state_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

graph_result <- graph_test(
  cds, neighbor_graph = "principal_graph", cores = 8L
)
graph_result$gene <- rownames(graph_result)
graph_result <- graph_result[order(graph_result$q_value), , drop = FALSE]
graph_connection <- gzfile(
  file.path(output_dir, "M09_principal_graph_gene_test.tsv.gz"), "wt"
)
write.table(
  graph_result, graph_connection,
  sep = "\t", quote = FALSE, row.names = FALSE
)
close(graph_connection)
saveRDS(cds, file.path(output_dir, "M09_epithelial_monocle3_ordered.rds"))
checkpoint <- data.frame(
  module = "M09", status = "PASS", cells = ncol(cds), genes = nrow(cds),
  root_label = "KC_Basal",
  terminal_labels = "KC_Basal_Mig;KC_Spinous_Mig",
  stored_umap_reused = TRUE, global_embedding_overwritten = FALSE,
  k = 20L, resolution = 1e-3, minimal_branch_len = 10,
  geodesic_distance_ratio = 0.5, random_seed = 20260802
)
write.table(
  checkpoint, file.path(output_dir, "M09_monocle3_checkpoint.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_Monocle3.txt"))
