suppressPackageStartupMessages({
  library(monocle3)
  library(Matrix)
  library(hdf5r)
})

set.seed(20260811)
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
    dims = as.integer(h5[["shape"]][]), repr = "C",
    dimnames = list(as.character(h5[["genes"]][]), as.character(h5[["cells"]][]))
  )
}

counts <- read_sparse_h5(file.path(prepared_dir, "fibk14_kc_counts.h5"))
metadata <- read.delim(file.path(prepared_dir, "fibk14_kc_metadata.tsv"), check.names = FALSE)
stopifnot(identical(colnames(counts), metadata$cell))
rownames(metadata) <- metadata$cell
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
cds <- new_cell_data_set(counts, cell_metadata = metadata, gene_metadata = gene_metadata)

cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 50L, method = "PCA", norm_method = "log", verbose = TRUE)
cds <- reduce_dimension(
  cds, reduction_method = "UMAP", preprocess_method = "PCA",
  umap.metric = "cosine", umap.min_dist = 0.20, umap.n_neighbors = 30L,
  umap.fast_sgd = TRUE, cores = 12L, verbose = TRUE
)
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 30L, resolution = 1e-4,
                     random_seed = 20260811, verbose = TRUE)
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE,
                   learn_graph_control = list(minimal_branch_len = 15), verbose = TRUE)

fib_cells <- rownames(metadata)[metadata$sub_labels == "Fib_K14"]
if (!length(fib_cells)) stop("Fib_K14 root cells missing")
fib_umap <- reducedDims(cds)$UMAP[fib_cells, , drop = FALSE]
centroid <- colMeans(fib_umap)
dist2 <- rowSums((fib_umap - matrix(centroid, nrow(fib_umap), 2L, byrow=TRUE))^2)
root_cells <- names(sort(dist2))[seq_len(min(100L, length(dist2)))]
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

pt <- pseudotime(cds)
coords <- as.data.frame(reducedDims(cds)$UMAP)
colnames(coords) <- c("UMAP_1", "UMAP_2")
trajectory <- cbind(
  data.frame(cell=colnames(cds), pseudotime=as.numeric(pt), stringsAsFactors=FALSE),
  as.data.frame(colData(cds)[, c("bulk_labels","sub_labels","patient","Time_category","sample_id","orig.ident")]),
  coords
)
write.table(trajectory, gzfile(file.path(output_dir, "human_FibK14_allKC_monocle3_cells.tsv.gz"), "wt"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(cell=root_cells, root_definition="Fib_K14_UMAP_centroid"),
            file.path(output_dir, "human_FibK14_allKC_root_cells.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

focus <- intersect(c("KLF4","KRT5","KRT14","TP63","EPCAM","TACSTD2","PDGFRA","VIM","COL1A1"), rownames(cds))
expr <- normalized_counts(cds, norm_method = "log")
expr_focus <- as.matrix(expr[focus, , drop=FALSE])
expr_long <- data.frame(cell=rep(colnames(cds), each=length(focus)),
                        gene=rep(focus, times=ncol(cds)),
                        expression=as.vector(expr_focus), stringsAsFactors=FALSE)
write.table(expr_long, gzfile(file.path(output_dir, "human_FibK14_allKC_focus_expression.tsv.gz"), "wt"),
            sep="\t", quote=FALSE, row.names=FALSE)

saveRDS(cds, file.path(output_dir, "human_FibK14_allKC_monocle3_ordered.rds"), compress=FALSE)
checkpoint <- data.frame(
  status="PASS", species="human", genome_build="GRCh38", cells=ncol(cds), genes=nrow(cds),
  gate="Fib_K14 plus all bulk Keratinocyte cells", root_label="Fib_K14",
  all_cells_retained=TRUE, subset_umap_recomputed=TRUE, stored_global_umap_reused=FALSE,
  use_partition=FALSE, random_seed=20260811
)
write.table(checkpoint, file.path(output_dir, "monocle3_checkpoint.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
capture.output(sessionInfo(), file=file.path(output_dir, "sessionInfo_Monocle3.txt"))
