suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(hdf5r)
  library(future)
})

set.seed(20260801)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: run_cellchat.R <prepared_dir> <output_dir>")
prepared_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
groups <- c("Normal", "SPTDI", "DPTDI")
tracks <- c("bulk_labels", "sub_labels")
if (.Platform$OS.type == "unix") {
  # Forked workers share the large read-only sparse matrix by copy-on-write.
  future::plan("multicore", workers = 8L)
} else {
  future::plan("multisession", workers = 8L)
}
options(future.globals.maxSize = 80 * 1024^3)

read_sparse_h5 <- function(path) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  shape <- as.integer(h5[["shape"]][])
  data <- as.numeric(h5[["data"]][])
  indices <- as.integer(h5[["indices"]][])
  indptr <- as.integer(h5[["indptr"]][])
  genes <- as.character(h5[["genes"]][])
  cells <- as.character(h5[["cells"]][])
  matrix <- sparseMatrix(
    i = indices + 1L, p = indptr, x = data,
    dims = shape, repr = "C",
    dimnames = list(genes, cells)
  )
  matrix
}

write_tsv <- function(value, path) {
  write.table(
    value, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
  )
}

checkpoint <- list()
for (group in groups) {
  message("Loading ", group, " count matrix...")
  counts <- read_sparse_h5(file.path(
    prepared_dir, paste0("cellchat_", group, "_counts.h5")
  ))
  metadata <- read.delim(
    file.path(prepared_dir, paste0("cellchat_", group, "_metadata.tsv")),
    check.names = FALSE
  )
  stopifnot(identical(colnames(counts), metadata$cell))
  library_sizes <- Matrix::colSums(counts)
  if (any(library_sizes <= 0)) stop("Zero-library cell reached CellChat.")
  normalized <- counts %*% Diagonal(x = 10000 / library_sizes)
  normalized@x <- log1p(normalized@x)
  for (track in tracks) {
    run_dir <- file.path(output_dir, track, group)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    label_counts <- table(metadata[[track]])
    valid_labels <- names(label_counts[label_counts >= 5L])
    keep <- metadata[[track]] %in% valid_labels
    run_meta <- metadata[keep, , drop = FALSE]
    rownames(run_meta) <- run_meta$cell
    run_matrix <- normalized[, keep, drop = FALSE]
    result <- tryCatch({
      cellchat <- createCellChat(
        object = run_matrix, meta = run_meta, group.by = track
      )
      cellchat@DB <- CellChatDB.human
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(
        cellchat, type = "triMean", population.size = FALSE, raw.use = TRUE
      )
      cellchat <- filterCommunication(cellchat, min.cells = 5L)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      lr <- subsetCommunication(cellchat, slot.name = "net", thresh = 0.05)
      pathway <- subsetCommunication(
        cellchat, slot.name = "netP", thresh = 0.05
      )
      write_tsv(lr, file.path(run_dir, "ligand_receptor_interactions.tsv"))
      write_tsv(pathway, file.path(run_dir, "pathway_interactions.tsv"))
      count_matrix <- as.data.frame(cellchat@net$count, check.names = FALSE)
      count_matrix$source <- rownames(count_matrix)
      weight_matrix <- as.data.frame(cellchat@net$weight, check.names = FALSE)
      weight_matrix$source <- rownames(weight_matrix)
      write_tsv(count_matrix, file.path(run_dir, "interaction_count_matrix.tsv"))
      write_tsv(weight_matrix, file.path(run_dir, "interaction_weight_matrix.tsv"))
      write_tsv(
        data.frame(label = names(label_counts), cells = as.integer(label_counts)),
        file.path(run_dir, "cell_group_counts.tsv")
      )
      saveRDS(cellchat, file.path(run_dir, paste0(
        "cellchat_", track, "_", group, ".rds"
      )))
      data.frame(
        track = track, group = group, status = "PASS",
        cells = ncol(run_matrix), labels = length(valid_labels),
        lr_interactions = nrow(lr), pathway_interactions = nrow(pathway),
        compute_type = "triMean_default", population_size = FALSE,
        min_cells = 5L, normalization = "module_only_log1p_CP10K",
        prohibited_global_recomputation = FALSE, reason = ""
      )
    }, error = function(error) {
      data.frame(
        track = track, group = group, status = "FAILED",
        cells = ncol(run_matrix), labels = length(valid_labels),
        lr_interactions = NA_integer_, pathway_interactions = NA_integer_,
        compute_type = "triMean_default", population_size = FALSE,
        min_cells = 5L, normalization = "module_only_log1p_CP10K",
        prohibited_global_recomputation = FALSE,
        reason = conditionMessage(error)
      )
    })
    checkpoint[[length(checkpoint) + 1L]] <- result
    rm(run_matrix)
    gc(verbose = FALSE)
  }
  rm(counts, normalized)
  gc(verbose = FALSE)
}
checkpoint <- do.call(rbind, checkpoint)
write_tsv(checkpoint, file.path(output_dir, "M08_cellchat_checkpoint.tsv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_CellChat.txt"))
if (any(checkpoint$status != "PASS")) quit(save = "no", status = 1L)
