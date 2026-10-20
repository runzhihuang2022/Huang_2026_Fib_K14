suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

set.seed(20260801)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: run_gsea.R <formal_dir> <output_dir>")
formal_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
index <- read.delim(
  file.path(formal_dir, "M05_result_index.tsv"), check.names = FALSE
)
checkpoints <- list()
for (row_index in seq_len(nrow(index))) {
  record <- index[row_index, ]
  stem <- paste(
    gsub("[^A-Za-z0-9._-]+", "_", record$track),
    gsub("[^A-Za-z0-9._-]+", "_", record$label),
    record$comparison, record$analysis, "GSEA_GO_BP", sep = "__"
  )
  output_path <- file.path(output_dir, paste0(stem, ".tsv"))
  if (file.exists(output_path)) {
    existing <- read.delim(output_path, check.names = FALSE)
    checkpoints[[length(checkpoints) + 1L]] <- data.frame(
      track = record$track, label = record$label,
      comparison = record$comparison, analysis = record$analysis,
      status = "PASS", pathways = nrow(existing),
      reason = "resumed_existing_result"
    )
    next
  }
  result <- read.delim(record$result_file, check.names = FALSE)
  result <- result[is.finite(result$stat), , drop = FALSE]
  mapped <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = unique(result$gene),
    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  ))
  result$entrez <- unname(mapped[result$gene])
  result <- result[!is.na(result$entrez), , drop = FALSE]
  result <- result[order(abs(result$stat), decreasing = TRUE), , drop = FALSE]
  result <- result[!duplicated(result$entrez), , drop = FALSE]
  ranks <- result$stat
  names(ranks) <- result$entrez
  ranks <- sort(ranks, decreasing = TRUE)
  gsea <- tryCatch(
    suppressMessages(gseGO(
      geneList = ranks, OrgDb = org.Hs.eg.db, ont = "BP",
      keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500,
      pvalueCutoff = 0.05, pAdjustMethod = "BH",
      eps = 1e-10, by = "fgsea", seed = TRUE, verbose = FALSE
    )),
    error = function(error) error
  )
  if (inherits(gsea, "error")) {
    checkpoints[[length(checkpoints) + 1L]] <- data.frame(
      track = record$track, label = record$label,
      comparison = record$comparison, analysis = record$analysis,
      status = "FAILED", pathways = NA_integer_,
      reason = conditionMessage(gsea)
    )
    next
  }
  table <- as.data.frame(gsea)
  write.table(
    table, output_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  checkpoints[[length(checkpoints) + 1L]] <- data.frame(
    track = record$track, label = record$label,
    comparison = record$comparison, analysis = record$analysis,
    status = "PASS", pathways = nrow(table), reason = ""
  )
}
checkpoint <- do.call(rbind, checkpoints)
write.table(
  checkpoint, file.path(output_dir, "M06_GSEA_checkpoint.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_GSEA.txt"))
