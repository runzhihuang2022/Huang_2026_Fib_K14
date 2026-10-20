args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "data/external/local_f41d8e2a45cb"
.libPaths(c("data/external/local_75d7abf73244", .libPaths()))
groups <- c(D0_7="Mode 1", D15_28="Mode 2", D8_14="Mode 3", M1_2="Mode 4")
rows <- list()
for (g in names(groups)) {
  p <- file.path(root, "01_input_adapter", paste0("CellChat_", g, "_sub_labels.rds"))
  x <- readRDS(p)
  ids <- as.character(x@idents)
  tab <- sort(table(ids), decreasing=TRUE)
  rows[[g]] <- data.frame(group=g, mode=groups[[g]], cell_type=names(tab), n=as.integer(tab), eligible_n10=as.integer(tab)>=10)
}
out <- do.call(rbind, rows)
write.table(out, file.path(root, "00_manifest", "ALL_MODE_CELLCHAT_SUBTYPE_COUNTS.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat("wrote", nrow(out), "rows\n")
