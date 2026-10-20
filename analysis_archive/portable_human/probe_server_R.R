packages <- c(
  "DESeq2", "edgeR", "limma", "clusterProfiler", "fgsea", "ReactomePA",
  "org.Hs.eg.db", "hdf5r", "rhdf5", "CellChat", "monocle3"
)
result <- data.frame(
  package = packages,
  available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(
    packages,
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(packageVersion(package))
    },
    character(1)
  )
)
cat("libPaths\n")
writeLines(.libPaths())
cat("packages\n")
write.table(result, row.names = FALSE, sep = "\t", quote = FALSE)
