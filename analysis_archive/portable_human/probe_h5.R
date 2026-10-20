suppressPackageStartupMessages(library(hdf5r))
args <- commandArgs(trailingOnly = TRUE)
h5 <- H5File$new(args[[1]], mode = "r")
on.exit(h5$close_all())
dataset <- h5[["bulk_labels/counts"]]
cat("dataset_dims", paste(dataset$dims, collapse = "x"), "\n")
block <- dataset[c(1L, 2L, 3L), ]
cat("block_dims", paste(dim(block), collapse = "x"), "\n")
cat("genes", length(h5[["genes"]][]), "\n")
