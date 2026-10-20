# Explicit optional installation in the R library selected by the user.
if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","edgeR","limma","fgsea"),ask=FALSE,update=FALSE)
# For exact historical versions, provision the recorded R 4.4/Bioconductor 3.20 environment.
# The installation command above by itself does not guarantee historical package versions.
