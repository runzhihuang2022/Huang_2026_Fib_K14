out<-commandArgs(trailingOnly=TRUE)[1]
suppressPackageStartupMessages(library(fgsea))
sets<-gmtPathways(file.path(out,'data/GO_KEGG_version_locked.gmt'))
for (nm in c('KLF4','IGFL1')) {
 d<-read.csv(file.path(out,'tables',paste0(nm,'_DESeq2_all.csv')),fileEncoding='UTF-8-BOM')
 d<-d[is.finite(d$stat)&!is.na(d$gene_name)&nzchar(d$gene_name),];d<-d[order(-d$baseMean),];d<-d[!duplicated(d$gene_name),]
 ranks<-setNames(d$stat,d$gene_name);ranks<-sort(ranks,decreasing=TRUE)
 set.seed(20261017)
 r<-fgseaMultilevel(pathways=sets,stats=ranks,minSize=10,maxSize=500,eps=1e-6,nproc=1)
 r<-as.data.frame(r);r$leadingEdge<-vapply(r$leadingEdge,paste,collapse=';',FUN.VALUE=character(1))
 r<-r[order(r$padj),];write.csv(r,file.path(out,'tables',paste0(nm,'_GSEA_all_GO_KEGG.csv')),row.names=FALSE)
 focus<-grepl('::(MAPK signaling pathway|NF-kappa B signaling pathway|PI3K-Akt signaling pathway|JAK-STAT signaling pathway|TGF-beta signaling pathway|TNF signaling pathway|ErbB signaling pathway|Notch signaling pathway|Wnt signaling pathway|Cell cycle|Adherens junction|Tight junction|Focal adhesion|Apoptosis)$',r$pathway) |
 grepl('::(Epithelium Development|Epithelial Cell Differentiation|Keratinocyte Differentiation|Epithelial Cell Migration|Positive Regulation Of Epithelial Cell Proliferation|Positive Regulation Of Epithelial Cell Migration|Wound Healing) \\(GO:',r$pathway)
 write.csv(r[focus,],file.path(out,'tables',paste0(nm,'_GSEA_mechanism_focus.csv')),row.names=FALSE)
 cat('GSEA DONE',nm,nrow(r),'sets\n');print(r[focus,c('pathway','NES','pval','padj')],row.names=FALSE)
}
