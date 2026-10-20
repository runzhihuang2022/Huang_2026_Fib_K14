args <- commandArgs(trailingOnly=TRUE)
out <- args[1]
suppressPackageStartupMessages({library(DESeq2);library(edgeR);library(limma)})
ann <- read.delim(file.path(out,'tables/gene_annotation_source.tsv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
stopifnot(all(c('gene_id','gene_name') %in% colnames(ann)))
targets <- c('KLF4','IGFL1','IGFL2','IGFL3','IGFL4','IGFLR1','TMEM149','IGF1','IGF1R','KRT14','KRT5','KRT6B','KRT17','VIM','PDGFRA','EPCAM','TACSTD2','KLF5','SOX9','MMP3','CDKN1A','MKI67','PCNA','EGFR','AREG','HBEGF','TGFA','FGF7','FGF10','TGFB1','TGFBR1','STAT3','IL6','NFKBIA','DUSP1','DUSP6','FOS','JUN','EGR1','BCL2','BAX','BIRC3','CASP8','FADD','TRADD','TNF','TNFRSF1A')
for (name in c('KLF4','IGFL1')) {
  x <- read.csv(file.path(out,'data',paste0(name,'_counts.csv')),row.names=1,check.names=FALSE,fileEncoding='UTF-8-BOM')
  stopifnot(all(as.matrix(x)>=0),all(as.matrix(x)==round(as.matrix(x))))
  group <- if(name=='KLF4') c('Control','Control','Treated','Treated') else c('Control','Control','Control','Treated','Treated')
  cd <- data.frame(condition=factor(group,levels=c('Control','Treated')),row.names=colnames(x))
  write.csv(cbind(sample_id=rownames(cd),cd,library_size=colSums(x)),file.path(out,'tables',paste0(name,'_sample_audit.csv')),row.names=FALSE)
  keep <- rowSums(x)>=10
  dds <- DESeqDataSetFromMatrix(as.matrix(x[keep,,drop=FALSE]),cd,~condition)
  dds <- DESeq(dds,quiet=FALSE)
  res <- as.data.frame(results(dds,contrast=c('condition','Treated','Control'),alpha=0.05,independentFiltering=FALSE))
  res$gene_id <- rownames(res);res$gene_name <- ann$gene_name[match(res$gene_id,ann$gene_id)]
  write.csv(res,file.path(out,'tables',paste0(name,'_DESeq2_all.csv')),row.names=FALSE)
  n <- as.data.frame(counts(dds,normalized=TRUE));n$gene_id<-rownames(n);n$gene_name<-res$gene_name
  write.csv(n,file.path(out,'tables',paste0(name,'_normalized_counts.csv')),row.names=FALSE)
  z <- res[res$gene_name %in% targets,]
  z$contrast <- paste0(name,' vs control')
  write.csv(z,file.path(out,'tables',paste0(name,'_axis_genes.csv')),row.names=FALSE)
  y <- DGEList(counts=x[keep,]);y<-calcNormFactors(y);design<-model.matrix(~condition,cd)
  y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
  ed<-topTags(glmQLFTest(fit,coef=2),n=Inf,sort.by='none')$table
  ed$gene_id<-rownames(ed);ed$gene_name<-ann$gene_name[match(ed$gene_id,ann$gene_id)]
  write.csv(ed,file.path(out,'tables',paste0(name,'_edgeR_QL_sensitivity.csv')),row.names=FALSE)
  pdf(file.path(out,'figures',paste0(name,'_bulk_QC.pdf')),width=8,height=7)
  par(mfrow=c(2,2));plotDispEsts(dds);DESeq2::plotMA(dds,ylim=c(-5,5),main=name)
  v<-varianceStabilizingTransformation(dds,blind=TRUE)
  pc<-prcomp(t(assay(v)));plot(pc$x[,1:2],col=c('#3399bb','#f06455')[as.integer(cd$condition)],pch=19,xlab='PC1',ylab='PC2',main=paste(name,'VST PCA'));text(pc$x[,1:2],labels=colnames(x),pos=3,cex=.8)
  plotMDS(y,labels=colnames(x),col=c('#3399bb','#f06455')[as.integer(cd$condition)]);dev.off()
  saveRDS(dds,file.path(out,'analysis',paste0(name,'_DESeq2.rds')))
  cat('\nBULK DONE',name,'tested',nrow(res),'FDR05',sum(res$padj<.05,na.rm=TRUE),'\n')
  print(z[,c('gene_name','baseMean','log2FoldChange','pvalue','padj')],row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(out,'logs/R_sessionInfo.txt'))
