suppressPackageStartupMessages(library(Matrix))
out<-commandArgs(trailingOnly=TRUE)[1]
roots<-jsonlite::fromJSON(file.path(out,'analysis/roots.json'))$packages
cat('Reading rat RDS\n');flush.console()
x<-readRDS('data/external/local_360d3c8458b3')
print(class(x));print(names(attributes(x)))
print(lapply(attributes(x),function(y)class(y)))
md<-attr(x,'meta.data');print(dim(md));print(names(md));print(head(md,2))
a<-attr(x,'assays');print(names(a));print(lapply(a,function(y)names(attributes(y))))
writeLines(capture.output(list(class=class(x),attrs=names(attributes(x)),assays=lapply(a,function(y)names(attributes(y))),metadata=names(md))),file.path(out,'reports/rat_RDS_structure.txt'))
write.csv(md,file.path(out,'data/rat_scrna_obs.csv'),row.names=TRUE)
a1<-a[['RNA']];at<-attributes(a1);print(lapply(at,class))
ct<-at[['counts']]
if(is.null(ct)&&!is.null(at[['layers']])) {print(names(at[['layers']]));ct<-at[['layers']][['counts']]}
if(is.null(ct))stop('No recoverable RNA counts in raw object attributes')
print(dim(ct));print(head(rownames(ct)));cfg<-jsonlite::fromJSON(file.path(out,'analysis/defined_modules.json'))
fn<-attr(at[['features']],'dimnames')[[1]];cn<-attr(at[['cells']],'dimnames')[[1]]
stopifnot(length(fn)==nrow(ct),length(cn)==ncol(ct),identical(cn,rownames(md)))
dimnames(ct)<-list(fn,cn)
genes<-cfg$genes;ind<-match(genes,toupper(rownames(ct)));prs<-!is.na(ind)
panel<-matrix(0,ncol(ct),length(genes),dimnames=list(colnames(ct),genes))
panel[,prs]<-t(as.matrix(ct[ind[prs],,drop=FALSE]));tot<-Matrix::colSums(ct)
panel<-log1p(panel/pmax(tot,1)*1e4)
write.csv(panel,file.path(out,'data/rat_scrna_panel_logCP10K.csv'),row.names=TRUE)
write.csv(data.frame(gene=genes,present=prs,source_feature=rownames(ct)[ind],n_positive=colSums(panel>0)),file.path(out,'tables/rat_scrna_gene_availability.csv'),row.names=FALSE)
cat('RAT RNA EXTRACTION DONE\n')
