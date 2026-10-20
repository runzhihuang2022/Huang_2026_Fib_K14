#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ArchR);library(Seurat);library(Matrix);library(data.table);library(SummarizedExperiment)})
set.seed(20260811); addArchRThreads(threads=24)
ROOT <- 'data/external/local_6f6697f200d0'
OUT <- file.path(ROOT,'working','stage2_joint'); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
RNA <- file.path(ROOT,'working','stage1_rna','integration_reference')
ATAC_RDS <- 'data/external/local_3145112678a0'
LIBMAP <- file.path(ROOT,'inputs','atac_library_map.csv')
msg <- function(...) {cat(format(Sys.time()),...,"\n");flush.console()}

msg('READ RNA REFERENCE')
rna_ck <- file.path(OUT,'checkpoint_rna_reference.rds')
rmeta <- fread(file.path(RNA,'metadata.csv'),data.table=FALSE,check.names=FALSE); rownames(rmeta) <- rmeta$cell_barcode
if(file.exists(rna_ck)){
  rna <- readRDS(rna_ck); msg('RNA CHECKPOINT RESTORED')
} else {
  rm <- readMM(file.path(RNA,'counts_genes_by_cells.mtx'))
  genes <- fread(file.path(RNA,'genes.tsv'))$gene_symbol
  rownames(rm) <- make.unique(genes); colnames(rm) <- rmeta$cell_barcode
  rna <- CreateSeuratObject(rm,meta.data=rmeta,min.cells=0,min.features=0,project='rat_RNA_reference')
  rna <- NormalizeData(rna,verbose=FALSE); rna <- FindVariableFeatures(rna,nfeatures=3000,verbose=FALSE)
  rna <- ScaleData(rna,features=VariableFeatures(rna),verbose=FALSE); rna <- RunPCA(rna,npcs=40,verbose=FALSE)
  rna <- RunUMAP(rna,reduction='pca',dims=1:30,return.model=TRUE,reduction.name='umap.reference',seed.use=20260811,verbose=FALSE)
  saveRDS(rna,rna_ck)
}

msg('READ ARCHR CHILD INPUT')
proj <- readRDS(ATAC_RDS)
md <- as.data.frame(proj@cellColData)
md$cell_id <- rownames(md)
qc <- as.numeric(md$TSSEnrichment)>=4 & as.numeric(md$nFrags)>=1000
if ('PassQC' %in% colnames(md)) qc <- qc & as.logical(md$PassQC)
if ('analysis_qc_pass' %in% colnames(md)) qc <- qc & as.logical(md$analysis_qc_pass)
keep <- md$cell_id[qc]
if (length(keep)<nCells(proj)) proj <- proj[keep,]
md <- as.data.frame(proj@cellColData); md$cell_id <- rownames(md)
libmap <- fread(LIBMAP,data.table=FALSE)
md$library_id <- if ('library_id' %in% colnames(md)) as.character(md$library_id) else as.character(md$Sample)
md <- merge(md,libmap,by='library_id',all.x=TRUE,sort=FALSE)
rownames(md) <- md$cell_id
md <- md[getCellNames(proj),,drop=FALSE]
write.csv(md,file.path(OUT,'atac_qc_cell_metadata.csv'),row.names=FALSE)
qdt <- as.data.table(md)
qsum <- qdt[,.(n_cells=.N,median_TSS=median(as.numeric(TSSEnrichment),na.rm=TRUE),
               median_fragments=median(as.numeric(nFrags),na.rm=TRUE)),by=.(library_id,timepoint,donor_id)]
write.csv(qsum,file.path(OUT,'atac_qc_by_library.csv'),row.names=FALSE)

msg('EXTRACT GENE SCORE MATRIX')
gse <- getMatrixFromProject(proj,useMatrix='GeneScoreMatrix')
gs <- assays(gse)[[1]]
gn <- if ('name' %in% colnames(rowData(gse))) as.character(rowData(gse)$name) else rownames(gse)
rownames(gs) <- make.unique(gn); colnames(gs) <- colnames(gse)
common <- intersect(rownames(rna),rownames(gs))
rna2 <- subset(rna,features=common)
atac <- CreateSeuratObject(gs[common,,drop=FALSE],meta.data=md[colnames(gs),,drop=FALSE],min.cells=0,min.features=0,project='rat_ATAC_gene_activity')
atac <- NormalizeData(atac,verbose=FALSE)

msg('CCA ANCHORS AND LABEL TRANSFER')
map_meta_file <- file.path(OUT,'checkpoint_atac_mapped_metadata.csv.gz')
map_umap_file <- file.path(OUT,'checkpoint_atac_mapped_umap.csv.gz')
if(file.exists(map_meta_file) && file.exists(map_umap_file)){
  mm <- fread(map_meta_file,data.table=FALSE,check.names=FALSE); rownames(mm) <- mm$cell_id; mm$cell_id <- NULL
  atac@meta.data <- mm[colnames(atac),,drop=FALSE]
  uu <- fread(map_umap_file,data.table=FALSE); rownames(uu) <- uu$cell_id; uu$cell_id <- NULL
  atac[['ref.umap']] <- CreateDimReducObject(embeddings=as.matrix(uu[colnames(atac),,drop=FALSE]),key='refUMAP_',assay='RNA')
  anchors <- NULL
  msg('MAPPED CHECKPOINT RESTORED')
} else {
  anchors <- FindTransferAnchors(reference=rna2,query=atac,features=intersect(VariableFeatures(rna2),common),
                                 reduction='cca',dims=1:30,k.anchor=5,verbose=FALSE)
  atac <- MapQuery(anchorset=anchors,reference=rna2,query=atac,refdata=list(
    bulk_labels='bulk_labels',sub_labels='sub_labels_reannotated',Fib_K14_strict='Fib_K14_strict_gate'),
    reference.reduction='pca',reduction.model='umap.reference',verbose=FALSE)
  msg('MAPPED REDUCTIONS',paste(names(atac@reductions),collapse=','))
  fwrite(cbind(cell_id=rownames(atac@meta.data),atac@meta.data),map_meta_file)
  umap_red <- grep('umap',names(atac@reductions),value=TRUE)[1]
  fwrite(cbind(cell_id=rownames(Embeddings(atac,umap_red)),as.data.frame(Embeddings(atac,umap_red))),map_umap_file)
}

# KLF4/epithelial/mesenchymal accessibility supports an ATAC Fib_K14-like call;
# it is explicitly not described as RNA expression.
DefaultAssay(atac) <- 'RNA'
d <- GetAssayData(atac,slot='data')
epi <- intersect(c('Tp63','Krt5','Krt14','Epcam','Tacstd2'),rownames(d))
mes <- intersect(c('Vim','Pdgfra','Col1a1'),rownames(d))
apos <- function(g) if(length(g)) Matrix::colSums(d[g,,drop=FALSE]>0)>0 else rep(FALSE,ncol(d))
klf <- if ('Klf4' %in% rownames(d)) as.numeric(d['Klf4',])>0 else rep(FALSE,ncol(d))
atac$Fib_K14_ATAC_supported <- as.character(atac$predicted.sub_labels)=='Fib_K14' & atac$predicted.sub_labels.score>=0.5 & klf & apos(epi) & apos(mes)

msg('BUILD JOINT EMBEDDING')
if ('ref.pca' %in% names(atac@reductions)) {
  rp <- Embeddings(rna2,'pca')[,1:30,drop=FALSE]
  ap <- Embeddings(atac,'ref.pca')[,1:30,drop=FALSE]
  shared_geometry <- 'reference PCA plus projected ATAC PCA'
} else {
  ured <- grep('umap',names(atac@reductions),value=TRUE)[1]
  if(is.na(ured)) stop('No projected UMAP reduction found after MapQuery')
  rp <- Embeddings(rna2,'umap.reference')
  ap <- Embeddings(atac,ured)
  shared_geometry <- paste0('reference UMAP plus ATAC ',ured,' projection')
}
emb <- rbind(rp,ap)
jmeta <- as.data.frame(rbindlist(list(
  data.frame(cell_id=colnames(rna2),modality='scRNA-seq',rmeta[colnames(rna2),,drop=FALSE],check.names=FALSE),
  data.frame(cell_id=colnames(atac),modality='scATAC-seq',atac@meta.data[colnames(atac),,drop=FALSE],check.names=FALSE)
),fill=TRUE,use.names=TRUE))
rownames(jmeta) <- jmeta$cell_id; jmeta$cell_id <- NULL
rownames(jmeta) <- rownames(emb)
jmeta$provided_bulk <- ifelse(jmeta$modality=='scRNA-seq',as.character(jmeta$bulk_labels),as.character(jmeta$predicted.bulk_labels))
jmeta$provided_sub <- ifelse(jmeta$modality=='scRNA-seq',as.character(jmeta$sub_labels_reannotated),as.character(jmeta$predicted.sub_labels))

signatures <- list(
 Fibroblast=c('Col1a1','Col1a2','Pdgfra','Dcn','Lum','Col3a1'), Keratinocyte=c('Krt5','Krt1','Krt10','Krt17','Epcam','Krt19'),
 Endothelial=c('Pecam1','Kdr','Vwf','Emcn'), Pericyte=c('Rgs5','Cspg4','Pdgfrb','Acta2'),
 Myeloid=c('Ptprc','Lyz2','Cd68','Adgre1','Csf1r'), T_NK=c('Cd3d','Cd3e','Nkg7','Cd247'),
 B_cell=c('Cd79a','Ms4a1','Cd37'), Mast=c('Kit','Tpsab1','Cpa3'), Schwann=c('S100b','Sox10','Plp1'),
 Melanocyte=c('Pmel','Mitf','Tyr'), Muscle=c('Myod1','Myog','Des'))
score_mat <- matrix(NA_real_,nrow=nrow(jmeta),ncol=length(signatures),dimnames=list(rownames(jmeta),names(signatures)))
for(nm in names(signatures)){
  g <- intersect(signatures[[nm]],common); if(!length(g)) next
  score_mat[colnames(rna2),nm] <- Matrix::colMeans(GetAssayData(rna2,slot='data')[g,,drop=FALSE])
  score_mat[colnames(atac),nm] <- Matrix::colMeans(GetAssayData(atac,slot='data')[g,,drop=FALSE])
}
coords <- as.data.frame(emb); colnames(coords) <- paste0('UMAP_',seq_len(ncol(coords))); coords$cell_id <- rownames(coords)
outmeta <- cbind(data.frame(cell_id=rownames(jmeta)),jmeta)
fwrite(coords,file.path(OUT,'joint_umap_coordinates.csv'))
fwrite(coords,file.path(OUT,'joint_pca_coordinates.csv.gz'))
fwrite(outmeta,file.path(OUT,'joint_base_metadata.csv.gz'))
fwrite(cbind(cell_id=rownames(score_mat),as.data.frame(score_mat)),file.path(OUT,'joint_signature_scores.csv.gz'))
writeLines(if(is.null(anchors)) 'Restored from mapped metadata/UMAP checkpoint' else capture.output(anchors),file.path(OUT,'RNA_ATAC_CCA_anchor_summary.txt'))

summary <- list(seed=20260811,n_rna=ncol(rna2),n_atac=ncol(atac),n_joint=nrow(emb),n_common_genes=length(common),
                integration='Seurat anchor-based CCA with MapQuery; unpaired modalities',
                shared_geometry=shared_geometry,
                atac_primary='published-style uncorrected IterativeLSI retained in ArchR child',
                batch_correction='Harmony is sensitivity-only because library/sample is confounded with donor/time',
                n_atac_FibK14_supported=sum(atac$Fib_K14_ATAC_supported))
writeLines(jsonlite::toJSON(summary,pretty=TRUE,auto_unbox=TRUE),file.path(OUT,'stage2_R_summary.json'))
writeLines('PASS',file.path(OUT,'STAGE2_R_DONE'))
msg('STAGE2_R_DONE',jsonlite::toJSON(summary,auto_unbox=TRUE))
