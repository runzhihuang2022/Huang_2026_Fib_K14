suppressPackageStartupMessages({
  library(ArchR); library(BSgenome.Rnorvegicus.UCSC.rn7)
  library(grid); library(gtable)
})
addArchRThreads(12)
P <- 'data/external/local_7dd29a255152'
out <- file.path(P,'analysis','06_revisualized_peak_tracks_all_Fib')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
proj <- readRDS(file.path(P,'analysis','03_peaks_motifs','rat_burn_scatac_stage3_peaks_motifs.rds'))
corr <- read.csv(file.path(P,'analysis','04_fibk14_klf4','fibk14_klf4_target_peak_correlations_full.csv'))

# Figure contract: one gene per compact stacked coverage/peak panel.
width_mm <- 100; height_mm <- 95; dpi <- 300
groups <- c('Fib_K14','Fib_Inflama','Fib_Cd45')
pal <- c(Fib_K14='#D55E00',Fib_Inflama='#56B4E9',Fib_Cd45='#009E73',
         KC_Spinous='#0072B2',KC_Spinous_Mat='#CC79A7',KC_Spinous_Mig='#E69F00',
         SAC_Ductal='#8FBC8F',SAC_HF_Matrix='#6A5ACD',M1_Macrophage='#7F7F7F',M2_Macrophage='#00A6A6')
genes <- c('Klf4','Krt5','Krt15','Krt17','Krt1','Mmp3','Cd34','Grhl3','Col1a1','Areg','Krt10')

arialize <- function(g){
  if (inherits(g,'text')) {
    lab <- as.character(g$label)
    if (length(lab) && (any(grepl('Coverage|Norm\\. ATAC Signal',lab)) ||
                        all(grepl('^[0-9]{7,}$',lab)))) g$label <- ''
    if (is.null(g$gp)) g$gp <- gpar()
    g$gp$fontfamily <- 'Arial'
    fs <- g$gp$fontsize
    if (is.null(fs) || length(fs)==0 || any(!is.finite(fs))) fs <- 9
    g$gp$fontsize <- pmin(11,pmax(8,fs))
  }
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs,arialize)
  if (!is.null(g$children)) g$children <- do.call(gList,lapply(g$children,arialize))
  g
}

peak_gr <- function(row){
  z <- strcapture('([^:]+):([0-9]+)-([0-9]+)',row$peak,
                  proto=list(chr=character(),start=integer(),end=integer()))
  GRanges(z$chr,IRanges(z$start,z$end))
}

manifest <- list()
for (gene in genes){
  rr <- corr[corr$target_gene==gene,]
  if(nrow(rr)>0){ rr <- rr[order(rr$FDR),]; best <- rr[1,]; hi <- peak_gr(best)
    ttl <- sprintf('%s  |  KLF4–peak r = %.2f, FDR = %.1e',gene,best$partial_spearman,best$FDR)
    ttl <- sprintf('%s | KLF4-peak r=%.2f | FDR=%.1e',gene,best$partial_spearman,best$FDR)
    direction <- ifelse(best$partial_spearman>=0,'positive','negative')
  } else { best <- NULL; hi <- NULL; ttl <- gene; direction <- 'not_tested' }
  ttl <- gene
  pp <- plotBrowserTrack(ArchRProj=proj,groupBy='cell_type_cluster_final',useGroups=groups,
    geneSymbol=gene,upstream=25000,downstream=25000,
    plotSummary=c('bulkTrack','featureTrack','geneTrack'),sizes=c(10,1.8,3.2),
    maxCells=500,minCells=25,normMethod='ReadsInTSS',highlight=hi,
    highlightFill='#F4A6A6',pal=pal,baseSize=8,facetbaseSize=8,
    borderWidth=0.45,tickWidth=0.45,title=ttl,threads=12,verbose=FALSE)
  g <- arialize(pp[[1]])
  stem <- file.path(out,paste0('PeakTrack_',gene))
  cairo_pdf(paste0(stem,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,
            family='Arial',onefile=TRUE)
  grid.newpage(); grid.draw(g); dev.off()
  png(paste0(stem,'.png'),width=width_mm,height=height_mm,units='mm',res=dpi,
      type='cairo-png',bg='white')
  grid.newpage(); grid.draw(g); dev.off()
  manifest[[gene]] <- data.frame(gene=gene,width_mm=width_mm,height_mm=height_mm,dpi=dpi,
    highlighted_peak=if(is.null(best)) NA else best$peak,
    partial_spearman=if(is.null(best)) NA else best$partial_spearman,
    FDR=if(is.null(best)) NA else best$FDR,direction=direction)
}
write.csv(do.call(rbind,manifest),file.path(out,'peak_track_manifest.csv'),row.names=FALSE)
ct <- as.character(proj@cellColData[['cell_type_cluster_final']])
ctab <- as.data.frame(table(ct)); ctab <- ctab[ctab$ct %in% groups,]
write.csv(ctab,
          file.path(out,'source_cell_counts.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('REVIS_DONE',length(genes),'figures\n')
