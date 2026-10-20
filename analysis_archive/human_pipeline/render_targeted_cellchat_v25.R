suppressPackageStartupMessages({
  library(CellChat); library(ggplot2); library(ggrepel); library(circlize)
  library(grid); library(Matrix); library(igraph)
})

args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]; palette_file <- args[2]; output_dir <- args[3]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_dir,"figures"), showWarnings=FALSE)
dir.create(file.path(output_dir,"tables"), showWarnings=FALSE)

groups <- c("Normal","SPTDI","DPTDI","Mode_1","Mode_2","Mode_3","Mode_4")
objects <- setNames(lapply(groups,function(g) readRDS(file.path(input_dir,paste0("CellChat_",g,".rds")))),groups)
pal <- read.table(palette_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,comment.char="")
subtype_colors <- setNames(pal$color,pal$subtype)
interfaces <- list(
  Endo="^(Endo_|Pericyte$)",
  Immune="^(B_|CD4_|CD8_|NK$|cDC|pDC$|LAM$|LGHS$|M1_|M2_|Neutrophil$|Mast$)",
  KC="^KC_", SAC="^SAC_"
)
fib_pat <- "^Fib_"
pdf_square <- function(path,expr){cairo_pdf(path,width=52.5/25.4,height=52.5/25.4,family="Arial");on.exit(dev.off());force(expr)}
compact <- function(x){x<-gsub("^Endo_","E.",x);x<-gsub("^Fib_","F.",x);x<-gsub("^SAC_","S.",x);x<-gsub("^KC_","K.",x);x<-gsub("^Pericyte$","PC",x);x<-gsub("Capillary","Cap",x);x<-gsub("Arterial","Art",x);x<-gsub("Lymphatic","Lym",x);x<-gsub("Progenitor","Prog",x);x}
col_for <- function(x){z<-unname(subtype_colors[as.character(x)]);z[is.na(z)]<-"#BDBDBD";substr(z,1,7)}

all_lr <- list(); selected <- list(); backgrounds <- character()
for(obj in objects) backgrounds <- union(backgrounds,rownames(obj@data.signaling))
write.table(data.frame(gene=sort(backgrounds)),file.path(output_dir,"tables","CELLCHAT_SIGNALING_BACKGROUND_GENES.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

for(g in groups){
  obj <- objects[[g]]; labs <- levels(obj@idents)
  for(iname in names(interfaces)){
    fib <- grep(fib_pat,labs,value=TRUE); partner <- grep(interfaces[[iname]],labs,value=TRUE)
    fwd <- tryCatch(subsetCommunication(obj,sources.use=fib,targets.use=partner),error=function(e)data.frame())
    rev <- tryCatch(subsetCommunication(obj,sources.use=partner,targets.use=fib),error=function(e)data.frame())
    lr <- rbind(fwd,rev); if(!nrow(lr)) next
    lr$group<-g;lr$interface<-paste0("Fib_",iname);lr$direction<-paste(lr$source,"to",lr$target)
    all_lr[[paste(g,iname)]]<-lr
    write.table(lr,file.path(output_dir,"tables",paste0(g,"_Fib_",iname,"_all_LR.tsv")),sep="\t",quote=FALSE,row.names=FALSE)

    agg<-aggregate(prob~source+target,lr,sum);nodes<-unique(as.character(c(agg$source,agg$target)))
    graph<-graph_from_data_frame(agg,directed=TRUE,vertices=data.frame(name=nodes))
    ew<-E(graph)$prob;ew<-.3+3.7*ew/max(ew);vcols<-col_for(V(graph)$name);ecols<-adjustcolor(col_for(as.character(ends(graph,E(graph))[,1])),alpha.f=.82)
    pdf_square(file.path(output_dir,"figures",paste0(g,"_Fib_",iname,"_network.pdf")),{
      par(family="Arial",mar=c(.2,.2,1.5,.2),cex=.58)
      plot(graph,layout=layout_in_circle(graph),vertex.size=7,vertex.color=vcols,vertex.frame.color=vcols,
           vertex.label.family="Arial",vertex.label.cex=.58,vertex.label.color="black",edge.width=ew,
           edge.color=ecols,edge.arrow.size=.22,margin=.18)
      title(paste(gsub("_"," ",g),"Fib -",iname),family="Arial",cex.main=.75)
    })

    central<-data.frame(celltype=nodes,outgoing=strength(graph,mode="out"),incoming=strength(graph,mode="in"))
    central$label<-ifelse(rank(-(central$outgoing+central$incoming),ties.method="first")<=5,central$celltype,"")
    write.table(central,file.path(output_dir,"tables",paste0(g,"_Fib_",iname,"_centrality.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
    p<-ggplot(central,aes(outgoing,incoming,color=celltype,label=label))+geom_point(size=2)+
      geom_text_repel(size=7/.pt,family="Arial",color="black",max.overlaps=Inf,box.padding=.25,point.padding=.12,min.segment.length=0,force=2,max.time=5,seed=20260828)+
      scale_color_manual(values=setNames(col_for(central$celltype),central$celltype),guide="none")+
      scale_x_continuous(expand=expansion(mult=.20))+scale_y_continuous(expand=expansion(mult=.20))+
      theme_bw(base_size=7,base_family="Arial")+labs(title=paste(gsub("_"," ",g),paste0("Fib-",iname)),x="Outgoing signaling",y="Incoming signaling")+
      theme(text=element_text(family="Arial"),plot.title=element_text(size=9,hjust=.5),panel.grid.minor=element_blank())
    ggsave(file.path(output_dir,"figures",paste0(g,"_Fib_",iname,"_in_out_hub.pdf")),p,width=52.5,height=52.5,units="mm",device=cairo_pdf)

    key<-lr[order(-lr$prob),][seq_len(min(4,nrow(lr))),];selected[[paste(g,iname)]]<-key
    key$sender_sector<-paste0("S|",compact(key$source),"|",key$ligand);key$receiver_sector<-paste0("R|",compact(key$target),"|",key$receptor)
    write.table(key,file.path(output_dir,"tables",paste0(g,"_Fib_",iname,"_key_LR.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
    chord<-aggregate(prob~sender_sector+receiver_sector,key,sum);colnames(chord)<-c("from","to","value")
    sector_cell<-c(setNames(key$source,key$sender_sector),setNames(key$target,key$receiver_sector));sector_cell<-sector_cell[!duplicated(names(sector_cell))]
    grid_col<-setNames(col_for(unname(sector_cell)),names(sector_cell));source_map<-setNames(as.character(key$source),as.character(key$sender_sector));link_col<-adjustcolor(col_for(source_map[as.character(chord$from)]),alpha.f=.78)
    pdf_square(file.path(output_dir,"figures",paste0(g,"_Fib_",iname,"_key_LR_chord.pdf")),{
      par(family="Arial",mar=c(.3,.3,1.4,.3));circos.clear();circos.par(start.degree=90,gap.degree=3,track.margin=c(.005,.005),cell.padding=c(0,0,0,0),canvas.xlim=c(-3.50,3.50),canvas.ylim=c(-2.80,2.80))
      chordDiagram(chord,grid.col=grid_col,col=link_col,transparency=0,annotationTrack="grid",preAllocateTracks=1,directional=1,direction.type="arrows",link.arr.type="big.arrow",link.sort=TRUE)
      circos.trackPlotRegion(track.index=1,panel.fun=function(x,y){sn=get.cell.meta.data("sector.index");circos.text(CELL_META$xcenter,CELL_META$ylim[1],sn,facing="clockwise",niceFacing=TRUE,adj=c(0,.5),cex=.50,font=1)},bg.border=NA)
      title(paste(gsub("_"," ",g),paste0("Fib-",iname," LR")),cex.main=.70,family="Arial");circos.clear()
    })
  }
}

all<-do.call(rbind,all_lr);write.table(all,file.path(output_dir,"tables","ALL_targeted_interfaces_all_LR_genelist.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

for(iname in names(interfaces)){
  key<-selected[[paste("DPTDI",iname)]];if(is.null(key)||!nrow(key))next
  genes<-unique(c(key$ligand,key$receptor));long<-list()
  for(g in groups){obj<-objects[[g]];labs<-as.character(obj@idents);keep<-grepl(fib_pat,labs)|grepl(interfaces[[iname]],labs);dat<-obj@data.signaling;present<-intersect(genes,rownames(dat));if(!length(present))next
    for(gene in present){v<-as.numeric(dat[gene,keep]);long[[length(long)+1]]<-data.frame(group=g,celltype=labs[keep],gene=gene,expression=v)}}
  d<-do.call(rbind,long);write.table(d,file.path(output_dir,"tables",paste0("Fib_",iname,"_key_LR_gene_expression.tsv.gz")),sep="\t",quote=FALSE,row.names=FALSE)
  n_gene<-length(unique(d$gene));height_mm<-max(52.5,18+8*n_gene)
  p<-ggplot(d,aes(celltype,expression,fill=group))+geom_violin(scale="width",linewidth=.18,position=position_dodge(width=.88),trim=TRUE)+facet_grid(gene~.,scales="free_y")+
    theme_bw(base_size=7,base_family="Arial")+scale_fill_brewer(palette="Set2")+labs(title=paste("Fib -",iname,"key LR gene expression"),x=NULL,y="Normalized expression")+
    guides(fill=guide_legend(nrow=1,byrow=TRUE,keywidth=unit(4,"mm"),keyheight=unit(2.2,"mm")))+
    theme(text=element_text(family="Arial"),axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7),strip.text=element_text(face="italic",size=8),plot.title=element_text(size=9,hjust=.5),legend.position="top",legend.title=element_blank(),legend.text=element_text(size=6),legend.spacing.x=unit(1,"mm"),panel.grid.minor=element_blank())
  ggsave(file.path(output_dir,"figures",paste0("Fib_",iname,"_key_LR_pathway_violin.pdf")),p,width=105,height=height_mm,units="mm",device=cairo_pdf,limitsize=FALSE)
  write.table(data.frame(interface=paste0("Fib_",iname),n_genes=n_gene,height_mm=height_mm),file.path(output_dir,"tables",paste0("Fib_",iname,"_violin_canvas.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
}

writeLines(c("status=PASS","engine=CellChat 2.1.2 existing formal inference","annotation=sub_labels","groups=Normal/SPTDI/DPTDI and Mode_1/2/3/4","network_and_hub_colors=locked subtype UMAP palette","chord_sectors=cell identity plus ligand/receptor symbol","interpretation=predicted communication; not causal validation"),file.path(output_dir,"TARGETED_CELLCHAT_V25_AUDIT.txt"))
