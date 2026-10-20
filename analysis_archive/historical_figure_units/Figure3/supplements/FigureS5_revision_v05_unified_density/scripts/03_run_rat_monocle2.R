#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle)
  library(Matrix)
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 03_run_rat_monocle2.R <FigureS5 package>")
# Keep an ASCII junction path unresolved because this Windows R build cannot
# reliably reopen paths containing CJK characters after normalizePath().
root <- args[[1]]
tables <- file.path(root, "tables")
figures <- file.path(root, "figures", "panels")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

counts_cells_genes <- readMM(file.path(tables, "S5E_monocle2_counts_cells_by_genes.mtx"))
genes <- read.csv(file.path(tables, "S5E_monocle2_genes.csv"), check.names = FALSE)
pheno <- read.csv(file.path(tables, "S5E_monocle2_phenodata.csv"), check.names = FALSE)
counts <- t(counts_cells_genes)
rownames(counts) <- genes$gene
colnames(counts) <- pheno$cell_id
rownames(pheno) <- pheno$cell_id
pd <- new("AnnotatedDataFrame", data = pheno)
fd <- new("AnnotatedDataFrame", data = data.frame(
  gene_short_name = genes$gene, row.names = genes$gene
))

cds <- newCellDataSet(
  counts, phenoData = pd, featureData = fd,
  expressionFamily = negbinomial.size(), lowerDetectionLimit = 0.5
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, relative_expr = TRUE)
cds <- detectGenes(cds, min_expr = 0.1)
disp <- dispersionTable(cds)
ordering <- subset(
  disp, mean_expression >= 0.1 & dispersion_empirical >= dispersion_fit
)$gene_id
if (length(ordering) > 2000) {
  ordering <- head(ordering[order(-disp$dispersion_empirical[
    match(ordering, disp$gene_id)
  ])], 2000)
}
cds <- setOrderingFilter(cds, ordering)
cds <- reduceDimension(
  cds, max_components = 2, method = "DDRTree", norm_method = "log",
  pseudo_expr = 1, relative_expr = TRUE, verbose = TRUE
)

# Monocle2 2.30.1 uses igraph::nei(), removed in igraph >=2.1.
if (utils::packageVersion("igraph") >= "2.1.0") {
  monocle_ns <- asNamespace("monocle")
  project2MST_compat <- function(cds, Projection_Method) {
    dp_mst <- minSpanningTree(cds)
    Z <- reducedDimS(cds); Y <- reducedDimK(cds)
    cds <- findNearestPointOnMST(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex_names <- colnames(Y)[closest_vertex]
    closest_vertex_df <- as.matrix(closest_vertex)
    row.names(closest_vertex_df) <- row.names(closest_vertex)
    tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
    if (!is.function(Projection_Method)) {
      P <- Y[, closest_vertex]
    } else {
      P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
      for (i in seq_along(closest_vertex)) {
        neighbors_i <- names(igraph::neighbors(
          dp_mst, closest_vertex_names[i], mode = "all"
        ))
        projection <- NULL; distance <- NULL; Z_i <- Z[, i]
        for (neighbor in neighbors_i) {
          if (closest_vertex_names[i] %in% tip_leaves) {
            tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
          } else {
            tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
          }
          projection <- rbind(projection, tmp)
          distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
        }
        if (!methods::is(projection, "matrix")) projection <- as.matrix(projection)
        P[, i] <- projection[which(distance == min(distance))[1], ]
      }
    }
    colnames(P) <- colnames(Z)
    dp <- as.matrix(stats::dist(t(P)))
    min_dist <- min(dp[dp != 0]); dp <- dp + min_dist; diag(dp) <- 0
    cellPairwiseDistances(cds) <- dp
    gp <- igraph::graph_from_adjacency_matrix(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- igraph::mst(gp)
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
    cds
  }
  environment(project2MST_compat) <- monocle_ns
  unlockBinding("project2MST", monocle_ns)
  assign("project2MST", project2MST_compat, envir = monocle_ns)
  lockBinding("project2MST", monocle_ns)
}
cds <- orderCells(cds)

coords <- as.data.frame(t(reducedDimS(cds)))
colnames(coords) <- c("Component_1", "Component_2")
out <- cbind(pData(cds), coords)
write.csv(out, file.path(tables, "S5E_monocle2_coordinates_pseudotime.csv"),
          row.names = FALSE)
saveRDS(cds, file.path(tables, "S5E_monocle2_cell_data_set.rds"))

state_colors <- c(
  Fib_K14 = "#E31A1C", SAC_Progenitor = "#0099CC",
  KC_Basal = "#008B45", KC_Basal_Mig = "#FF8C00",
  KC_Spinous_Mig = "#00008B", KC_Spinous = "#00A000"
)
time_levels <- c(
  "Normal", "0dpi", "0.5dpi", "1dpi", "2dpi", "3dpi", "5dpi",
  "7dpi", "10dpi", "11dpi", "15dpi", "19dpi"
)
out$time_point <- factor(out$time_point, levels = time_levels)
time_colors <- setNames(
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(
    length(time_levels)
  ), time_levels
)
out$display_state <- ifelse(out$sub_labels == "Fib_K14", "Fib_K14-like", out$sub_labels)
names(state_colors)[names(state_colors) == "Fib_K14"] <- "Fib_K14-like"

theme_pub <- theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 7), axis.text = element_text(size = 6),
    legend.title = element_text(size = 7), legend.text = element_text(size = 6)
  )
p1 <- ggplot(out, aes(Component_1, Component_2, color = Pseudotime)) +
  geom_point(size = 0.34, alpha = 0.85) +
  scale_color_viridis_c(option = "plasma") + coord_equal() +
  labs(title = "Rat Monocle2 pseudotime", color = "Pseudotime") + theme_pub
p2 <- ggplot(out, aes(Component_1, Component_2, color = display_state)) +
  geom_point(size = 0.34, alpha = 0.85) +
  scale_color_manual(values = state_colors) + coord_equal() +
  labs(title = "Cell-state mapping", color = NULL) + theme_pub
p3 <- ggplot(out, aes(Component_1, Component_2, color = time_point)) +
  geom_point(size = 0.34, alpha = 0.85) +
  scale_color_manual(values = time_colors, na.value = "grey70", drop = FALSE) +
  coord_equal() + labs(title = "Observed wound time", color = NULL) +
  guides(color = guide_legend(ncol = 4, byrow = TRUE,
                              override.aes = list(size = 1.4))) +
  theme_pub + theme(legend.position = "bottom")

combined <- arrangeGrob(p1, p2, p3, nrow = 1)
ggsave(file.path(figures, "FigureS5_E_rat_monocle2_trajectory.pdf"), combined,
       width = 7.08, height = 2.35, units = "in", device = cairo_pdf)
ggsave(file.path(figures, "FigureS5_E_rat_monocle2_trajectory.png"), combined,
       width = 7.08, height = 2.35, units = "in", dpi = 600, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(figures, "FigureS5_E_rat_monocle2_trajectory.svg"), combined,
         width = 7.08, height = 2.35, units = "in", device = "svg")
}
