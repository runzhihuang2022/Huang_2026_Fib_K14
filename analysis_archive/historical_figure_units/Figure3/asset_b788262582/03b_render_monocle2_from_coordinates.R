#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 03b_render_monocle2_from_coordinates.R <FigureS4 package>")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
tables <- file.path(root, "tables")
figures <- file.path(root, "figures", "panels")
out <- read.csv(
  file.path(tables, "S4E_monocle2_coordinates_pseudotime.csv"),
  check.names = FALSE
)

state_colors <- c(
  Fib_K14 = "#E31A1C", SAC_SG_Progenitor = "#0099CC",
  KC_Basal = "#008B45", KC_Basal_Mig = "#FF8C00",
  KC_Spinous_Mig = "#00008B", KC_Spinous = "#00A000"
)
time_levels <- c(
  "Normal", "5dpb", "6dpb", "7dpb", "10dpb", "11dpb", "12dpb",
  "19dpb", "22dpb", "26dpb", "1mph", "2mph"
)
out$time_point <- factor(out$time_point, levels = time_levels)
time_colors <- setNames(
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(
    length(time_levels)
  ),
  time_levels
)

theme_pub <- theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  )

p1 <- ggplot(out, aes(Component_1, Component_2, color = Pseudotime)) +
  geom_point(size = 0.35, alpha = 0.85) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() + labs(title = "Monocle2 pseudotime", color = "Pseudotime") +
  theme_pub
p2 <- ggplot(out, aes(Component_1, Component_2, color = sub_labels)) +
  geom_point(size = 0.35, alpha = 0.85) +
  scale_color_manual(values = state_colors) +
  coord_equal() + labs(title = "Cell-state mapping", color = NULL) +
  theme_pub
p3 <- ggplot(out, aes(Component_1, Component_2, color = time_point)) +
  geom_point(size = 0.35, alpha = 0.85) +
  scale_color_manual(values = time_colors, na.value = "grey70", drop = FALSE) +
  coord_equal() + labs(title = "Observed healing time", color = NULL) +
  guides(color = guide_legend(
    ncol = 4, byrow = TRUE, override.aes = list(size = 1.4)
  )) +
  theme_pub +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(0.28, "cm"),
    legend.spacing.x = grid::unit(0.03, "cm")
  )

combined <- arrangeGrob(p1, p2, p3, nrow = 1)
ggsave(file.path(figures, "FigureS4_E_monocle2_trajectory.pdf"), combined,
       width = 7.08, height = 2.35, units = "in", device = cairo_pdf)
ggsave(file.path(figures, "FigureS4_E_monocle2_trajectory.png"), combined,
       width = 7.08, height = 2.35, units = "in", dpi = 600, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(figures, "FigureS4_E_monocle2_trajectory.svg"), combined,
         width = 7.08, height = 2.35, units = "in", device = "svg")
}
