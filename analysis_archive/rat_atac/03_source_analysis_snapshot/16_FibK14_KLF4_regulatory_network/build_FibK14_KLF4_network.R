suppressPackageStartupMessages({
  library(grid)
})

root <- "data/external/local_b7a9e07281f7"
in_dir <- file.path(root, "04_FibK14_KLF4")
val_dir <- file.path(root, "05_GSE193564_validation")
motif_dir <- file.path(root, "03_peaks_motifs")
out_dir <- file.path(root, "16_FibK14_KLF4_regulatory_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

W_MM <- 120
H_MM <- 80
DPI <- 600

COL <- c(
  ink = "#202426",
  grey = "#697177",
  light_grey = "#D7DEE2",
  pale = "#F7F9FA",
  coral = "#D9504A",
  coral_light = "#F6C2BD",
  blue = "#236C9D",
  blue_light = "#BBD9EA",
  teal = "#258F86",
  teal_light = "#CBE7E2",
  orange = "#E4872E",
  orange_light = "#F4D4AE",
  amber = "#B98116",
  amber_light = "#F3E1B7"
)

num_cols <- c(
  "peak2gene_cor", "peak2gene_FDR", "detection_fraction",
  "partial_spearman", "p_value", "FDR", "sample_spearman",
  "sample_p", "n_samples", "sample_FDR"
)

x <- read.csv(
  file.path(in_dir, "fibk14_klf4_target_peak_correlations_full.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
for (nm in intersect(num_cols, names(x))) {
  x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
}

role_map <- data.frame(
  target_gene = c(
    "Krt1", "Krt5", "Krt10", "Krt15", "Krt17", "Col17a1",
    "Grhl3", "Areg", "Mmp3", "Tacstd2", "Col1a1", "Cd34"
  ),
  display_gene = c(
    "KRT1", "KRT5", "KRT10", "KRT15", "KRT17", "COL17A1",
    "GRHL3", "AREG", "MMP3", "TACSTD2", "COL1A1", "CD34"
  ),
  biological_role = c(
    "Differentiation/barrier",
    "Basal epithelial identity",
    "Differentiation/barrier",
    "Progenitor/barrier",
    "Wound epithelialization",
    "Epithelial adhesion",
    "Epithelial transcription",
    "EGFR signaling",
    "Wound remodeling",
    "Epithelial identity",
    "Fibrosis/ECM context",
    "Fibroblast state context"
  ),
  epithelial_class = c(
    "Differentiation", "Basal/progenitor", "Differentiation",
    "Basal/progenitor", "Wound repair", "Adhesion",
    "Transcription", "Signaling", "Remodeling", "Identity",
    "Context", "Context"
  ),
  stringsAsFactors = FALSE
)

summarize_gene <- function(z) {
  pos <- z[!is.na(z$FDR) & z$FDR < 0.05 & z$partial_spearman > 0, , drop = FALSE]
  neg <- z[!is.na(z$FDR) & z$FDR < 0.05 & z$partial_spearman < 0, , drop = FALSE]
  sample_pos <- pos[
    !is.na(pos$sample_p) & pos$sample_p < 0.05 &
      !is.na(pos$sample_spearman) & pos$sample_spearman > 0,
    , drop = FALSE
  ]
  best_pos <- if (nrow(pos)) pos[which.max(pos$partial_spearman), , drop = FALSE] else NULL
  data.frame(
    target_gene = z$target_gene[1],
    tested_peaks = nrow(z),
    positive_cell_FDR = nrow(pos),
    negative_cell_FDR = nrow(neg),
    max_positive_partial_rho = if (nrow(pos)) max(pos$partial_spearman) else NA_real_,
    min_positive_FDR = if (nrow(pos)) min(pos$FDR) else NA_real_,
    positive_sample_rawP = nrow(sample_pos),
    best_positive_sample_rho = if (nrow(pos)) max(pos$sample_spearman, na.rm = TRUE) else NA_real_,
    min_positive_sample_p = if (nrow(pos)) min(pos$sample_p, na.rm = TRUE) else NA_real_,
    best_positive_peak = if (nrow(pos)) best_pos$peak[1] else NA_character_,
    best_positive_peak2gene_cor = if (nrow(pos)) best_pos$peak2gene_cor[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

gene_summary <- do.call(rbind, lapply(split(x, x$target_gene), summarize_gene))
gene_summary <- merge(role_map, gene_summary, by = "target_gene", all.x = TRUE)

for (nm in c("tested_peaks", "positive_cell_FDR", "negative_cell_FDR", "positive_sample_rawP")) {
  gene_summary[[nm]][is.na(gene_summary[[nm]])] <- 0
}

gene_summary$evidence_class <- with(
  gene_summary,
  ifelse(
    positive_cell_FDR > negative_cell_FDR & positive_cell_FDR > 0,
    "positive-dominant candidate",
    ifelse(
      positive_cell_FDR > 0 & negative_cell_FDR > 0,
      "mixed-direction candidate",
      ifelse(
        positive_cell_FDR == 0 & negative_cell_FDR > 0,
        "negative-only candidate",
        "no cell-level FDR edge"
      )
    )
  )
)

gene_summary$sample_evidence <- ifelse(
  gene_summary$positive_sample_rawP > 0,
  "nominal positive sample support",
  "not sample-replicated"
)

gene_order <- role_map$target_gene
gene_summary$order <- match(gene_summary$target_gene, gene_order)
gene_summary <- gene_summary[order(gene_summary$order), ]

write.csv(
  gene_summary,
  file.path(out_dir, "FibK14_KLF4_gene_evidence_summary.csv"),
  row.names = FALSE
)

# Positive-dominant epithelialization candidate arm.
main_genes <- c("Krt5", "Krt15", "Krt17", "Krt1", "Mmp3")
main_summary <- gene_summary[match(main_genes, gene_summary$target_gene), ]

positive_edges_all <- x[
  x$target_gene %in% main_genes &
    !is.na(x$FDR) & x$FDR < 0.05 & x$partial_spearman > 0,
  ,
  drop = FALSE
]
positive_edges_all <- positive_edges_all[
  order(match(positive_edges_all$target_gene, main_genes),
        -positive_edges_all$partial_spearman),
]
write.csv(
  positive_edges_all,
  file.path(out_dir, "FibK14_KLF4_positive_peak_edges_all.csv"),
  row.names = FALSE
)

network_edges <- data.frame(
  from = "KLF4",
  to = main_summary$display_gene,
  edge_type = "KLF4 expression–peak accessibility association",
  positive_peaks = main_summary$positive_cell_FDR,
  negative_peaks = main_summary$negative_cell_FDR,
  max_partial_rho = main_summary$max_positive_partial_rho,
  min_FDR = main_summary$min_positive_FDR,
  sample_evidence = main_summary$sample_evidence,
  biological_role = main_summary$biological_role,
  stringsAsFactors = FALSE
)
write.csv(
  network_edges,
  file.path(out_dir, "FibK14_KLF4_network_edges.csv"),
  row.names = FALSE
)

network_nodes <- rbind(
  data.frame(
    node = "KLF4", node_type = "transcription factor",
    biological_role = "candidate regulator", stringsAsFactors = FALSE
  ),
  data.frame(
    node = main_summary$display_gene, node_type = "candidate target gene",
    biological_role = main_summary$biological_role, stringsAsFactors = FALSE
  ),
  data.frame(
    node = "Pro-epithelial candidate program", node_type = "functional outcome",
    biological_role = "annotation-supported outcome", stringsAsFactors = FALSE
  )
)
write.csv(
  network_nodes,
  file.path(out_dir, "FibK14_KLF4_network_nodes.csv"),
  row.names = FALSE
)

stage4 <- read.csv(file.path(in_dir, "stage4_statistical_summary.csv"))
stage4_value <- function(metric) {
  suppressWarnings(as.numeric(stage4$value[stage4$metric == metric][1]))
}

validation <- read.csv(
  file.path(val_dir, "GSE193564_FibK14_Klf4_epithelial_module_correlation.csv"),
  check.names = FALSE
)
external_rho <- as.numeric(validation$partial_spearman_sample_qc_controlled[1])
external_p <- as.numeric(validation$p_value[1])

chromvar <- read.csv(file.path(in_dir, "KLF4_chromVAR_by_cell_type.csv"))
fib_chromvar <- chromvar[
  chromvar$cell_type_cluster_final == "Fib_K14",
  ,
  drop = FALSE
]
fib_chromvar_mean <- as.numeric(fib_chromvar$KLF4_chromVAR_z.mean[1])

motifs <- read.csv(
  file.path(motif_dir, "marker_motifs_chromvar_by_cell_type.csv"),
  stringsAsFactors = FALSE
)
klf4_marker_fib <- motifs[
  grepl("^Klf4\\.", motifs$name) & motifs$cell_type == "Fib_K14",
  ,
  drop = FALSE
]

audit <- data.frame(
  metric = c(
    "Fib_K14 cells", "Samples", "Candidate target peaks",
    "Tested peaks", "Positive cell-level FDR<0.05",
    "Negative cell-level FDR<0.05", "Positive sample-level FDR<0.05",
    "Fib_K14 KLF4 chromVAR mean z", "Fib_K14 KLF4 marker-motif rows",
    "External Fib_K14-like Klf4-epithelial module partial rho",
    "External validation p"
  ),
  value = c(
    stage4_value("Fib_K14_cells"),
    stage4_value("samples"),
    stage4_value("candidate_target_peaks"),
    stage4_value("tested_peaks"),
    stage4_value("positive_FDR_lt_0.05"),
    stage4_value("negative_FDR_lt_0.05"),
    0,
    fib_chromvar_mean,
    nrow(klf4_marker_fib),
    external_rho,
    external_p
  ),
  stringsAsFactors = FALSE
)
write.csv(
  audit,
  file.path(out_dir, "FibK14_KLF4_network_evidence_audit.csv"),
  row.names = FALSE
)

gp_txt <- function(size = 8, face = "plain", col = COL[["ink"]]) {
  gpar(
    fontfamily = "Arial", fontsize = size,
    fontface = face, col = col
  )
}

txt <- function(label, x, y, size = 8, face = "plain",
                col = COL[["ink"]], just = "centre") {
  grid.text(label, x, y, just = just, gp = gp_txt(size, face, col))
}

line_gp <- function(col = COL[["ink"]], lwd = 0.7, lty = 1) {
  gpar(
    col = col, fill = col, lwd = lwd, lty = lty,
    lineend = "round", linejoin = "round"
  )
}

arrow_line <- function(x0, y0, x1, y1, col, lwd = 0.8, lty = 1) {
  grid.lines(
    c(x0, x1), c(y0, y1),
    arrow = arrow(type = "closed", length = unit(1.0, "mm")),
    gp = line_gp(col, lwd, lty)
  )
}

node_box <- function(x, y, w, h, fill, col, label,
                     label_col = COL[["ink"]], italic = FALSE) {
  grid.roundrect(
    x, y, width = w, height = h, r = unit(1.1, "mm"),
    gp = gpar(fill = fill, col = col, lwd = 0.75)
  )
  if (italic) {
    grid.text(
      parse(text = paste0("italic(", label, ")")),
      x, y, gp = gp_txt(8, "italic", label_col)
    )
  } else {
    txt(label, x, y, 8, "bold", label_col)
  }
}

role_color <- function(role) {
  if (grepl("Basal|Progenitor", role)) {
    c(fill = COL[["teal_light"]], edge = COL[["teal"]])
  } else if (grepl("Wound epithelialization", role)) {
    c(fill = COL[["coral_light"]], edge = COL[["coral"]])
  } else if (grepl("Differentiation", role)) {
    c(fill = COL[["blue_light"]], edge = COL[["blue"]])
  } else {
    c(fill = COL[["orange_light"]], edge = COL[["orange"]])
  }
}

draw_network <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  txt("a", 0.015, 0.960, 10, "bold", just = c("left", "centre"))
  txt(
    "Fib_K14 KLF4 epithelialization candidate network",
    0.060, 0.960, 9, "bold", just = c("left", "centre")
  )
  txt(
    "Association network; causal direction requires perturbation",
    0.060, 0.915, 8, "plain", COL[["grey"]],
    just = c("left", "centre")
  )

  # Evidence-side panel
  grid.roundrect(
    0.875, 0.525, width = 0.225, height = 0.700,
    r = unit(1.5, "mm"),
    gp = gpar(fill = COL[["pale"]], col = COL[["light_grey"]], lwd = 0.55)
  )
  txt("Evidence audit", 0.875, 0.835, 8, "bold", COL[["blue"]])
  evidence_lines <- c(
    sprintf("Cells  %d", stage4_value("Fib_K14_cells")),
    sprintf("Samples  %d", stage4_value("samples")),
    sprintf("Positive  %d", stage4_value("positive_FDR_lt_0.05")),
    sprintf("Negative  %d", stage4_value("negative_FDR_lt_0.05")),
    "Sample-FDR  0/124",
    sprintf("chromVAR  %.2f", fib_chromvar_mean)
  )
  evidence_cols <- c(
    COL[["ink"]], COL[["ink"]], COL[["coral"]],
    COL[["blue"]], COL[["blue"]], COL[["grey"]]
  )
  for (i in seq_along(evidence_lines)) {
    txt(
      evidence_lines[i], 0.785, 0.770 - (i - 1) * 0.078,
      8, if (i %in% c(3, 5)) "bold" else "plain",
      evidence_cols[i], just = c("left", "centre")
    )
  }
  txt(
    sprintf("External RNA\nr = %.3f", external_rho),
    0.875, 0.285, 8, "plain", COL[["blue"]]
  )
  txt(
    "Candidate network",
    0.875, 0.205, 8, "bold", COL[["amber"]]
  )

  # Network area
  grid.circle(
    0.105, 0.545, r = unit(5.5, "mm"),
    gp = gpar(fill = COL[["coral_light"]], col = COL[["coral"]], lwd = 1.0)
  )
  txt("KLF4", 0.105, 0.545, 10, "bold", COL[["coral"]])
  txt("Fib_K14", 0.105, 0.460, 8, "bold", COL[["teal"]])

  gene_y <- c(0.785, 0.665, 0.545, 0.425, 0.305)
  gene_x <- 0.355
  outcome_x <- 0.650

  node_box(
    outcome_x, 0.545, 0.180, 0.145,
    COL[["teal_light"]], COL[["teal"]],
    "Epithelialization\ncandidate", COL[["teal"]]
  )

  for (i in seq_len(nrow(main_summary))) {
    rr <- main_summary[i, ]
    cc <- role_color(rr$biological_role)
    grid.roundrect(
      gene_x, gene_y[i], width = 0.200, height = 0.090,
      r = unit(1.1, "mm"),
      gp = gpar(fill = cc[["fill"]], col = cc[["edge"]], lwd = 0.75)
    )
    grid.text(
      parse(text = paste0("italic(", rr$display_gene, ")")),
      gene_x, gene_y[i] + 0.014, gp = gp_txt(8, "italic")
    )
    txt(
      sprintf("%d peaks; r = %.3f",
              rr$positive_cell_FDR, rr$max_positive_partial_rho),
      gene_x, gene_y[i] - 0.020, 8, "plain", COL[["grey"]]
    )

    strength <- 0.65 + 4.5 * rr$max_positive_partial_rho
    arrow_line(
      0.155, 0.545, gene_x - 0.105, gene_y[i],
      COL[["coral"]], strength
    )

    # Dotted functional-annotation link to the epithelialization outcome.
    grid.lines(
      c(gene_x + 0.105, outcome_x - 0.095),
      c(gene_y[i], 0.545),
      gp = line_gp(COL[["teal"]], 0.55, 2)
    )
  }

  # Orthogonal evidence and legend
  grid.roundrect(
    0.370, 0.125, 0.625, 0.105,
    r = unit(1.2, "mm"),
    gp = gpar(fill = "#FFFDF8", col = COL[["light_grey"]], lwd = 0.55)
  )
  txt(
    "Solid coral: positive cell-level FDR edge",
    0.075, 0.145, 8, "plain", COL[["coral"]],
    just = c("left", "centre")
  )
  txt(
    "Dotted teal: functional annotation, not causal evidence",
    0.075, 0.100, 8, "plain", COL[["teal"]],
    just = c("left", "centre")
  )
}

draw_balance <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  txt("b", 0.015, 0.960, 10, "bold", just = c("left", "centre"))
  txt(
    "Complete KLF4–target peak evidence in Fib_K14",
    0.060, 0.960, 9, "bold", just = c("left", "centre")
  )

  plot_data <- gene_summary[gene_summary$epithelial_class != "Context", ]
  plot_data <- plot_data[order(
    -(plot_data$positive_cell_FDR - plot_data$negative_cell_FDR),
    -plot_data$positive_cell_FDR
  ), ]

  max_count <- max(
    plot_data$positive_cell_FDR,
    plot_data$negative_cell_FDR,
    na.rm = TRUE
  )
  yvals <- seq(0.835, 0.220, length.out = nrow(plot_data))
  x0 <- 0.400
  scale <- 0.260 / max_count

  txt("Negative peaks", 0.265, 0.890, 8, "bold", COL[["blue"]])
  txt("Positive peaks", 0.535, 0.890, 8, "bold", COL[["coral"]])
  grid.lines(c(x0, x0), c(0.180, 0.860), gp = line_gp(COL[["grey"]], 0.45))

  for (i in seq_len(nrow(plot_data))) {
    rr <- plot_data[i, ]
    grid.text(
      parse(text = paste0("italic(", rr$display_gene, ")")),
      0.155, yvals[i], just = c("right","centre"),
      gp = gp_txt(8, "italic")
    )
    if (rr$negative_cell_FDR > 0) {
      grid.rect(
        x0 - rr$negative_cell_FDR * scale/2, yvals[i],
        width = rr$negative_cell_FDR * scale, height = 0.040,
        gp = gpar(fill = COL[["blue"]], col = NA)
      )
    }
    if (rr$positive_cell_FDR > 0) {
      grid.rect(
        x0 + rr$positive_cell_FDR * scale/2, yvals[i],
        width = rr$positive_cell_FDR * scale, height = 0.040,
        gp = gpar(fill = COL[["coral"]], col = NA)
      )
    }
    txt(
      sprintf("%d-", rr$negative_cell_FDR),
      0.205, yvals[i], 8, "plain", COL[["blue"]]
    )
    txt(
      sprintf("%d+", rr$positive_cell_FDR),
      0.595, yvals[i], 8, "plain", COL[["coral"]]
    )
    if (rr$target_gene %in% main_genes) {
      grid.circle(
        0.675, yvals[i], r = unit(0.75, "mm"),
        gp = gpar(fill = COL[["teal"]], col = NA)
      )
    }
  }

  txt(
    "Teal dot: positive-dominant gene retained in the main candidate network",
    0.090, 0.135, 8, "plain", COL[["teal"]],
    just = c("left", "centre")
  )
  txt(
    sprintf("External RNA module: partial r = %.3f, p = %.2g; opposite direction.",
            external_rho, external_p),
    0.090, 0.080, 8, "plain", COL[["blue"]],
    just = c("left", "centre")
  )
}

export_figure <- function(draw_fun, stem) {
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"),
    width = W_MM/25.4, height = H_MM/25.4,
    family = "Arial", onefile = TRUE
  )
  draw_fun()
  dev.off()

  grDevices::svg(
    paste0(stem, ".svg"),
    width = W_MM/25.4, height = H_MM/25.4,
    family = "Arial", onefile = TRUE
  )
  draw_fun()
  dev.off()

  grDevices::png(
    paste0(stem, ".png"),
    width = W_MM, height = H_MM, units = "mm",
    res = DPI, type = "cairo", bg = "white"
  )
  draw_fun()
  dev.off()

  grDevices::tiff(
    paste0(stem, ".tiff"),
    width = W_MM, height = H_MM, units = "mm",
    res = DPI, compression = "lzw", type = "cairo", bg = "white"
  )
  draw_fun()
  dev.off()
}

export_figure(
  draw_network,
  file.path(out_dir, "FibK14_KLF4_pro_epithelial_candidate_network")
)
export_figure(
  draw_balance,
  file.path(out_dir, "FibK14_KLF4_complete_evidence_balance")
)

report_lines <- c(
  "# Fib_K14 KLF4 regulatory-network analysis",
  "",
  "## Main result",
  "",
  paste0(
    "A positive-dominant candidate arm links KLF4 to ",
    paste(main_summary$display_gene, collapse = ", "),
    ". These genes were retained because positive cell-level FDR edges ",
    "outnumbered negative cell-level FDR edges and they map to epithelial ",
    "identity, wound keratinization, differentiation or remodeling."
  ),
  "",
  "## Evidence limits",
  "",
  sprintf(
    "- %d of %d tested peaks were positive at cell-level FDR<0.05, while %d were negative.",
    stage4_value("positive_FDR_lt_0.05"),
    stage4_value("tested_peaks"),
    stage4_value("negative_FDR_lt_0.05")
  ),
  "- No peak passed BH FDR<0.05 in the 19-sample pseudobulk sensitivity analysis.",
  sprintf(
    "- Fib_K14 mean KLF4 chromVAR z-score was %.3f, and KLF4 was not a Fib_K14 marker motif in the exported marker-motif table.",
    fib_chromvar_mean
  ),
  sprintf(
    "- In GSE193564 Fib_K14-like cells, Klf4 and the epithelialization expression module were negatively correlated after sample/QC control (r=%.3f, p=%.3g).",
    external_rho, external_p
  ),
  "",
  "## Interpretation",
  "",
  paste0(
    "The defensible conclusion is that KLF4 is associated with a target-specific, ",
    "state-dependent pro-epithelial candidate branch in Fib_K14. The current ",
    "multi-omics data do not establish that KLF4 globally promotes epithelialization. ",
    "Causal wording requires KLF4 perturbation plus orthogonal target-occupancy validation."
  ),
  "",
  "## Recommended causal validation",
  "",
  "1. Fib_K14-specific KLF4 knockdown/overexpression with biological replicates.",
  "2. Measure KRT5, KRT15, KRT17, KRT1 and MMP3 expression and target-peak accessibility.",
  "3. KLF4 CUT&Tag/ChIP-qPCR at the selected peaks.",
  "4. Functional epithelialization assays: KRT14 induction, epithelial sheet formation and wound closure.",
  "5. Mediation analysis only after perturbation confirms KLF4 -> peak -> gene direction."
)
writeLines(
  report_lines,
  file.path(out_dir, "FibK14_KLF4_network_interpretation.md"),
  useBytes = TRUE
)

message("FIBK14_KLF4_NETWORK_DONE")
