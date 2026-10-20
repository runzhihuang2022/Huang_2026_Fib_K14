.libPaths("data/external/local_835d389f6169")
suppressPackageStartupMessages({
  library(CellChat)
  library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript render_cellchat_default_chord_v26.R <rds_dir> <palette.tsv> <figure_dir> <table_dir>")
}
rds_dir <- args[[1]]
palette_file <- args[[2]]
figure_dir <- args[[3]]
table_dir <- args[[4]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

groups <- c("D0_7", "D15_28", "D8_14", "M1_2")
interfaces <- list(
  Endo = "^(Endo_|Pericyte$)",
  Immune = "^(B_|CD4_|CD8_|NK$|cDC|pDC$|LAM$|LGHS$|M1_|M2_|Neutrophil$|Mast$)",
  KC = "^KC_",
  SAC = "^SAC_"
)
fib_pattern <- "^Fib_"

palette <- read.table(
  palette_file, sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, comment.char = "", check.names = FALSE
)
subtype_colors <- setNames(substr(palette$color, 1, 7), palette$subtype)

open_pdf <- function(path) {
  cairo_pdf(path, width = 52.5 / 25.4, height = 52.5 / 25.4,
            family = "Arial", pointsize = 8)
  par(family = "Arial", ps = 8, mar = c(0.25, 0.25, 1.2, 0.25))
}
open_png <- function(path) {
  png(path, width = 1240, height = 1240, res = 600,
      type = "cairo", family = "Arial", pointsize = 8)
  par(family = "Arial", ps = 8, mar = c(0.25, 0.25, 1.2, 0.25))
}

render_default_chord <- function(object, key, colors, title, path) {
  draw <- function() {
    circos.clear()
    netVisual_chord_gene(
      object,
      slot.name = "net",
      net = key,
      color.use = colors,
      lab.cex = 0.75,
      small.gap = 1,
      big.gap = 8,
      annotationTrackHeight = 0.025,
      link.visible = TRUE,
      scale = FALSE,
      directional = 1,
      link.target.prop = TRUE,
      reduce = -1,
      transparency = 0.25,
      link.border = NA,
      title.name = NULL,
      show.legend = FALSE,
      thresh = 0.05
    )
    circos.clear()
  }
  open_pdf(paste0(path, ".pdf")); draw(); dev.off()
  open_png(paste0(path, ".png")); draw(); dev.off()
}

records <- list()
selection <- list()
for (group in groups) {
  object_path <- file.path(rds_dir, paste0("CellChat_", group, "_sub_labels.rds"))
  object <- readRDS(object_path)
  labels <- levels(object@idents)
  fib <- grep(fib_pattern, labels, value = TRUE)
  for (interface_name in names(interfaces)) {
    partner <- grep(interfaces[[interface_name]], labels, value = TRUE)
    forward <- tryCatch(
      subsetCommunication(object, sources.use = fib, targets.use = partner),
      error = function(e) data.frame()
    )
    reverse <- tryCatch(
      subsetCommunication(object, sources.use = partner, targets.use = fib),
      error = function(e) data.frame()
    )
    lr <- rbind(forward, reverse)
    stem <- paste0(group, "_Fib_", interface_name, "_CellChat_default_LR_chord")
    if (!nrow(lr)) {
      records[[length(records) + 1]] <- data.frame(
        group = group, interface = paste0("Fib_", interface_name),
        status = "NOT_ESTIMABLE", n_all_lr = 0, n_display_lr = 0,
        reason = "No significant ligand-receptor pair"
      )
      next
    }

    lr <- lr[order(-lr$prob, lr$pval, lr$source, lr$target, lr$ligand, lr$receptor), ]
    lr$key_id <- paste(lr$source, lr$ligand, lr$target, lr$receptor, sep = "|")
    key <- lr[!duplicated(lr$key_id), ][seq_len(min(8, nrow(lr))), ]
    write.table(
      lr, file.path(table_dir, paste0(group, "_Fib_", interface_name, "_all_LR.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(
      key, file.path(table_dir, paste0(group, "_Fib_", interface_name, "_displayed_LR.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    involved <- unique(c(as.character(key$source), as.character(key$target)))
    color_use <- unname(subtype_colors[involved])
    color_use[is.na(color_use)] <- "#BDBDBD"
    names(color_use) <- involved
    render_default_chord(
      object = object,
      key = key,
      colors = color_use,
      title = paste(gsub("_", " ", group), paste0("Fib-", interface_name)),
      path = file.path(figure_dir, paste0(stem, "_raw_default"))
    )
    selection[[length(selection) + 1]] <- transform(
      key,
      group = group,
      interface = paste0("Fib_", interface_name),
      display_rank = seq_len(nrow(key))
    )
    records[[length(records) + 1]] <- data.frame(
      group = group, interface = paste0("Fib_", interface_name),
      status = "PASS", n_all_lr = nrow(lr), n_display_lr = nrow(key),
      reason = ""
    )
  }
}

record_table <- do.call(rbind, records)
write.table(record_table, file.path(table_dir, "CELLCHAT_DEFAULT_CHORD_INSTANCE_AUDIT.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (length(selection)) {
  write.table(do.call(rbind, selection), file.path(table_dir, "ALL_DISPLAYED_DEFAULT_CHORD_LR.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

writeLines(
  c(
    paste0("status=", ifelse(all(record_table$status %in% c("PASS", "NOT_ESTIMABLE")), "PASS_WITH_EXPLICIT_BLOCKS", "FAIL")),
    "engine=CellChat::netVisual_chord_gene",
    "canvas_mm=52.5x52.5",
    "font=Arial",
    "font_range_pt=6-8",
    "selection=top 8 unique significant LR pairs by probability then p-value; every estimable interface has at least 8",
    "labels=CellChat default cell subtype plus ligand/receptor gene labels",
    "colors=locked subtype UMAP palette",
    "interpretation=predicted communication; not causal validation"
  ),
  file.path(table_dir, "CELLCHAT_DEFAULT_CHORD_V26_AUDIT.txt")
)
