suppressPackageStartupMessages({
  library(grid)
  library(png)
})

out_dir <- "data/external/local_0c446d9e1891"
base_img <- png::readPNG(file.path(out_dir, "beautiful_base_no_text.png"))

W_MM <- 70
H_MM <- 59.4
DPI <- 600

COL <- c(
  ink = "#222222",
  coral = "#D94F47",
  teal = "#238F87",
  grey = "#6E7579"
)

txt <- function(label, x, y, size = 8, face = "plain",
                col = COL[["ink"]], just = "centre") {
  grid.text(
    label, x = x, y = y, just = just,
    gp = gpar(
      fontfamily = "Arial",
      fontsize = size,
      fontface = face,
      col = col
    )
  )
}

draw_beautiful <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.raster(
    base_img,
    x = 0.5, y = 0.5,
    width = 1, height = 1,
    interpolate = TRUE
  )
  # Remove the generated divider; whitespace provides the panel separation.
  grid.rect(
    x = 0.575, y = 0.5,
    width = 0.010, height = 1,
    gp = gpar(fill = "white", col = NA)
  )

  # Panel title
  txt("a", 0.017, 0.968, 10, "bold", just = c("left", "centre"))
  txt(
    "Rat burn multi-omics", 0.065, 0.968, 10, "bold",
    just = c("left", "centre")
  )

  # One-row chronological sampling design
  txt("Burn", 0.090, 0.825, 8, "bold")
  txt("Time points", 0.305, 0.825, 8, "bold")
  txt("Biopsy", 0.525, 0.825, 8, "bold")
  txt("0/12/24 h  |  3/7/11/15/19 d", 0.325, 0.615, 8)

  # Integrative analysis label
  txt(
    "Multi-omics analysis", 0.285, 0.535, 9, "bold",
    COL[["teal"]]
  )

  # Mechanistic readout
  txt("Fib_K14", 0.785, 0.965, 9, "bold")
  txt("KLF4–peak links", 0.785, 0.915, 8.5, "bold", COL[["coral"]])

  genes <- c("KRT14", "KRT5", "TACSTD2", "TP63", "EPCAM")
  gene_y <- c(0.735, 0.610, 0.485, 0.360, 0.235)

  # Clear the generated label-capsule column and use journal-style direct labels.
  grid.rect(
    x = 0.647, y = 0.455,
    width = 0.150, height = 0.620,
    gp = gpar(fill = "white", col = NA)
  )

  for (i in seq_along(genes)) {
    grid.text(
      parse(text = paste0("italic(", genes[i], ")")),
      x = 0.700, y = gene_y[i],
      just = c("right", "centre"),
      gp = gpar(
        fontfamily = "Arial",
        fontsize = 8,
        fontface = "italic",
        col = COL[["ink"]]
      )
    )
  }
}

stem <- file.path(out_dir, "Rat_burn_multiomics_beautiful")

grDevices::png(
  paste0(stem, ".png"),
  width = W_MM, height = H_MM, units = "mm",
  res = DPI, type = "cairo", bg = "white"
)
draw_beautiful()
dev.off()

grDevices::tiff(
  paste0(stem, ".tiff"),
  width = W_MM, height = H_MM, units = "mm",
  res = DPI, compression = "lzw", type = "cairo",
  bg = "white"
)
draw_beautiful()
dev.off()

message("BEAUTIFUL_RASTER_DONE")
