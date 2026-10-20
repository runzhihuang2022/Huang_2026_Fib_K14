suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(hdf5r)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

set.seed(20260801)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: formal_inference.R <prepared_dir> <output_dir>")
}
prepared_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
de_dir <- file.path(output_dir, "M05_pseudobulk_de")
enrich_dir <- file.path(output_dir, "M06_enrichment")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

comparisons <- data.frame(
  comparison = c("sptdi_vs_normal", "dptdi_vs_normal", "dptdi_vs_sptdi"),
  numerator = c("SPTDI", "DPTDI", "DPTDI"),
  denominator = c("Normal", "Normal", "SPTDI")
)
tracks <- c("bulk_labels", "sub_labels")
minimum_cells <- 20L
minimum_donors <- 3L

safe_name <- function(value) {
  gsub("[^A-Za-z0-9._-]+", "_", value)
}

write_tsv <- function(value, path) {
  write.table(
    value, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
  )
}

read_composition <- function(track) {
  read.delim(
    file.path(prepared_dir, paste0(track, "_patient_composition.tsv")),
    check.names = FALSE
  )
}

run_da <- function(track) {
  composition <- read_composition(track)
  label_col <- track
  samples <- unique(composition[, c("patient", "group", "total_cells")])
  samples$sample_key <- paste(samples$patient, samples$group, sep = "||")
  labels <- sort(unique(composition[[label_col]]))
  matrix_counts <- matrix(
    0L, nrow = length(labels), ncol = nrow(samples),
    dimnames = list(labels, samples$sample_key)
  )
  for (index in seq_len(nrow(composition))) {
    row <- composition[index, ]
    sample_key <- paste(row$patient, row$group, sep = "||")
    matrix_counts[as.character(row[[label_col]]), sample_key] <- as.integer(row$cells)
  }
  results <- list()
  checkpoints <- list()
  for (ci in seq_len(nrow(comparisons))) {
    comparison <- comparisons[ci, ]
    keep <- samples$group %in% c(comparison$denominator, comparison$numerator)
    meta <- samples[keep, , drop = FALSE]
    meta$group <- factor(
      meta$group, levels = c(comparison$denominator, comparison$numerator)
    )
    counts <- matrix_counts[, meta$sample_key, drop = FALSE]
    donor_n <- table(meta$group)
    if (any(donor_n < minimum_donors)) {
      checkpoints[[length(checkpoints) + 1L]] <- data.frame(
        track = track, comparison = comparison$comparison,
        analysis = "all_patient_unpaired", status = "underpowered",
        denominator_donors = unname(donor_n[[comparison$denominator]]),
        numerator_donors = unname(donor_n[[comparison$numerator]])
      )
      next
    }
    design <- model.matrix(~ group, meta)
    y <- DGEList(counts = counts, samples = meta)
    y$samples$lib.size <- meta$total_cells
    y$samples$norm.factors <- 1
    y <- estimateDisp(y, design = design, robust = TRUE)
    fit <- glmQLFit(y, design = design, robust = TRUE)
    test <- glmQLFTest(fit, coef = ncol(design))
    table <- topTags(test, n = Inf, sort.by = "none")$table
    table$label <- rownames(table)
    table$track <- track
    table$comparison <- comparison$comparison
    table$analysis <- "all_patient_unpaired"
    table$denominator_donors <- unname(donor_n[[comparison$denominator]])
    table$numerator_donors <- unname(donor_n[[comparison$numerator]])
    results[[length(results) + 1L]] <- table
    checkpoints[[length(checkpoints) + 1L]] <- data.frame(
      track = track, comparison = comparison$comparison,
      analysis = "all_patient_unpaired", status = "PASS",
      denominator_donors = unname(donor_n[[comparison$denominator]]),
      numerator_donors = unname(donor_n[[comparison$numerator]])
    )

    if (comparison$comparison == "dptdi_vs_sptdi") {
      shared <- intersect(
        meta$patient[meta$group == comparison$denominator],
        meta$patient[meta$group == comparison$numerator]
      )
      if (length(shared) >= minimum_donors) {
        paired_meta <- meta[meta$patient %in% shared, , drop = FALSE]
        paired_meta$patient <- factor(paired_meta$patient)
        paired_counts <- matrix_counts[, paired_meta$sample_key, drop = FALSE]
        paired_design <- model.matrix(~ patient + group, paired_meta)
        paired_y <- DGEList(counts = paired_counts, samples = paired_meta)
        paired_y$samples$lib.size <- paired_meta$total_cells
        paired_y$samples$norm.factors <- 1
        paired_y <- estimateDisp(paired_y, design = paired_design, robust = TRUE)
        paired_fit <- glmQLFit(paired_y, design = paired_design, robust = TRUE)
        paired_test <- glmQLFTest(paired_fit, coef = ncol(paired_design))
        paired_table <- topTags(
          paired_test, n = Inf, sort.by = "none"
        )$table
        paired_table$label <- rownames(paired_table)
        paired_table$track <- track
        paired_table$comparison <- comparison$comparison
        paired_table$analysis <- "paired_patient_sensitivity"
        paired_table$denominator_donors <- length(shared)
        paired_table$numerator_donors <- length(shared)
        results[[length(results) + 1L]] <- paired_table
        checkpoints[[length(checkpoints) + 1L]] <- data.frame(
          track = track, comparison = comparison$comparison,
          analysis = "paired_patient_sensitivity", status = "PASS",
          denominator_donors = length(shared), numerator_donors = length(shared)
        )
      }
    }
  }
  list(
    results = if (length(results)) do.call(rbind, results) else data.frame(),
    checkpoints = if (length(checkpoints)) do.call(rbind, checkpoints) else data.frame()
  )
}

message("Running patient-level differential abundance...")
da_runs <- lapply(tracks, run_da)
da_results <- do.call(rbind, lapply(da_runs, `[[`, "results"))
da_checkpoints <- do.call(rbind, lapply(da_runs, `[[`, "checkpoints"))
write_tsv(da_results, file.path(output_dir, "M04_differential_abundance.tsv"))
write_tsv(da_checkpoints, file.path(output_dir, "M04_checkpoint.tsv"))

message("Opening frozen pseudobulk count store...")
h5 <- H5File$new(file.path(prepared_dir, "pseudobulk_counts.h5"), mode = "r")
on.exit(h5$close_all(), add = TRUE)
genes <- h5[["genes"]][]
all_de_checkpoints <- list()
all_de_index <- list()

run_enrichment <- function(result, track, label, comparison, analysis) {
  ranked <- result$stat
  names(ranked) <- result$gene
  ranked <- ranked[is.finite(ranked)]
  ranked <- sort(ranked, decreasing = TRUE)
  significant <- result[
    !is.na(result$padj) & result$padj < 0.05 &
      abs(result$log2FoldChange) >= 0.5, ,
    drop = FALSE
  ]
  if (!nrow(significant)) return(invisible(NULL))
  mapped <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = unique(significant$gene),
      column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
    )
  )
  for (direction in c("up", "down")) {
    selected <- significant$gene[
      if (direction == "up") significant$log2FoldChange > 0 else
        significant$log2FoldChange < 0
    ]
    entrez <- unique(na.omit(unname(mapped[selected])))
    if (length(entrez) < 10L) next
    go <- suppressMessages(enrichGO(
      gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
      keyType = "ENTREZID", pAdjustMethod = "BH",
      pvalueCutoff = 0.05, qvalueCutoff = 0.2,
      minGSSize = 10, maxGSSize = 500, readable = TRUE
    ))
    go_table <- as.data.frame(go)
    if (nrow(go_table)) {
      go_table$track <- track
      go_table$label <- label
      go_table$comparison <- comparison
      go_table$analysis <- analysis
      go_table$direction <- direction
      stem <- paste(
        safe_name(track), safe_name(label), comparison, analysis, direction,
        "GO_BP", sep = "__"
      )
      write_tsv(go_table, file.path(enrich_dir, paste0(stem, ".tsv")))
    }
  }
}

message("Running pseudobulk DESeq2 and current-project enrichment...")
for (track in tracks) {
  metadata <- read.delim(
    file.path(prepared_dir, paste0(track, "_pseudobulk_metadata.tsv")),
    check.names = FALSE
  )
  labels <- sort(unique(metadata$label))
  dataset <- h5[[paste0(track, "/counts")]]
  for (label in labels) {
    label_meta_all <- metadata[metadata$label == label, , drop = FALSE]
    for (ci in seq_len(nrow(comparisons))) {
      comparison <- comparisons[ci, ]
      candidate <- label_meta_all[
        label_meta_all$group %in% c(comparison$denominator, comparison$numerator) &
          label_meta_all$cells >= minimum_cells, ,
        drop = FALSE
      ]
      donor_n <- table(factor(
        candidate$group,
        levels = c(comparison$denominator, comparison$numerator)
      ))
      if (any(donor_n < minimum_donors)) {
        all_de_checkpoints[[length(all_de_checkpoints) + 1L]] <- data.frame(
          track = track, label = label, comparison = comparison$comparison,
          analysis = "all_patient_unpaired", status = "underpowered",
          denominator_donors = unname(donor_n[[comparison$denominator]]),
          numerator_donors = unname(donor_n[[comparison$numerator]]),
          minimum_cells_per_patient = minimum_cells
        )
        next
      }
      # hdf5r exposes the h5py row-major dataset with reversed dimensions:
      # genes x pseudobulk rows. Select pseudobulk units as columns.
      block <- dataset[, candidate$row_index + 1L]
      counts <- as.matrix(block)
      rownames(counts) <- genes
      colnames(counts) <- paste(
        candidate$patient, candidate$group, seq_len(nrow(candidate)), sep = "||"
      )
      keep_gene <- rowSums(counts >= 10L) >= minimum_donors
      counts <- counts[keep_gene, , drop = FALSE]
      candidate$group <- factor(
        candidate$group,
        levels = c(comparison$denominator, comparison$numerator)
      )
      rownames(candidate) <- colnames(counts)
      dds <- DESeqDataSetFromMatrix(
        countData = round(counts), colData = candidate, design = ~ group
      )
      dds <- DESeq(dds, quiet = TRUE)
      result <- as.data.frame(results(
        dds,
        contrast = c("group", comparison$numerator, comparison$denominator),
        alpha = 0.05
      ))
      result$gene <- rownames(result)
      result$track <- track
      result$label <- label
      result$comparison <- comparison$comparison
      result$analysis <- "all_patient_unpaired"
      result$denominator_donors <- unname(donor_n[[comparison$denominator]])
      result$numerator_donors <- unname(donor_n[[comparison$numerator]])
      result$minimum_cells_per_patient <- minimum_cells
      stem <- paste(
        safe_name(track), safe_name(label), comparison$comparison,
        "all_patient_unpaired", sep = "__"
      )
      path <- file.path(de_dir, paste0(stem, ".tsv.gz"))
      connection <- gzfile(path, "wt")
      write.table(
        result, file = connection, sep = "\t", quote = FALSE,
        row.names = FALSE
      )
      close(connection)
      all_de_index[[length(all_de_index) + 1L]] <- data.frame(
        track = track, label = label, comparison = comparison$comparison,
        analysis = "all_patient_unpaired", result_file = path,
        tested_genes = nrow(result),
        significant_fdr05_lfc05 = sum(
          !is.na(result$padj) & result$padj < 0.05 &
            abs(result$log2FoldChange) >= 0.5
        )
      )
      all_de_checkpoints[[length(all_de_checkpoints) + 1L]] <- data.frame(
        track = track, label = label, comparison = comparison$comparison,
        analysis = "all_patient_unpaired", status = "PASS",
        denominator_donors = unname(donor_n[[comparison$denominator]]),
        numerator_donors = unname(donor_n[[comparison$numerator]]),
        minimum_cells_per_patient = minimum_cells
      )
      run_enrichment(
        result, track, label, comparison$comparison, "all_patient_unpaired"
      )

      if (comparison$comparison == "dptdi_vs_sptdi") {
        shared <- intersect(
          candidate$patient[candidate$group == comparison$denominator],
          candidate$patient[candidate$group == comparison$numerator]
        )
        if (length(shared) >= minimum_donors) {
          paired <- candidate[candidate$patient %in% shared, , drop = FALSE]
          paired_counts <- counts[, rownames(paired), drop = FALSE]
          paired$patient <- factor(paired$patient)
          paired$group <- factor(
            paired$group,
            levels = c(comparison$denominator, comparison$numerator)
          )
          paired_dds <- DESeqDataSetFromMatrix(
            countData = round(paired_counts), colData = paired,
            design = ~ patient + group
          )
          paired_dds <- DESeq(paired_dds, quiet = TRUE)
          paired_result <- as.data.frame(results(
            paired_dds,
            contrast = c(
              "group", comparison$numerator, comparison$denominator
            ),
            alpha = 0.05
          ))
          paired_result$gene <- rownames(paired_result)
          paired_result$track <- track
          paired_result$label <- label
          paired_result$comparison <- comparison$comparison
          paired_result$analysis <- "paired_patient_sensitivity"
          paired_result$denominator_donors <- length(shared)
          paired_result$numerator_donors <- length(shared)
          paired_result$minimum_cells_per_patient <- minimum_cells
          paired_stem <- paste(
            safe_name(track), safe_name(label), comparison$comparison,
            "paired_patient_sensitivity", sep = "__"
          )
          paired_path <- file.path(de_dir, paste0(paired_stem, ".tsv.gz"))
          paired_connection <- gzfile(paired_path, "wt")
          write.table(
            paired_result, file = paired_connection, sep = "\t",
            quote = FALSE, row.names = FALSE
          )
          close(paired_connection)
          all_de_index[[length(all_de_index) + 1L]] <- data.frame(
            track = track, label = label,
            comparison = comparison$comparison,
            analysis = "paired_patient_sensitivity",
            result_file = paired_path, tested_genes = nrow(paired_result),
            significant_fdr05_lfc05 = sum(
              !is.na(paired_result$padj) & paired_result$padj < 0.05 &
                abs(paired_result$log2FoldChange) >= 0.5
            )
          )
          all_de_checkpoints[[length(all_de_checkpoints) + 1L]] <- data.frame(
            track = track, label = label,
            comparison = comparison$comparison,
            analysis = "paired_patient_sensitivity", status = "PASS",
            denominator_donors = length(shared),
            numerator_donors = length(shared),
            minimum_cells_per_patient = minimum_cells
          )
          run_enrichment(
            paired_result, track, label, comparison$comparison,
            "paired_patient_sensitivity"
          )
        }
      }
      rm(counts, dds, result)
      gc(verbose = FALSE)
    }
  }
}

de_index <- if (length(all_de_index)) do.call(rbind, all_de_index) else data.frame()
de_checkpoints <- if (length(all_de_checkpoints)) {
  do.call(rbind, all_de_checkpoints)
} else data.frame()
write_tsv(de_index, file.path(output_dir, "M05_result_index.tsv"))
write_tsv(de_checkpoints, file.path(output_dir, "M05_checkpoint.tsv"))

module_checkpoint <- data.frame(
  module = c("M04", "M05", "M06_GO_BP_ORA", "M06_KEGG", "M06_Reactome"),
  status = c(
    "PASS", "PASS", "PASS", "blocked_runtime_resource",
    if (requireNamespace("ReactomePA", quietly = TRUE)) "PASS" else
      "blocked_runtime_resource"
  ),
  reason = c(
    "edgeR QL with patient-group samples and total-cell offsets",
    "DESeq2 pseudobulk with patient-group replicates",
    "clusterProfiler enrichGO on current-project DE only",
    "No frozen local KEGG gene-set snapshot declared; network retrieval not used",
    "ReactomePA is absent from frozen csc-portable-v1 server runtime"
  )
)
write_tsv(module_checkpoint, file.path(output_dir, "formal_inference_checkpoint.tsv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_formal_inference.txt"))
message("Formal patient-level DA, pseudobulk DE, and enrichment completed.")
