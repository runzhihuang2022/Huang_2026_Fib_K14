.libPaths("data/external/local_835d389f6169")
suppressPackageStartupMessages(library(CellChat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript build_required_cellchat_chord_adapters.R <input_rds_dir> <output_rds_dir> <table_dir>")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
table_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

groups <- c("D0_7", "D15_28", "D8_14", "M1_2")
interfaces <- list(
  Endo = list(fib = c("Fib_Prolif", "Fib_Myo"), partner = c("Endo_Capillary", "Endo_Venous", "Endo_Arterial")),
  Immune = list(fib = c("Fib_Inflama", "Fib_CD45"), partner = c("M1_Macrophage", "M2_Macrophage", "Mast")),
  KC = list(fib = c("Fib_Papi", "Fib_K14"), partner = c("KC_Basal_Mig", "KC_Basal", "KC_Spinous_Mig")),
  SAC = list(fib = c("Fib_Papi", "Fib_K14"), partner = c("SAC_SG_Progenitor", "SAC_SG_Ductal", "SAC_SG_Dark"))
)

rank_lr <- function(lr) {
  if (!nrow(lr)) return(lr)
  lr <- lr[order(-lr$prob, lr$pval, lr$source, lr$target, lr$ligand, lr$receptor), , drop = FALSE]
  lr$key_id <- paste(lr$source, lr$ligand, lr$target, lr$receptor, sep = "|")
  lr[!duplicated(lr$key_id), , drop = FALSE]
}

select_with_required_coverage <- function(lr, required, n_display = 8L) {
  lr <- rank_lr(lr)
  if (!nrow(lr)) return(lr)
  picked <- integer(0)
  for (node in required) {
    idx <- which((lr$source == node | lr$target == node) & !(seq_len(nrow(lr)) %in% picked))
    if (length(idx)) picked <- c(picked, idx[[1]])
  }
  fill <- setdiff(seq_len(nrow(lr)), picked)
  picked <- c(picked, fill)[seq_len(min(n_display, nrow(lr)))]
  out <- lr[picked, , drop = FALSE]
  out$display_rank <- seq_len(nrow(out))
  out$selection_basis <- ifelse(
    seq_len(nrow(out)) <= length(unique(picked[seq_len(min(length(required), length(picked)))])),
    "required-cell coverage then locked probability rank",
    "locked probability rank"
  )
  out
}

records <- list()
all_selected <- list()
for (group in groups) {
  src <- file.path(input_dir, paste0("CellChat_", group, "_sub_labels.rds"))
  object <- readRDS(src)
  original_prob <- object@net$prob
  original_pval <- object@net$pval
  masked_prob <- original_prob * 0
  masked_pval <- original_pval * 0 + 1
  counts <- table(object@idents)

  for (interface_name in names(interfaces)) {
    spec <- interfaces[[interface_name]]
    required <- c(spec$fib, spec$partner)
    present <- required %in% names(counts)
    n_required <- setNames(rep(0L, length(required)), required)
    n_required[present] <- as.integer(counts[required[present]])
    eligible <- all(n_required >= 10L)
    lr <- data.frame()
    selected <- data.frame()
    reason <- ""
    if (eligible) {
      forward <- tryCatch(subsetCommunication(object, sources.use = spec$fib, targets.use = spec$partner), error = function(e) data.frame())
      reverse <- tryCatch(subsetCommunication(object, sources.use = spec$partner, targets.use = spec$fib), error = function(e) data.frame())
      lr <- rank_lr(rbind(forward, reverse))
      selected <- select_with_required_coverage(lr, required, 8L)
      seen <- unique(c(as.character(selected$source), as.character(selected$target)))
      missing <- setdiff(required, seen)
      if (nrow(selected) < 8L) reason <- paste0("Only ", nrow(selected), " significant LR pairs; minimum requested is 8")
      if (length(missing)) reason <- paste0(ifelse(nchar(reason), paste0(reason, "; "), ""), "No significant displayed pair for: ", paste(missing, collapse = ","))
      if (nrow(selected) >= 8L && !length(missing)) {
        for (i in seq_len(nrow(selected))) {
          s <- as.character(selected$source[[i]])
          t <- as.character(selected$target[[i]])
          k <- as.character(selected$interaction_name[[i]])
          if (s %in% dimnames(masked_prob)[[1]] && t %in% dimnames(masked_prob)[[2]] && k %in% dimnames(masked_prob)[[3]]) {
            masked_prob[s, t, k] <- original_prob[s, t, k]
            masked_pval[s, t, k] <- original_pval[s, t, k]
          }
        }
      }
    } else {
      low <- names(n_required)[n_required < 10L]
      reason <- paste0("Required subtype count <10: ", paste(paste0(low, " n=", n_required[low]), collapse = "; "))
    }

    if (nrow(lr)) write.table(lr, file.path(table_dir, paste0(group, "_Fib_", interface_name, "_all_significant_LR.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    if (nrow(selected)) {
      selected$group <- group
      selected$interface <- paste0("Fib_", interface_name)
      selected$required_cells <- paste(required, collapse = ";")
      write.table(selected, file.path(table_dir, paste0(group, "_Fib_", interface_name, "_displayed_LR.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
      all_selected[[length(all_selected) + 1L]] <- selected
    }
    status <- if (eligible && nrow(selected) >= 8L && !nchar(reason)) "PASS" else "NOT_ESTIMABLE"
    records[[length(records) + 1L]] <- data.frame(
      group = group,
      interface = paste0("Fib_", interface_name),
      required_cells = paste(required, collapse = ";"),
      required_counts = paste(paste0(required, "=", n_required[required]), collapse = ";"),
      status = status,
      n_all_significant_lr = nrow(lr),
      n_display_lr = nrow(selected),
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  object@net$prob <- masked_prob
  object@net$pval <- masked_pval
  saveRDS(object, file.path(output_dir, paste0("CellChat_", group, "_sub_labels.rds")), compress = TRUE)
}

audit <- do.call(rbind, records)
write.table(audit, file.path(table_dir, "REQUIRED_CHORD_INSTANCE_AUDIT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(all_selected)) {
  write.table(do.call(rbind, all_selected), file.path(table_dir, "ALL_DISPLAYED_REQUIRED_CHORD_LR.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
