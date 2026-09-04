# =============================================================================
# Module 06: Shared, Core, and Unique Taxa Analysis (Prevalence, UpSet)
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(UpSetR)
})

run_shared_taxa <- function(context) {
  cfg <- context$config
  shared_dir <- cfg$output$dirs$shared_taxa
  dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)

  all_outputs <- character(0)

  # Cohort gate: Requires at least 2 samples
  if (context$mode != "cohort" || length(context$samples) < 2) {
    skip_file <- file.path(shared_dir, "shared_taxa_skipped.tsv")
    skip_df <- data.frame(
      Status = "Skipped",
      Reason = sprintf("Shared taxa module requires >= 2 samples. Current mode is '%s' with %d sample(s).",
                       context$mode, length(context$samples)),
      stringsAsFactors = FALSE
    )
    write.table(skip_df, skip_file, sep = "\t", row.names = FALSE, quote = FALSE)
    return(list(
      status = "skipped",
      reason = skip_df$Reason[1],
      outputs = skip_file
    ))
  }

  samples <- context$samples
  unclass_idx <- context$unclass_index
  count_matrix <- context$count_matrix
  class_counts <- count_matrix[-unclass_idx, , drop = FALSE]
  tax_df <- context$taxonomy[-unclass_idx, , drop = FALSE]

  meta <- context$metadata
  min_count <- cfg$shared_taxa$minimum_count %||% 1L
  group_prev_thresh <- cfg$shared_taxa$group_prevalence %||% 0.5
  analysis_rank <- cfg$shared_taxa$rank %||% "species"

  # Group by analysis rank
  rk_idx <- which(colnames(context$taxonomy) == analysis_rank)
  if (length(rk_idx) != 1L) {
    stop(sprintf("Unsupported shared-taxa rank: '%s'", analysis_rank), call. = FALSE)
  }

  rk_prefixes <- apply(tax_df[, 1:rk_idx, drop = FALSE], 1, paste, collapse = ";")
  tax_names <- tax_df[[rk_idx]]

  agg_counts <- data.frame(TaxonPath = rk_prefixes, Taxon = tax_names, class_counts, check.names = FALSE) %>%
    group_by(TaxonPath, Taxon) %>%
    summarise(across(all_of(samples), sum), .groups = "drop")

  tax_labels <- agg_counts$Taxon
  mat_counts <- as.matrix(agg_counts[, samples, drop = FALSE])
  rownames(mat_counts) <- agg_counts$TaxonPath

  # 1. Sample Presence/Absence Matrix
  presence_mat <- (mat_counts >= min_count) * 1L
  presence_df <- data.frame(
    TaxonPath = agg_counts$TaxonPath,
    Taxon = tax_labels,
    presence_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  presence_file <- file.path(shared_dir, "sample_presence.tsv")
  write.table(presence_df, presence_file, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, presence_file)

  # 2. Group Prevalence Table & Membership Matrix
  groups <- unique(meta$Group)
  prev_list <- list()
  membership_list <- list()

  for (grp in groups) {
    grp_samples <- meta$SampleID[meta$Group == grp]
    n_grp <- length(grp_samples)

    if (n_grp > 0) {
      grp_pres <- presence_mat[, grp_samples, drop = FALSE]
      pos_count <- rowSums(grp_pres)
      prevalence <- pos_count / n_grp
      prev_list[[grp]] <- prevalence
      membership_list[[grp]] <- (prevalence >= group_prev_thresh) * 1L
    }
  }

  prev_mat <- do.call(cbind, prev_list)
  prev_df <- data.frame(
    TaxonPath = agg_counts$TaxonPath,
    Taxon = tax_labels,
    prev_mat,
    stringsAsFactors = FALSE
  )
  prev_file <- file.path(shared_dir, "group_prevalence.tsv")
  write.table(prev_df, prev_file, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, prev_file)

  mem_mat <- do.call(cbind, membership_list)
  mem_df <- data.frame(
    TaxonPath = agg_counts$TaxonPath,
    Taxon = tax_labels,
    Threshold = sprintf("Prevalence >= %.2f within group", group_prev_thresh),
    mem_mat,
    stringsAsFactors = FALSE
  )
  mem_file <- file.path(shared_dir, "group_membership.tsv")
  write.table(mem_df, mem_file, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, mem_file)

  # 3. Core & Unique Taxa
  n_groups_present <- rowSums(mem_mat)
  is_core <- (n_groups_present == length(groups))
  is_unique <- (n_groups_present == 1)

  core_df <- data.frame(
    TaxonPath = agg_counts$TaxonPath[is_core],
    Taxon = tax_labels[is_core],
    Threshold = sprintf("count >= %d per sample; present in all %d groups at prevalence >= %.2f",
                        min_count, length(groups), group_prev_thresh),
    stringsAsFactors = FALSE
  )
  core_file <- file.path(shared_dir, "core_taxa.tsv")
  write.table(core_df, core_file, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, core_file)

  unique_group_name <- apply(mem_mat[is_unique, , drop = FALSE], 1, function(row) {
    names(row)[which(row == 1)[1]]
  })
  unique_df <- data.frame(
    TaxonPath = agg_counts$TaxonPath[is_unique],
    Taxon = tax_labels[is_unique],
    ExclusiveGroup = unique_group_name,
    Threshold = sprintf("count >= %d per sample; present exclusively in 1 group at prevalence >= %.2f",
                        min_count, group_prev_thresh),
    stringsAsFactors = FALSE
  )
  unique_file <- file.path(shared_dir, "unique_taxa.tsv")
  write.table(unique_df, unique_file, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, unique_file)

  # 4. UpSet Plot (When >= 2 non-empty groups exist)
  if (length(groups) >= 2 && sum(colSums(mem_mat) > 0) >= 2) {
    upset_df <- as.data.frame(mem_mat)
    upset_path <- file.path(shared_dir, "06_upset_plot.png")

    png(upset_path, width = 8, height = 5.5, units = "in", res = 150)
    suppressWarnings(print(UpSetR::upset(
      upset_df,
      sets = groups,
      order.by = "freq",
      mainbar.y.label = "Shared Taxa Intersections",
      sets.x.label = sprintf("Taxa per Group (prev >= %.1f)", group_prev_thresh),
      text.scale = 1.2
    )))
    dev.off()
    all_outputs <- c(all_outputs, upset_path)
  }

  list(
    status = "completed",
    outputs = all_outputs
  )
}
