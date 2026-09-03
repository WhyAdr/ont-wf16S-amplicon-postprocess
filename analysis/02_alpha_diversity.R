# =============================================================================
# Module 02: Alpha Diversity & Rarefaction Analysis
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(vegan)
})

run_alpha <- function(context) {
  cfg <- context$config
  alpha_dir <- cfg$output$dirs$alpha
  dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_outputs <- character(0)
  
  # Parameters
  seed <- cfg$seed %||% 42L
  n_points <- cfg$alpha$rarefaction_points %||% 25L
  resample_depth_cfg <- cfg$alpha$resample_depth %||% 50000L
  fraction_cap <- cfg$alpha$resample_fraction_cap %||% 0.90
  n_iterations <- cfg$alpha$resample_iterations %||% 100L
  
  unclass_idx <- context$unclass_index
  count_matrix <- context$count_matrix
  # Extract classified counts only
  class_matrix <- count_matrix[-unclass_idx, , drop = FALSE]
  
  samples <- context$samples
  
  # 1. Compute Alpha Diversity Indices per sample
  alpha_records <- list()
  for (s in samples) {
    counts_s <- class_matrix[, s]
    idx_df <- calc_alpha_indices(counts_s)
    idx_df$SampleID <- s
    alpha_records[[s]] <- idx_df
  }
  alpha_combined <- do.call(rbind, alpha_records)
  
  # Pivot to wide table for export: SampleID, Observed_Richness, Chao1, etc.
  alpha_wide <- alpha_combined %>%
    tidyr::pivot_wider(names_from = Metric, values_from = Value)
  
  # Attach metadata if available
  if (!is.null(context$metadata)) {
    alpha_wide <- merge(context$metadata, alpha_wide, by = "SampleID", sort = FALSE)
  }
  
  alpha_tsv <- file.path(alpha_dir, "alpha_diversity.tsv")
  write.table(alpha_wide, alpha_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, alpha_tsv)
  
  # 2. Analytical Rarefaction Curves
  rare_curve_records <- list()
  for (s in samples) {
    counts_s <- class_matrix[, s]
    curve_df <- calc_analytical_rarefaction(counts_s, n_points = n_points)
    if (nrow(curve_df) > 0) {
      curve_df$SampleID <- s
      rare_curve_records[[s]] <- curve_df
    }
  }
  rare_curve_combined <- do.call(rbind, rare_curve_records)
  
  rare_tsv <- file.path(alpha_dir, "rarefaction_curve.tsv")
  write.table(rare_curve_combined, rare_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, rare_tsv)
  
  # Figure 2a: Rarefaction Curve Plot
  p_curve <- ggplot(rare_curve_combined, aes(x = depth, y = mean_richness, group = SampleID, color = SampleID)) +
    geom_ribbon(aes(ymin = mean_richness - sd_richness, ymax = mean_richness + sd_richness, fill = SampleID),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0.01, 0.08))) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Analytical Rarefaction Curves",
      subtitle = "Expected species richness vs. sequencing depth (vegan::rarefy \u00b1 SE)",
      x = "Reads Subsampled", y = "Expected Species Richness",
      color = "Sample", fill = "Sample"
    ) +
    theme_amplicon()
  
  # For single sample, add actual classified depth dashed line
  if (length(samples) == 1) {
    single_depth <- sum(class_matrix[, samples[1]])
    p_curve <- p_curve +
      geom_vline(xintercept = single_depth, linetype = "dashed", color = "grey50") +
      annotate("text", x = single_depth, y = min(rare_curve_combined$mean_richness),
               label = "actual depth ", hjust = 1, vjust = 0, color = "grey40", size = 3) +
      theme(legend.position = "none")
  }
  
  p_curve_path <- file.path(alpha_dir, "02a_rarefaction_curve.png")
  save_plot(p_curve_path, p_curve, width = 8, height = 5.5)
  all_outputs <- c(all_outputs, p_curve_path)
  
  # 3. Rarefaction Resamples
  # Determine realized depth
  sample_depths <- colSums(class_matrix)
  min_depth <- min(sample_depths)
  
  realized_depth <- if (context$mode == "single") {
    min(resample_depth_cfg, floor(min_depth * fraction_cap))
  } else {
    floor(min_depth * fraction_cap)
  }
  
  resample_records <- list()
  for (s in samples) {
    counts_s <- class_matrix[, s]
    if (sum(counts_s) >= realized_depth && realized_depth > 0) {
      res_df <- calc_rarefaction_resamples(
        counts = counts_s,
        subsample_depth = realized_depth,
        n_iterations = n_iterations,
        seed = seed
      )
      res_df$SampleID <- s
      resample_records[[s]] <- res_df
    }
  }
  
  if (length(resample_records) > 0) {
    resamples_combined <- do.call(rbind, resample_records)
    resample_tsv <- file.path(alpha_dir, "rarefaction_resamples.tsv")
    write.table(resamples_combined, resample_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, resample_tsv)
    
    # Figure 2b: Resample Boxplots (Sensitivity of alpha metrics across resamples)
    resamples_long <- resamples_combined %>%
      select(SampleID, richness, shannon, ens, simpson, invsimpson) %>%
      tidyr::pivot_longer(cols = c(richness, shannon, ens, simpson, invsimpson),
                          names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = factor(
        Metric,
        levels = c("richness", "shannon", "ens", "simpson", "invsimpson"),
        labels = c("Richness (S)", "Shannon (H)", "ENS (e^H)", "Simpson (D)", "Inv. Simpson")
      ))
    
    p_box <- ggplot(resamples_long, aes(x = SampleID, y = Value, fill = SampleID)) +
      geom_boxplot(alpha = 0.7, outlier.size = 1) +
      facet_wrap(~Metric, scales = "free_y", nrow = 2) +
      labs(
        title = sprintf("Alpha Diversity Sensitivity Across %d Rarefaction Resamples", n_iterations),
        subtitle = sprintf("Standardized subsampling depth: %s reads (seed = %d)",
                            scales::comma(realized_depth), seed),
        x = NULL, y = "Value"
      ) +
      theme_amplicon() +
      theme(
        axis.text.x = if (length(samples) > 5) element_text(angle = 45, hjust = 1) else element_text(),
        legend.position = if (length(samples) == 1) "none" else "right"
      )
    
    p_box_path <- file.path(alpha_dir, "02b_resample_boxplots.png")
    save_plot(p_box_path, p_box, width = 9, height = 6)
    all_outputs <- c(all_outputs, p_box_path)
  }
  
  # 4. Cohort Group Tests (Only when mode == cohort)
  if (context$mode == "cohort") {
    meta <- context$metadata
    alpha_with_group <- merge(meta, alpha_wide, by = "SampleID")
    groups <- unique(alpha_with_group$Group)
    
    # Gating: At least 2 groups with at least 3 biological samples each
    group_counts <- table(alpha_with_group$Group)
    valid_groups <- names(group_counts)[group_counts >= 3]
    
    diff_file <- file.path(alpha_dir, "group_differences.tsv")
    
    if (length(valid_groups) < 2) {
      skip_note <- data.frame(
        Status = "Skipped",
        Reason = sprintf(
          "Alpha group comparison requires at least 2 groups with >= 3 samples. Found groups: %s",
          paste(sprintf("%s (n=%d)", names(group_counts), as.integer(group_counts)), collapse = ", ")
        ),
        stringsAsFactors = FALSE
      )
      write.table(skip_note, diff_file, sep = "\t", row.names = FALSE, quote = FALSE)
    } else {
      # Perform non-parametric tests on raw metrics
      metrics_to_test <- c("Observed species richness (S)", "Chao1 (estimated richness)",
                           "Shannon (H)", "Effective number of species (e^H)",
                           "Simpson's D (1-sum p^2)", "Inverse Simpson", "Pielou's evenness (J)")
      
      test_rows <- list()
      for (m in metrics_to_test) {
        if (m %in% colnames(alpha_with_group)) {
          df_sub <- alpha_with_group[alpha_with_group$Group %in% valid_groups, ]
          y <- df_sub[[m]]
          grp <- factor(df_sub$Group)
          
          if (length(valid_groups) == 2) {
            wt <- suppressWarnings(wilcox.test(y ~ grp))
            test_rows[[m]] <- data.frame(
              Metric = m, Test = "Wilcoxon rank-sum",
              Statistic = wt$statistic, P_Value = wt$p.value,
              stringsAsFactors = FALSE
            )
          } else {
            kt <- suppressWarnings(kruskal.test(y ~ grp))
            test_rows[[m]] <- data.frame(
              Metric = m, Test = "Kruskal-Wallis",
              Statistic = kt$statistic, P_Value = kt$p.value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      diff_df <- do.call(rbind, test_rows)
      diff_df$FDR_BH <- p.adjust(diff_df$P_Value, method = "BH")
      write.table(diff_df, diff_file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
    all_outputs <- c(all_outputs, diff_file)
  }
  
  list(
    status = "completed",
    outputs = all_outputs
  )
}
