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

build_richness_overview <- function(class_matrix, samples) {
  do.call(rbind, lapply(samples, function(sample_id) {
    positive <- class_matrix[, sample_id][class_matrix[, sample_id] > 0]
    classified_reads <- sum(positive)
    singleton_taxa <- sum(positive == 1)
    low_count <- positive[positive <= 10]
    data.frame(
      SampleID = sample_id,
      ClassifiedReads = classified_reads,
      PositiveTaxa = length(positive),
      SingletonTaxa = singleton_taxa,
      SingletonPct = 100 * singleton_taxa / length(positive),
      TaxaLeq10 = length(low_count),
      ReadsInTaxaLeq10 = sum(low_count),
      ReadsInTaxaLeq10Pct = 100 * sum(low_count) / classified_reads,
      stringsAsFactors = FALSE
    )
  }))
}

write_alpha_tsv <- function(value, path) {
  write.table(value, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  path
}

alpha_plot_subtitle <- function(depth, iterations, seed) {
  sprintf(
    paste0("Standardized depth: %s reads; %d iterations; base seed %d with deterministic sample-derived seeds\n",
           "Rarefaction sensitivity only; iterations are not biological replicates"),
    scales::comma(depth), iterations, seed
  )
}

build_alpha_resample_plot <- function(long_table, registry, figure_key, title,
                                      samples, depth, iterations, seed) {
  contract <- registry[registry$FigureKey == figure_key, , drop = FALSE]
  contract <- contract[order(contract$FacetOrder), , drop = FALSE]
  valid_keys <- contract$MetricKey[vapply(contract$MetricKey, function(key) {
    any(long_table$MetricKey == key & long_table$Valid)
  }, logical(1))]
  plot_data <- long_table[long_table$MetricKey %in% valid_keys & long_table$Valid, , drop = FALSE]
  if (!nrow(plot_data)) return(NULL)
  labels <- stats::setNames(contract$MetricLabel, contract$MetricKey)
  plot_data$MetricFacet <- factor(plot_data$MetricKey, levels = valid_keys,
                                  labels = unname(labels[valid_keys]))
  ggplot(plot_data, aes(x = SampleID, y = Value, fill = SampleID)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    facet_wrap(~MetricFacet, scales = "free_y", nrow = 2) +
    labs(
      title = title,
      subtitle = alpha_plot_subtitle(depth, iterations, seed),
      x = NULL, y = "Value"
    ) +
    theme_amplicon() +
    theme(
      axis.text.x = if (length(samples) > 5) element_text(angle = 45, hjust = 1) else element_text(),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = if (length(samples) == 1) "none" else "right"
    )
}

run_alpha <- function(context) {
  cfg <- context$config
  alpha_dir <- cfg$output$dirs$alpha
  dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)
  all_outputs <- character(0)

  seed <- cfg$seed %||% 42L
  n_points <- cfg$alpha$rarefaction_points %||% 25L
  resample_depth_cfg <- cfg$alpha$resample_depth %||% 50000L
  fraction_cap <- cfg$alpha$resample_fraction_cap %||% 0.90
  n_iterations <- cfg$alpha$resample_iterations %||% 100L
  hill_orders <- cfg$alpha$hill_orders %||% c(0, 1, 2)
  renyi_orders <- cfg$alpha$renyi_orders %||% c(1)
  registry <- alpha_metric_registry(hill_orders, renyi_orders)

  count_matrix <- context$count_matrix
  class_matrix <- count_matrix[-context$unclass_index, , drop = FALSE]
  taxon_paths <- context$taxonomy$TaxonPath[-context$unclass_index]
  rownames(class_matrix) <- taxon_paths
  samples <- context$samples
  positive_paths <- taxon_paths[rowSums(class_matrix) > 0]
  phylogeny <- prepare_alpha_phylogeny(cfg, positive_paths)

  definitions_path <- file.path(alpha_dir, "alpha_metric_definitions.tsv")
  write_alpha_tsv(registry, definitions_path)
  all_outputs <- c(all_outputs, definitions_path)

  phylo_status_path <- file.path(alpha_dir, "alpha_phylogeny_status.tsv")
  write_alpha_tsv(alpha_phylogeny_provenance(phylogeny), phylo_status_path)
  all_outputs <- c(all_outputs, phylo_status_path)

  alpha_records <- list()
  alpha_status_records <- list()
  for (sample_id in samples) {
    counts <- class_matrix[, sample_id]
    non_phylo <- calc_alpha_metric_long(counts, hill_orders, renyi_orders)
    status <- if (isTRUE(phylogeny$enabled)) {
      rbind(non_phylo, calc_phylogenetic_alpha(counts, phylogeny, registry))
    } else non_phylo
    status$SampleID <- sample_id
    alpha_status_records[[sample_id]] <- status[, c(
      "SampleID", "MetricKey", "MetricLabel", "Parameter", "Value", "Valid", "Reason"
    )]

    conventional <- calc_alpha_indices(counts, hill_orders, renyi_orders)
    if (isTRUE(phylogeny$enabled)) {
      phylo_rows <- status[status$MetricKey %in% c("faith_pd", "psr", "pse"), ]
      conventional <- rbind(conventional, data.frame(
        Metric = phylo_rows$MetricLabel, Value = phylo_rows$Value,
        stringsAsFactors = FALSE
      ))
    }
    conventional$SampleID <- sample_id
    alpha_records[[sample_id]] <- conventional
  }
  alpha_status <- do.call(rbind, alpha_status_records)
  alpha_status_path <- file.path(alpha_dir, "alpha_metric_status.tsv")
  write_alpha_tsv(alpha_status, alpha_status_path)
  all_outputs <- c(all_outputs, alpha_status_path)

  alpha_combined <- do.call(rbind, alpha_records)
  alpha_wide <- alpha_combined[, c("SampleID", "Metric", "Value")] %>%
    tidyr::pivot_wider(names_from = Metric, values_from = Value)
  alpha_wide <- alpha_wide[match(samples, alpha_wide$SampleID), , drop = FALSE]
  if (!is.null(context$metadata)) {
    collisions <- intersect(setdiff(names(context$metadata), "SampleID"),
                            setdiff(names(alpha_wide), "SampleID"))
    if (length(collisions)) stop(sprintf("Metadata/output column collision before alpha join: %s",
                                         paste(collisions, collapse = ", ")), call. = FALSE)
    alpha_wide <- dplyr::left_join(context$metadata, alpha_wide, by = "SampleID")
  }
  alpha_tsv <- file.path(alpha_dir, "alpha_diversity.tsv")
  write_alpha_tsv(alpha_wide, alpha_tsv)
  all_outputs <- c(all_outputs, alpha_tsv)

  richness_overview <- build_richness_overview(class_matrix, samples)
  richness_tsv <- file.path(alpha_dir, "02_richness_overview.tsv")
  write_alpha_tsv(richness_overview, richness_tsv)
  all_outputs <- c(all_outputs, richness_tsv)

  rare_curve_records <- list()
  for (sample_id in samples) {
    curve <- calc_analytical_rarefaction(class_matrix[, sample_id], n_points = n_points)
    if (nrow(curve)) {
      curve$SampleID <- sample_id
      rare_curve_records[[sample_id]] <- curve
    }
  }
  rare_curve_combined <- do.call(rbind, rare_curve_records)
  rare_tsv <- file.path(alpha_dir, "rarefaction_curve.tsv")
  write_alpha_tsv(rare_curve_combined, rare_tsv)
  all_outputs <- c(all_outputs, rare_tsv)

  p_curve <- ggplot(rare_curve_combined,
                    aes(x = depth, y = mean_richness, group = SampleID, color = SampleID)) +
    geom_ribbon(aes(ymin = mean_richness - sd_richness,
                    ymax = mean_richness + sd_richness, fill = SampleID),
                alpha = 0.2, color = NA) +
    geom_point(size = 1.5) +
    scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0.01, 0.08))) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Analytical Rarefaction Curves",
         subtitle = "Expected species richness vs. sequencing depth (vegan::rarefy +/- SE)",
         x = "Reads Subsampled", y = "Expected Species Richness", color = "Sample", fill = "Sample") +
    theme_amplicon()
  if (any(table(rare_curve_combined$SampleID) > 1L)) {
    p_curve <- p_curve + geom_line(linewidth = 1)
  }
  if (length(samples) == 1L) {
    single_depth <- sum(class_matrix[, samples[[1]]])
    p_curve <- p_curve +
      geom_vline(xintercept = single_depth, linetype = "dashed", color = "grey50") +
      annotate("text", x = single_depth, y = min(rare_curve_combined$mean_richness),
               label = "actual depth ", hjust = 1, vjust = 0, color = "grey40", size = 3) +
      theme(legend.position = "none")
  }
  p_curve_path <- file.path(alpha_dir, "02a_rarefaction_curve.png")
  save_plot(p_curve_path, p_curve, width = 8, height = 5.5, dpi = 150)
  all_outputs <- c(all_outputs, p_curve_path)

  sample_depths <- colSums(class_matrix)
  # A one-read classified sample still has a valid depth-one sensitivity
  # distribution. Avoid rounding the fractional cap down to an unusable zero.
  realized_depth <- min(resample_depth_cfg,
                        max(1L, floor(min(sample_depths) * fraction_cap)))
  resample_records <- list()
  resample_long_records <- list()
  for (sample_id in samples) {
    counts <- class_matrix[, sample_id]
    if (sum(counts) >= realized_depth && realized_depth > 0) {
      resampled <- calc_rarefaction_resamples(
        counts, realized_depth, n_iterations, derive_sample_seed(seed, sample_id),
        hill_orders, renyi_orders, phylogeny
      )
      long <- attr(resampled, "alpha_long")
      attr(resampled, "alpha_long") <- NULL
      resampled$SampleID <- sample_id
      long$SampleID <- sample_id
      resample_records[[sample_id]] <- resampled
      resample_long_records[[sample_id]] <- long[, c(
        "SampleID", "iteration", "subsample_depth", "MetricKey", "MetricLabel",
        "Parameter", "Value", "Valid", "Reason"
      )]
    }
  }

  if (length(resample_records)) {
    resamples_combined <- do.call(rbind, resample_records)
    resamples_long <- do.call(rbind, resample_long_records)
    resample_tsv <- file.path(alpha_dir, "rarefaction_resamples.tsv")
    resample_long_tsv <- file.path(alpha_dir, "rarefaction_resamples_long.tsv")
    write_alpha_tsv(resamples_combined, resample_tsv)
    write_alpha_tsv(resamples_long, resample_long_tsv)
    all_outputs <- c(all_outputs, resample_tsv, resample_long_tsv)

    figure_specs <- list(
      `02b` = list(file = "02b_resample_boxplots.png", width = 10, height = 7,
                   title = sprintf("Alpha Diversity Sensitivity Across %d Rarefaction Resamples", n_iterations)),
      `02c` = list(file = "02c_resample_boxplots_estimators_evenness.png", width = 10, height = 7,
                   title = "Richness Estimator and Evenness Sensitivity"),
      `02d` = list(file = "02d_resample_boxplots_evenness_entropy.png", width = 10, height = 7,
                   title = "Evenness, Diversity, and Entropy Sensitivity")
    )
    if (isTRUE(phylogeny$enabled)) {
      figure_specs$`02e` <- list(file = "02e_resample_boxplots_phylogenetic.png",
                                width = 10, height = 7,
                                title = "Hill and Phylogenetic Diversity Sensitivity")
    }
    extra_renyi <- registry$MetricKey[registry$FigureKey == "02f"]
    if (length(extra_renyi)) {
      figure_specs$`02f` <- list(file = "02f_resample_renyi_profile.png",
                                width = 10, height = 7,
                                title = "R\u00e9nyi Entropy Profile Sensitivity")
    }

    plot_registry <- list()
    for (figure_key in names(figure_specs)) {
      spec <- figure_specs[[figure_key]]
      plot <- build_alpha_resample_plot(
        resamples_long, registry, figure_key, spec$title, samples,
        realized_depth, n_iterations, seed
      )
      if (is.null(plot)) next
      path <- file.path(alpha_dir, spec$file)
      save_plot(path, plot, width = spec$width, height = spec$height, dpi = 150)
      all_outputs <- c(all_outputs, path)
      contract <- registry[registry$FigureKey == figure_key, , drop = FALSE]
      contract <- contract[order(contract$FacetOrder), , drop = FALSE]
      contract$FigureFile <- spec$file
      contract$WidthIn <- spec$width
      contract$HeightIn <- spec$height
      contract$DPI <- 150L
      contract$Status <- vapply(contract$MetricKey, function(key) {
        if (any(resamples_long$MetricKey == key & resamples_long$Valid)) "Completed" else "Ineligible"
      }, character(1))
      plot_registry[[figure_key]] <- contract[, c(
        "FigureKey", "FigureFile", "FacetOrder", "MetricKey", "MetricLabel",
        "Parameter", "WidthIn", "HeightIn", "DPI", "Status"
      )]
    }
    if (!isTRUE(phylogeny$enabled)) {
      contract <- registry[registry$FigureKey == "02e", , drop = FALSE]
      contract$FigureFile <- "02e_resample_boxplots_phylogenetic.png"
      contract$WidthIn <- 10
      contract$HeightIn <- 7
      contract$DPI <- 150L
      contract$Status <- "Skipped"
      plot_registry$`02e` <- contract[, c(
        "FigureKey", "FigureFile", "FacetOrder", "MetricKey", "MetricLabel",
        "Parameter", "WidthIn", "HeightIn", "DPI", "Status"
      )]
      legacy_skip_path <- file.path(alpha_dir, "alpha_phylogeny_skipped.tsv")
      write_alpha_tsv(phylogeny$skip, legacy_skip_path)
      all_outputs <- c(all_outputs, legacy_skip_path)
    }
    plot_registry_path <- file.path(alpha_dir, "alpha_plot_registry.tsv")
    write_alpha_tsv(do.call(rbind, plot_registry), plot_registry_path)
    all_outputs <- c(all_outputs, plot_registry_path)
  }

  if (context$mode == "cohort") {
    alpha_with_group <- alpha_wide
    group_counts <- table(alpha_with_group$Group)
    all_groups_replicated <- length(group_counts) >= 2L && all(group_counts >= 3L)
    diff_file <- file.path(alpha_dir, "group_differences.tsv")
    if (!all_groups_replicated) {
      skip_note <- data.frame(
        Status = "Skipped",
        Reason = sprintf(
          "Alpha group comparison requires at least 2 groups with >= 3 samples. Found groups: %s",
          paste(sprintf("%s (n=%d)", names(group_counts), as.integer(group_counts)), collapse = ", ")
        ), stringsAsFactors = FALSE
      )
      write_alpha_tsv(skip_note, diff_file)
    } else {
      metrics_to_test <- c(
        "Observed species richness (S)", "Chao1 (estimated richness)", "Shannon (H)",
        "Effective number of species (e^H)", "Simpson's D (1-sum p^2)",
        "Inverse Simpson", "Pielou's evenness (J)"
      )
      test_rows <- list()
      for (metric in metrics_to_test) {
        if (!metric %in% colnames(alpha_with_group)) next
        subset <- alpha_with_group[is.finite(alpha_with_group[[metric]]), , drop = FALSE]
        y <- subset[[metric]]
        group <- factor(subset$Group)
        counts <- table(group)
        if (length(counts) < 2L || any(counts < 3L)) next
        if (nlevels(group) == 2L) {
          test <- wilcox.test(y ~ group)
          test_name <- "Wilcoxon rank-sum"
        } else {
          test <- kruskal.test(y ~ group)
          test_name <- "Kruskal-Wallis"
        }
        test_rows[[metric]] <- data.frame(
          Metric = metric, Test = test_name, Statistic = unname(test$statistic),
          P_Value = test$p.value, stringsAsFactors = FALSE
        )
      }
      if (!length(test_rows)) {
        write_alpha_tsv(data.frame(
          Status = "Skipped",
          Reason = "No metric retained at least 3 finite observations in every group.",
          stringsAsFactors = FALSE
        ), diff_file)
      } else {
        differences <- do.call(rbind, test_rows)
        differences$FDR_BH <- p.adjust(differences$P_Value, method = "BH")
        write_alpha_tsv(differences, diff_file)
      }
    }
    all_outputs <- c(all_outputs, diff_file)
  }

  list(status = "completed", outputs = all_outputs)
}
