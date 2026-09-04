# =============================================================================
# Module 03: Beta Diversity (Bray-Curtis, Jaccard, PCoA, PERMANOVA, Betadisper)
# =============================================================================

suppressMessages({
  library(vegan)
  library(dplyr)
  library(ggplot2)
})

run_beta <- function(context) {
  cfg <- context$config
  beta_dir <- cfg$output$dirs$beta
  dir.create(beta_dir, recursive = TRUE, showWarnings = FALSE)

  all_outputs <- character(0)

  # Cohort gate
  if (context$mode != "cohort" || length(context$samples) < 2) {
    skip_file <- file.path(beta_dir, "beta_diversity_skipped.tsv")
    skip_df <- data.frame(
      Status = "Skipped",
      Reason = sprintf("Beta diversity is cohort-only. Current mode is '%s' with %d sample(s).",
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

  meta <- context$metadata
  seed <- cfg$seed %||% 42L
  set.seed(seed)

  # Samples as rows, taxa as columns
  otu_table <- t(class_counts)
  # Calculate classified relative abundances for Bray-Curtis
  sample_sums <- rowSums(otu_table)
  rel_otu <- sweep(otu_table, 1, sample_sums, "/")

  distances_cfg <- cfg$beta$distances %||% c("bray", "jaccard")
  n_perm <- cfg$beta$permutations %||% 999L
  minimum_count <- cfg$beta$minimum_count %||% 1L
  strata_col <- cfg$beta$strata_column

  dist_list <- list()

  for (d_name in distances_cfg) {
    dist_mat <- if (d_name == "bray") {
      vegan::vegdist(rel_otu, method = "bray")
    } else if (d_name == "jaccard") {
      vegan::vegdist(otu_table >= minimum_count, method = "jaccard", binary = TRUE)
    } else {
      stop(sprintf("Unsupported beta-diversity distance: '%s'", d_name), call. = FALSE)
    }
    dist_list[[d_name]] <- dist_mat

    # Save distance matrix TSV
    d_tsv <- file.path(beta_dir, sprintf("distance_%s.tsv", d_name))
    distance_df <- data.frame(
      SampleID = rownames(as.matrix(dist_mat)),
      as.matrix(dist_mat),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    write.table(distance_df, d_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, d_tsv)

    distance_values <- as.vector(dist_mat)
    if (!any(is.finite(distance_values) & distance_values > 0)) {
      skip_path <- file.path(beta_dir, sprintf("pcoa_skipped_%s.tsv", d_name))
      write.table(data.frame(
        Status = "Skipped",
        Reason = "All pairwise distances are zero; PCoA is undefined."
      ), skip_path, sep = "\t", row.names = FALSE, quote = FALSE)
      all_outputs <- c(all_outputs, skip_path)
      next
    }

    # PCoA via additive-corrected classical scaling.
    max_k <- max(1, min(nrow(otu_table) - 1, 2))
    pcoa <- tryCatch(
      stats::cmdscale(dist_mat, k = max_k, eig = TRUE, add = TRUE),
      error = function(e) structure(list(message = conditionMessage(e)), class = "pcoa_error")
    )
    if (inherits(pcoa, "pcoa_error")) {
      skip_path <- file.path(beta_dir, sprintf("pcoa_skipped_%s.tsv", d_name))
      write.table(data.frame(Status = "Skipped", Reason = pcoa$message), skip_path,
                  sep = "\t", row.names = FALSE, quote = FALSE)
      all_outputs <- c(all_outputs, skip_path)
      next
    }
    eig <- pcoa$eig
    pos_eig <- eig[eig > 0]
    total_pos <- if (length(pos_eig) > 0) sum(pos_eig) else 1

    var_exp <- c(
      if (length(eig) >= 1 && eig[1] > 0) round(100 * eig[1] / total_pos, 1) else 0,
      if (length(eig) >= 2 && eig[2] > 0) round(100 * eig[2] / total_pos, 1) else 0
    )

    pts <- as.matrix(pcoa$points)
    p1 <- if (ncol(pts) >= 1) pts[, 1] else rep(0, nrow(otu_table))
    p2 <- if (ncol(pts) >= 2) pts[, 2] else rep(0, nrow(otu_table))

    scores_df <- data.frame(
      SampleID = rownames(otu_table),
      PCoA1 = p1,
      PCoA2 = p2,
      stringsAsFactors = FALSE
    )
    if (!is.null(meta)) {
      scores_df <- dplyr::left_join(meta, scores_df, by = "SampleID")
    }

    scores_tsv <- file.path(beta_dir, sprintf("pcoa_scores_%s.tsv", d_name))
    write.table(scores_df, scores_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, scores_tsv)

    variance_tsv <- file.path(beta_dir, sprintf("pcoa_variance_%s.tsv", d_name))
    variance_df <- data.frame(
      Axis = paste0("PCoA", seq_along(eig)),
      Eigenvalue = eig,
      PositiveVariancePercent = ifelse(eig > 0, 100 * eig / total_pos, 0),
      stringsAsFactors = FALSE
    )
    write.table(variance_df, variance_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, variance_tsv)

    # PCoA 2D plot (only if at least 3 samples exist)
    n_unique_samples <- nrow(unique(as.data.frame(rel_otu)))
    if (length(samples) >= 3 && n_unique_samples >= 3 && ncol(pts) >= 2) {
      p_pcoa <- ggplot(scores_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
        geom_point(size = 3.5, alpha = 0.85) +
        labs(
          title = sprintf("PCoA of %s Dissimilarity", toupper(d_name)),
          subtitle = sprintf("Full-data deterministic ordination (seed = %d)", seed),
          x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]),
          y = sprintf("PCoA 2 (%.1f%%)", var_exp[2])
        ) +
        theme_amplicon()

      p_path <- file.path(beta_dir, sprintf("03_pcoa_%s.png", d_name))
      save_plot(p_path, p_pcoa, width = 7, height = 5.5)
      all_outputs <- c(all_outputs, p_path)
    }
  }

  # Optional rarefaction stability analysis. Full-data distances/PCoA above
  # remain primary; these coordinates are Procrustes-aligned sensitivity draws.
  if (isTRUE(cfg$beta$resampling$enabled)) {
    stability_file <- file.path(beta_dir, "pcoa_rarefaction_stability.tsv")
    stability_diag <- file.path(beta_dir, "pcoa_rarefaction_diagnostics.tsv")
    stability_depth <- floor(min(sample_sums) * cfg$beta$resampling$depth_fraction_of_minimum)
    reference_dist <- vegan::vegdist(rel_otu, method = "bray")
    reference_fit <- tryCatch(
      stats::cmdscale(reference_dist, k = 2L, eig = TRUE, add = TRUE),
      error = function(e) NULL
    )

    if (nrow(otu_table) < 3L || stability_depth < 1L || is.null(reference_fit) ||
        ncol(as.matrix(reference_fit$points)) < 2L) {
      write.table(data.frame(
        Status = "Skipped",
        Reason = "Rarefaction stability requires >=3 samples and a two-axis full-data PCoA solution."
      ), stability_diag, sep = "\t", row.names = FALSE, quote = FALSE)
      all_outputs <- c(all_outputs, stability_diag)
    } else {
      set.seed(seed)
      stability_rows <- list()
      failed_iterations <- integer(0)
      reference_points <- as.matrix(reference_fit$points)[, 1:2, drop = FALSE]
      for (iteration in seq_len(cfg$beta$resampling$iterations)) {
        rare_counts <- vegan::rrarefy(otu_table, sample = stability_depth)
        rare_rel <- sweep(rare_counts, 1, rowSums(rare_counts), "/")
        rare_dist <- vegan::vegdist(rare_rel, method = "bray")
        rare_fit <- tryCatch(
          stats::cmdscale(rare_dist, k = 2L, eig = TRUE, add = TRUE),
          error = function(e) NULL
        )
        if (is.null(rare_fit) || ncol(as.matrix(rare_fit$points)) < 2L) {
          failed_iterations <- c(failed_iterations, iteration)
          next
        }
        aligned <- vegan::procrustes(reference_points,
                                     as.matrix(rare_fit$points)[, 1:2, drop = FALSE])$Yrot
        stability_rows[[length(stability_rows) + 1L]] <- data.frame(
          Iteration = iteration,
          SampleID = rownames(aligned),
          PCoA1 = aligned[, 1],
          PCoA2 = aligned[, 2],
          Depth = stability_depth,
          stringsAsFactors = FALSE
        )
      }
      if (length(stability_rows) == 0L) {
        stop("All beta-diversity rarefaction stability iterations failed.", call. = FALSE)
      }
      write.table(do.call(rbind, stability_rows), stability_file,
                  sep = "\t", row.names = FALSE, quote = FALSE)
      write.table(data.frame(
        Status = "Completed",
        RequestedIterations = cfg$beta$resampling$iterations,
        SuccessfulIterations = length(stability_rows),
        FailedIterations = paste(failed_iterations, collapse = ","),
        Depth = stability_depth,
        Seed = seed
      ), stability_diag, sep = "\t", row.names = FALSE, quote = FALSE)
      all_outputs <- c(all_outputs, stability_file, stability_diag)
    }
  }

  # PERMANOVA & Betadisper (Primary distance: Bray-Curtis)
  primary_dist <- dist_list[["bray"]] %||% dist_list[[1]]

  permanova_file <- file.path(beta_dir, "permanova.tsv")
  betadisper_file <- file.path(beta_dir, "betadisper.tsv")

  # Gating: At least 2 groups with at least 2 samples per group
  group_counts <- table(meta$Group)
  primary_values <- as.vector(primary_dist)
  residual_df <- nrow(meta) - length(group_counts)
  can_run_permanova <- length(group_counts) >= 2 && all(group_counts >= 2) &&
    residual_df > 0 && any(is.finite(primary_values) & primary_values > 0)

  if (can_run_permanova) {
    # Check strata
    strata_vec <- if (!is.null(strata_col)) meta[[strata_col]] else NULL

    set.seed(seed)
    perm_res <- suppressWarnings(vegan::adonis2(
      primary_dist ~ Group,
      data = meta,
      permutations = n_perm,
      strata = strata_vec
    ))

    perm_df <- as.data.frame(perm_res)
    perm_df$Term <- rownames(perm_df)
    perm_df <- perm_df[, c("Term", setdiff(colnames(perm_df), "Term"))]

    write.table(perm_df, permanova_file, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, permanova_file)

    # Betadisper
    disp_res <- vegan::betadisper(primary_dist, meta$Group)
    disp_perm <- vegan::permutest(disp_res, permutations = n_perm)

    disp_df <- data.frame(
      Analysis = "Betadisper (Homogeneity of Multivariate Dispersions)",
      F_Statistic = disp_perm$tab$F[1],
      P_Value = disp_perm$tab$`Pr(>F)`[1],
      Permutations = n_perm,
      Warning = if (any(group_counts < 3)) "Sample size < 3 in at least one group; low statistical power" else "None",
      stringsAsFactors = FALSE
    )
    write.table(disp_df, betadisper_file, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, betadisper_file)
  } else {
    skip_perm <- data.frame(
      Status = "Skipped",
      Reason = sprintf(
        "PERMANOVA requires non-zero distances, residual degrees of freedom, and at least 2 groups with >= 2 samples each. Found: %s",
        paste(sprintf("%s (n=%d)", names(group_counts), as.integer(group_counts)), collapse = ", ")
      ),
      stringsAsFactors = FALSE
    )
    write.table(skip_perm, permanova_file, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(skip_perm, betadisper_file, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, permanova_file, betadisper_file)
  }

  list(
    status = "completed",
    outputs = all_outputs
  )
}
