# =============================================================================
# Module 05: Ordination (Hellinger PCA, Bray-Curtis NMDS)
# =============================================================================

suppressMessages({
  library(vegan)
  library(dplyr)
  library(ggplot2)
})

run_ordination <- function(context) {
  cfg <- context$config
  ord_dir <- cfg$output$dirs$ordination
  dir.create(ord_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_outputs <- character(0)
  
  # Cohort gate
  if (context$mode != "cohort" || length(context$samples) < 2) {
    skip_file <- file.path(ord_dir, "ordination_skipped.tsv")
    skip_df <- data.frame(
      Status = "Skipped",
      Reason = sprintf("Ordination is cohort-only. Current mode is '%s' with %d sample(s).",
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
  
  otu_table <- t(class_counts)
  sample_sums <- rowSums(otu_table)
  rel_otu <- sweep(otu_table, 1, sample_sums, "/")
  
  # 1. Hellinger PCA
  hel_otu <- vegan::decostand(rel_otu, method = "hellinger")
  pca_res <- vegan::rda(hel_otu)
  
  eig <- pca_res$CA$eig
  var_exp <- round(100 * eig / sum(eig), 2)
  var_df <- data.frame(
    PC = paste0("PC", seq_along(var_exp)),
    Eigenvalue = as.numeric(eig),
    VarianceExplained = as.numeric(var_exp),
    CumulativeVariance = cumsum(as.numeric(var_exp)),
    stringsAsFactors = FALSE
  )
  var_tsv <- file.path(ord_dir, "pca_variance.tsv")
  write.table(var_df, var_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, var_tsv)
  
  scores_mat <- as.matrix(scores(pca_res, display = "sites"))
  scores_df <- data.frame(
    SampleID = rownames(scores_mat),
    PC1 = if (ncol(scores_mat) >= 1) scores_mat[, 1] else rep(0, nrow(scores_mat)),
    PC2 = if (ncol(scores_mat) >= 2) scores_mat[, 2] else rep(0, nrow(scores_mat)),
    stringsAsFactors = FALSE
  )
  if (!is.null(meta)) {
    scores_df <- merge(meta, scores_df, by = "SampleID", sort = FALSE)
  }
  scores_tsv <- file.path(ord_dir, "pca_scores.tsv")
  write.table(scores_df, scores_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, scores_tsv)
  
  # Species loadings
  loadings_mat <- as.matrix(scores(pca_res, display = "species"))
  loadings_df <- data.frame(
    Taxon = rownames(loadings_mat),
    PC1 = if (ncol(loadings_mat) >= 1) loadings_mat[, 1] else rep(0, nrow(loadings_mat)),
    PC2 = if (ncol(loadings_mat) >= 2) loadings_mat[, 2] else rep(0, nrow(loadings_mat)),
    stringsAsFactors = FALSE
  )
  load_tsv <- file.path(ord_dir, "pca_loadings.tsv")
  write.table(loadings_df, load_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, load_tsv)
  
  # PCA plot (when at least 2 PCs available)
  if (ncol(scores_mat) >= 2) {
    p_pca <- ggplot(scores_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3.5, alpha = 0.85) +
      labs(
        title = "Hellinger PCA (Classified Relative Abundance)",
        x = sprintf("PC1 (%.1f%%)", var_exp[1]),
        y = sprintf("PC2 (%.1f%%)", var_exp[2])
      ) +
      theme_amplicon()
    
    pca_plot_path <- file.path(ord_dir, "05a_pca_plot.png")
    save_plot(pca_plot_path, p_pca, width = 7, height = 5.5)
    all_outputs <- c(all_outputs, pca_plot_path)
  }
  
  # 2. Bray-Curtis NMDS (requires at least 3 non-identical samples)
  if (length(samples) >= 3) {
    set.seed(seed)
    nmds_res <- suppressWarnings(tryCatch({
      vegan::metaMDS(rel_otu, distance = "bray", k = 2, trymax = 50, trace = 0)
    }, error = function(e) NULL))
    
    if (!is.null(nmds_res) && !is.null(nmds_res$points)) {
      nmds_df <- data.frame(
        SampleID = rownames(nmds_res$points),
        NMDS1 = nmds_res$points[, 1],
        NMDS2 = nmds_res$points[, 2],
        Stress = nmds_res$stress,
        Converged = nmds_res$converged,
        stringsAsFactors = FALSE
      )
      if (!is.null(meta)) {
        nmds_df <- merge(meta, nmds_df, by = "SampleID", sort = FALSE)
      }
      
      nmds_tsv <- file.path(ord_dir, "nmds_scores.tsv")
      write.table(nmds_df, nmds_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
      all_outputs <- c(all_outputs, nmds_tsv)
      
      p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group)) +
        geom_point(size = 3.5, alpha = 0.85) +
        labs(
          title = "Bray-Curtis NMDS Ordination",
          subtitle = sprintf("Stress: %.3f %s (k=2, seed=%d)",
                             nmds_res$stress,
                             if (nmds_res$stress > 0.2) "(Caution: High Stress)" else "",
                             seed),
          x = "NMDS1", y = "NMDS2"
        ) +
        theme_amplicon()
      
      nmds_plot_path <- file.path(ord_dir, "05b_nmds_plot.png")
      save_plot(nmds_plot_path, p_nmds, width = 7, height = 5.5)
      all_outputs <- c(all_outputs, nmds_plot_path)
    }
  }
  
  list(
    status = "completed",
    outputs = all_outputs
  )
}
