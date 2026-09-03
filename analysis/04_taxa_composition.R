# =============================================================================
# Module 04: Taxonomic Composition Analysis
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  library(pheatmap)
})

RANKS_TO_ANALYZE <- c("phylum", "class", "order", "family", "genus", "species")

run_taxa_composition <- function(context) {
  cfg <- context$config
  comp_dir <- cfg$output$dirs$composition
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_outputs <- character(0)
  
  top_n <- cfg$composition$top_n_taxa %||% 15L
  heatmap_rank <- cfg$composition$heatmap_rank %||% "genus"
  heatmap_transform <- cfg$composition$heatmap_transform %||% "log10_relative"
  
  samples <- context$samples
  unclass_idx <- context$unclass_index
  count_matrix <- context$count_matrix
  taxonomy_df <- context$taxonomy
  
  # Total and classified read denominators per sample
  sample_totals <- colSums(count_matrix)
  unclass_counts <- count_matrix[unclass_idx, ]
  class_totals <- sample_totals - unclass_counts
  
  # Classified rows
  class_counts <- count_matrix[-unclass_idx, , drop = FALSE]
  class_tax <- taxonomy_df[-unclass_idx, , drop = FALSE]
  
  # Analyze each rank
  rank_tables <- list()
  
  for (rk in RANKS_TO_ANALYZE) {
    rk_idx <- which(colnames(taxonomy_df) == rk)
    
    # Prefix-based aggregation up to rank rk
    rk_prefixes <- apply(class_tax[, 1:rk_idx, drop = FALSE], 1, paste, collapse = ";")
    leaf_names <- class_tax[[rk]]
    
    # Contextualize ambiguous labels (e.g. "Unknown" gets parent context)
    contextualized_names <- ifelse(
      leaf_names %in% c("Unknown", "unclassified", "uncultured"),
      sprintf("%s (%s)", leaf_names, class_tax[[rk_idx - 1]]),
      leaf_names
    )
    
    # Aggregate counts by prefix
    agg_df <- data.frame(
      TaxonPath = rk_prefixes,
      Taxon = contextualized_names,
      class_counts,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      group_by(TaxonPath, Taxon) %>%
      summarise(across(all_of(samples), sum), .groups = "drop")
    
    # Calculate relative abundances (classified-only denominator)
    rel_df <- agg_df
    for (s in samples) {
      denom <- class_totals[s]
      rel_df[[s]] <- if (denom > 0) rel_df[[s]] / denom else 0
    }
    
    # Write TSVs
    count_file <- file.path(comp_dir, sprintf("count_%s.tsv", rk))
    rel_file <- file.path(comp_dir, sprintf("rel_abundance_%s.tsv", rk))
    
    write.table(agg_df, count_file, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(rel_df, rel_file, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, count_file, rel_file)
    
    rank_tables[[rk]] <- list(counts = agg_df, rel = rel_df)
  }
  
  # Single-Sample Bar Plots
  if (length(samples) == 1) {
    s_col <- samples[1]
    
    # 1. Phylum bar plot
    phylum_rel <- rank_tables[["phylum"]]$rel %>%
      select(Taxon, all_of(s_col)) %>%
      rename(rel = all_of(s_col)) %>%
      arrange(desc(rel))
    
    n_keep_phylum <- min(7L, nrow(phylum_rel))
    phylum_top <- phylum_rel %>%
      mutate(
        DisplayTaxon = if_else(row_number() <= n_keep_phylum, Taxon, "Other phyla")
      ) %>%
      group_by(DisplayTaxon) %>%
      summarise(rel = sum(rel), .groups = "drop") %>%
      arrange(desc(rel))
    
    phylum_top$DisplayTaxon <- factor(phylum_top$DisplayTaxon, levels = rev(phylum_top$DisplayTaxon))
    
    p_phylum <- ggplot(phylum_top, aes(x = DisplayTaxon, y = rel, fill = DisplayTaxon)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = sprintf("%.1f%%", 100 * rel)), hjust = -0.15, size = 3.3) +
      scale_fill_manual(values = get_phylum_colors(levels(phylum_top$DisplayTaxon))) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.18))) +
      coord_flip() +
      labs(
        title = sprintf("Phylum-Level Composition (%s)", s_col),
        subtitle = sprintf("Relative to %s classified reads", scales::comma(class_totals[s_col])),
        x = NULL, y = "Relative Abundance"
      ) +
      theme_amplicon() +
      theme(legend.position = "none")
    
    p_phylum_path <- file.path(comp_dir, "04a_phylum_composition.png")
    save_plot(p_phylum_path, p_phylum, width = 7.5, height = 5)
    all_outputs <- c(all_outputs, p_phylum_path)
    
    # 2. Top Family bar plot
    family_counts <- rank_tables[["family"]]$counts %>%
      select(Taxon, all_of(s_col)) %>%
      rename(count = all_of(s_col)) %>%
      arrange(desc(count)) %>%
      slice_head(n = top_n)
    family_counts$Taxon <- factor(family_counts$Taxon, levels = rev(family_counts$Taxon))
    
    p_family <- ggplot(family_counts, aes(x = Taxon, y = count)) +
      geom_col(width = 0.7, fill = "#1b9e77") +
      geom_text(aes(label = scales::comma(count)), hjust = -0.15, size = 3) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
      coord_flip() +
      labs(
        title = sprintf("Top %d Families by Read Count (%s)", top_n, s_col),
        x = NULL, y = "Read count"
      ) +
      theme_amplicon()
    
    p_family_path <- file.path(comp_dir, "04b_family_composition.png")
    save_plot(p_family_path, p_family, width = 8.5, height = 6)
    all_outputs <- c(all_outputs, p_family_path)
    
    # 3. Top Genus bar plot
    genus_counts <- rank_tables[["genus"]]$counts %>%
      select(Taxon, all_of(s_col)) %>%
      rename(count = all_of(s_col)) %>%
      arrange(desc(count)) %>%
      slice_head(n = top_n)
    genus_counts$Taxon <- factor(genus_counts$Taxon, levels = rev(genus_counts$Taxon))
    
    p_genus <- ggplot(genus_counts, aes(x = Taxon, y = count)) +
      geom_col(width = 0.7, fill = "#d95f02") +
      geom_text(aes(label = scales::comma(count)), hjust = -0.15, size = 3) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
      coord_flip() +
      labs(
        title = sprintf("Top %d Genera by Read Count (%s)", top_n, s_col),
        x = NULL, y = "Read count"
      ) +
      theme_amplicon()
    
    p_genus_path <- file.path(comp_dir, "04c_genus_composition.png")
    save_plot(p_genus_path, p_genus, width = 8.5, height = 6)
    all_outputs <- c(all_outputs, p_genus_path)
    
    # 4. Top Species bar plot
    species_rel <- rank_tables[["species"]]$rel %>%
      select(Taxon, all_of(s_col)) %>%
      rename(rel = all_of(s_col)) %>%
      arrange(desc(rel)) %>%
      slice_head(n = top_n)
    species_rel$Taxon <- factor(species_rel$Taxon, levels = rev(species_rel$Taxon))
    
    p_species <- ggplot(species_rel, aes(x = Taxon, y = rel)) +
      geom_col(width = 0.7, fill = "#7570b3") +
      geom_text(aes(label = sprintf("%.2f%%", 100 * rel)), hjust = -0.1, size = 3) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.25))) +
      coord_flip() +
      labs(
        title = sprintf("Top %d Species by Relative Abundance (%s)", top_n, s_col),
        x = NULL, y = "Relative Abundance (Classified reads)"
      ) +
      theme_amplicon() +
      theme(axis.text.y = element_text(face = "italic"))
    
    p_species_path <- file.path(comp_dir, "04d_species_composition.png")
    save_plot(p_species_path, p_species, width = 9, height = 6)
    all_outputs <- c(all_outputs, p_species_path)
  }
  
  # Cohort Stacked Bars & Heatmap
  if (length(samples) >= 2) {
    # 1. Phylum Stacked Bar
    phylum_long <- rank_tables[["phylum"]]$rel %>%
      tidyr::pivot_longer(cols = all_of(samples), names_to = "SampleID", values_to = "rel") %>%
      group_by(Taxon) %>%
      mutate(mean_rel = mean(rel)) %>%
      ungroup()
    
    top_phyla <- phylum_long %>%
      distinct(Taxon, mean_rel) %>%
      arrange(desc(mean_rel)) %>%
      slice_head(n = 8) %>%
      pull(Taxon)
    
    phylum_stacked <- phylum_long %>%
      mutate(DisplayTaxon = if_else(Taxon %in% top_phyla, Taxon, "Other")) %>%
      group_by(SampleID, DisplayTaxon) %>%
      summarise(rel = sum(rel), .groups = "drop")
    
    phylum_levels <- c(setdiff(unique(phylum_stacked$DisplayTaxon), "Other"), "Other")
    phylum_stacked$DisplayTaxon <- factor(phylum_stacked$DisplayTaxon, levels = rev(phylum_levels))
    
    p_stack <- ggplot(phylum_stacked, aes(x = SampleID, y = rel, fill = DisplayTaxon)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = get_phylum_colors(levels(phylum_stacked$DisplayTaxon)), name = "Phylum") +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
      labs(
        title = "Phylum-Level Community Composition Across Samples",
        x = NULL, y = "Relative Abundance"
      ) +
      theme_amplicon() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
    
    p_stack_path <- file.path(comp_dir, "04_phylum_stacked.png")
    save_plot(p_stack_path, p_stack, width = 8.5, height = 6.5)
    all_outputs <- c(all_outputs, p_stack_path)
    
    # 2. Heatmap
    if (heatmap_rank %in% names(rank_tables)) {
      rk_data <- rank_tables[[heatmap_rank]]$rel
      
      # Select top N taxa by mean relative abundance
      mat_data <- as.matrix(rk_data[, samples, drop = FALSE])
      rownames(mat_data) <- rk_data$Taxon
      mean_rel <- rowMeans(mat_data)
      
      top_idx <- order(mean_rel, decreasing = TRUE)[seq_len(min(top_n, nrow(mat_data)))]
      mat_top <- mat_data[top_idx, , drop = FALSE]
      
      # Apply transform
      mat_transformed <- if (heatmap_transform == "log10_relative") {
        log10(mat_top + 1e-4)
      } else {
        mat_top
      }
      
      # Annotate columns with metadata if available
      anno_col <- NA
      if (!is.null(context$metadata) && "Group" %in% colnames(context$metadata)) {
        anno_df <- data.frame(Group = context$metadata$Group, row.names = context$metadata$SampleID)
        anno_col <- anno_df[samples, , drop = FALSE]
      }
      
      heatmap_path <- file.path(comp_dir, sprintf("04_heatmap_%s.png", heatmap_rank))
      
      # Save pheatmap
      png(heatmap_path, width = 8, height = 7, units = "in", res = 150)
      pheatmap::pheatmap(
        mat_transformed,
        annotation_col = if (is.data.frame(anno_col)) anno_col else NA,
        color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
        main = sprintf("Top %d %s Relative Abundance (%s)", nrow(mat_top), heatmap_rank, heatmap_transform),
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        fontsize = 9
      )
      dev.off()
      all_outputs <- c(all_outputs, heatmap_path)
    }
  }
  
  list(
    status = "completed",
    outputs = all_outputs
  )
}
