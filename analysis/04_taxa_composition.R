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
OTHER_TAXON_PATH <- "__OTHER__"

taxon_path_parts <- function(path) {
  parts <- strsplit(as.character(path), ";", fixed = TRUE)[[1]]
  if (!length(parts) || anyNA(parts) || any(!nzchar(parts))) {
    stop("TaxonPath values must contain at least one non-empty component.", call. = FALSE)
  }
  parts
}

make_taxon_display_map <- function(rank_table) {
  if (!is.data.frame(rank_table) || !"TaxonPath" %in% names(rank_table)) {
    stop("A rank table with a TaxonPath column is required.", call. = FALSE)
  }
  input_paths <- as.character(rank_table$TaxonPath)
  if (!length(input_paths) || anyNA(input_paths) || any(!nzchar(input_paths)) ||
      anyDuplicated(input_paths)) {
    stop("TaxonPath values must be unique, non-empty strings.", call. = FALSE)
  }

  # Build from an immutable sorted path set. Display labels are presentation
  # only; sorting here makes collision resolution independent of input-row
  # order while TaxonPath remains the stable join key.
  paths <- sort(input_paths, method = "radix")
  parts <- lapply(paths, taxon_path_parts)
  leaves <- vapply(parts, function(x) x[length(x)], character(1))
  parents <- vapply(parts, function(x) if (length(x) > 1L) x[length(x) - 1L] else "", character(1))
  base_display <- leaves

  unknown <- tolower(leaves) %in% c("unknown", "unclassified", "uncultured")
  base_display[unknown] <- ifelse(
    nzchar(parents[unknown]),
    sprintf("%s [%s]", leaves[unknown], parents[unknown]),
    sprintf("%s [%s]", leaves[unknown], paths[unknown])
  )
  biological_other <- leaves == "Other"
  base_display[biological_other] <- ifelse(
    nzchar(parents[biological_other]),
    sprintf("Other [%s]", parents[biological_other]),
    sprintf("Other [%s]", paths[biological_other])
  )

  candidate_labels <- function(depth) {
    vapply(seq_along(paths), function(i) {
      if (depth[i] <= 1L && nzchar(base_display[i])) return(base_display[i])
      width <- min(as.integer(depth[i]), length(parts[[i]]))
      suffix_parts <- head(tail(parts[[i]], width), -1L)
      if (!length(suffix_parts)) return(sprintf("%s [%s]", leaves[i], paths[i]))
      sprintf("%s [%s]", leaves[i], paste(suffix_parts, collapse = " / "))
    }, character(1))
  }

  depth <- rep(1L, length(paths))
  repeat {
    display <- candidate_labels(depth)
    collision <- duplicated(display) | duplicated(display, fromLast = TRUE) |
      display == "Other"
    if (!any(collision)) break
    can_expand <- collision & depth < vapply(parts, length, integer(1))
    if (!any(can_expand)) break
    depth[can_expand] <- depth[can_expand] + 1L
  }

  # If the full lineage still collides with a literal biological label or a
  # punctuation-shaped generated label, use the canonical path as the final
  # deterministic suffix. Paths are unique, so this fallback is stable.
  final_collision <- duplicated(display) | duplicated(display, fromLast = TRUE) |
    display == "Other"
  if (any(final_collision)) {
    display[final_collision] <- sprintf("%s [%s]", display[final_collision], paths[final_collision])
  }
  if (anyDuplicated(display) || any(display == "Other")) {
    display <- sprintf("%s {%s}", display, paths)
  }
  if (anyDuplicated(display) || any(display == "Other")) {
    stop("Could not derive unique deterministic composition display labels.", call. = FALSE)
  }
  result <- data.frame(TaxonPath = paths, DisplayTaxon = display, stringsAsFactors = FALSE)
  result[match(input_paths, result$TaxonPath), , drop = FALSE]
}

rank_relative_matrix <- function(rel_table, samples) {
  if (!is.data.frame(rel_table) || !"TaxonPath" %in% names(rel_table)) {
    stop("A relative-abundance rank table with TaxonPath is required.", call. = FALSE)
  }
  missing_samples <- setdiff(samples, names(rel_table))
  if (length(missing_samples)) {
    stop(sprintf("Relative-abundance table is missing sample column(s): %s",
                 paste(missing_samples, collapse = ", ")), call. = FALSE)
  }
  if (anyDuplicated(rel_table$TaxonPath)) {
    stop("Relative-abundance rank tables must have unique TaxonPath values.", call. = FALSE)
  }
  values <- as.matrix(rel_table[, samples, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "numeric")
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("Relative-abundance values must be finite and non-negative.", call. = FALSE)
  }
  rownames(values) <- as.character(rel_table$TaxonPath)
  values
}

select_display_taxa <- function(rel_table, samples, valid_samples = samples,
                                min_mean = 0, min_n = 1L, max_n = 15L) {
  values <- rank_relative_matrix(rel_table, samples)
  valid_samples <- intersect(samples, as.character(valid_samples))
  if (!length(valid_samples)) {
    return(data.frame(
      TaxonPath = character(0), MeanRelativeAbundance = numeric(0),
      SelectionOrder = integer(0), stringsAsFactors = FALSE
    ))
  }

  valid_values <- values[, valid_samples, drop = FALSE]
  positive <- rowSums(valid_values > 0) > 0
  candidate_paths <- rownames(values)[positive]
  if (!length(candidate_paths)) {
    return(data.frame(
      TaxonPath = character(0), MeanRelativeAbundance = numeric(0),
      SelectionOrder = integer(0), stringsAsFactors = FALSE
    ))
  }

  means <- rowMeans(valid_values[candidate_paths, , drop = FALSE])
  candidates <- data.frame(
    TaxonPath = candidate_paths,
    MeanRelativeAbundance = as.numeric(means),
    stringsAsFactors = FALSE
  )
  candidates <- candidates[order(-candidates$MeanRelativeAbundance,
                                 candidates$TaxonPath, method = "radix"), , drop = FALSE]

  eligible <- candidates[candidates$MeanRelativeAbundance >= min_mean, , drop = FALSE]
  selected_n <- min(as.integer(max_n), nrow(eligible))
  selected <- if (selected_n > 0L) eligible[seq_len(selected_n), , drop = FALSE] else eligible[FALSE, , drop = FALSE]
  minimum_n <- min(as.integer(min_n), nrow(candidates), as.integer(max_n))
  if (nrow(selected) < minimum_n) {
    remaining <- candidates[!candidates$TaxonPath %in% selected$TaxonPath, , drop = FALSE]
    fill_n <- min(minimum_n - nrow(selected), nrow(remaining))
    if (fill_n > 0L) selected <- rbind(selected, remaining[seq_len(fill_n), , drop = FALSE])
  }
  selected <- selected[order(-selected$MeanRelativeAbundance,
                             selected$TaxonPath, method = "radix"), , drop = FALSE]
  selected$SelectionOrder <- seq_len(nrow(selected))
  rownames(selected) <- NULL
  selected
}

plot_metadata <- function(samples, metadata = NULL, valid_samples = samples) {
  if (is.null(metadata)) {
    result <- data.frame(
      SampleID = samples,
      Group = rep("Single", length(samples)),
      stringsAsFactors = FALSE
    )
  } else {
    if (!is.data.frame(metadata) || !all(c("SampleID", "Group") %in% names(metadata))) {
      stop("Plot metadata must contain SampleID and Group columns.", call. = FALSE)
    }
    index <- match(samples, metadata$SampleID)
    if (anyNA(index)) stop("Plot metadata does not cover every sample.", call. = FALSE)
    result <- data.frame(
      SampleID = samples,
      Group = as.character(metadata$Group[index]),
      stringsAsFactors = FALSE
    )
  }
  if ("ValidDenominator" %in% names(metadata %||% data.frame())) {
    result$ValidDenominator <- as.logical(metadata$ValidDenominator[match(samples, metadata$SampleID)])
  } else {
    result$ValidDenominator <- samples %in% as.character(valid_samples)
  }
  if (anyNA(result$ValidDenominator)) result$ValidDenominator <- FALSE
  result$GroupOrder <- match(result$Group, unique(result$Group))
  result$SampleOrder <- match(result$SampleID, samples)
  result
}

ordered_cohort_samples <- function(samples, metadata) {
  meta <- plot_metadata(samples, metadata, valid_samples = samples)
  group_order <- unique(meta$Group)
  sample_order <- unlist(lapply(group_order, function(group) {
    meta$SampleID[meta$Group == group]
  }), use.names = FALSE)
  list(sample_order = sample_order, group_order = group_order)
}

collapse_rank_for_display <- function(rel_table, selected_paths, samples, metadata = NULL) {
  if (is.data.frame(selected_paths)) selected_paths <- selected_paths$TaxonPath
  selected_paths <- as.character(selected_paths)
  values <- rank_relative_matrix(rel_table, samples)
  if (length(selected_paths) && any(!selected_paths %in% rownames(values))) {
    stop("Selected TaxonPath values are absent from the rank table.", call. = FALSE)
  }

  valid_samples <- samples
  if (!is.null(metadata) && "ValidDenominator" %in% names(metadata)) {
    valid_samples <- metadata$SampleID[isTRUE(metadata$ValidDenominator) |
      as.logical(metadata$ValidDenominator)]
  }
  meta <- plot_metadata(samples, metadata, valid_samples = valid_samples)
  valid <- meta$ValidDenominator
  sample_order <- if (!is.null(metadata) && length(samples) > 1L) {
    ordered_cohort_samples(samples, metadata)$sample_order
  } else {
    samples
  }

  map <- make_taxon_display_map(rel_table)
  display_values <- values
  display_values[, !valid] <- 0
  mean_values <- if (any(valid) && length(selected_paths)) {
    rowMeans(values[selected_paths, samples[valid], drop = FALSE])
  } else {
    stats::setNames(numeric(length(selected_paths)), selected_paths)
  }
  if (length(selected_paths)) names(mean_values) <- selected_paths

  selected_sum <- if (length(selected_paths)) {
    colSums(display_values[selected_paths, samples, drop = FALSE])
  } else {
    stats::setNames(rep(0, length(samples)), samples)
  }
  other_values <- stats::setNames(pmax(0, 1 - selected_sum), samples)
  other_values[!valid] <- 0
  has_other <- any(other_values > 0)
  if (has_other) {
    mean_values <- c(mean_values, stats::setNames(
      if (any(valid)) mean(other_values[valid]) else NA_real_, OTHER_TAXON_PATH
    ))
  }
  path_order <- c(selected_paths, if (has_other) OTHER_TAXON_PATH else character(0))
  labels <- if (length(path_order)) {
    ifelse(path_order == OTHER_TAXON_PATH, "Other",
           map$DisplayTaxon[match(path_order, map$TaxonPath)])
  } else character(0)

  rows <- lapply(sample_order, function(sample_id) {
    sample_index <- match(sample_id, samples)
    values_for_sample <- if (length(selected_paths)) {
      as.numeric(display_values[selected_paths, sample_id])
    } else numeric(0)
    if (has_other) values_for_sample <- c(values_for_sample, other_values[sample_id])
    data.frame(
      TaxonPath = path_order,
      DisplayTaxon = labels,
      SampleID = sample_id,
      Group = meta$Group[sample_index],
      RelativeAbundance = values_for_sample,
      MeanRelativeAbundance = as.numeric(mean_values[path_order]),
      IsOther = path_order == OTHER_TAXON_PATH,
      ValidDenominator = valid[sample_index],
      StackOrder = seq_along(path_order),
      SampleOrder = match(sample_id, sample_order),
      GroupOrder = meta$GroupOrder[sample_index],
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

summarize_group_means <- function(display_long) {
  required <- c("Group", "GroupOrder", "SampleID", "TaxonPath", "DisplayTaxon",
                "RelativeAbundance", "IsOther", "StackOrder", "ValidDenominator")
  if (!is.data.frame(display_long) || length(setdiff(required, names(display_long)))) {
    stop("Display-long data is missing group-mean columns.", call. = FALSE)
  }
  sample_info <- unique(display_long[, c("Group", "GroupOrder", "SampleID", "ValidDenominator")])
  groups <- sample_info[!duplicated(sample_info$Group), c("Group", "GroupOrder"), drop = FALSE]
  taxa <- display_long[order(display_long$StackOrder),
                       c("TaxonPath", "DisplayTaxon", "IsOther", "StackOrder"), drop = FALSE]
  taxa <- taxa[!duplicated(taxa$TaxonPath), , drop = FALSE]
  rows <- list()
  for (g in groups$Group) {
    g_samples <- sample_info$SampleID[sample_info$Group == g]
    g_valid <- sample_info$ValidDenominator[sample_info$Group == g]
    for (i in seq_len(nrow(taxa))) {
      values <- display_long$RelativeAbundance[
        display_long$Group == g & display_long$TaxonPath == taxa$TaxonPath[i] &
          display_long$ValidDenominator
      ]
      rows[[length(rows) + 1L]] <- data.frame(
        Group = g,
        TaxonPath = taxa$TaxonPath[i],
        DisplayTaxon = taxa$DisplayTaxon[i],
        MeanRelativeAbundance = if (length(values)) mean(values) else NA_real_,
        SamplesTotal = length(g_samples),
        SamplesUsed = sum(g_valid),
        SamplesExcludedZeroClassified = sum(!g_valid),
        IsOther = taxa$IsOther[i],
        StackOrder = taxa$StackOrder[i],
        GroupOrder = groups$GroupOrder[match(g, groups$Group)],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  result <- do.call(rbind, rows)
  result[order(result$GroupOrder, result$StackOrder, result$TaxonPath), , drop = FALSE]
}

build_stacked_taxa_plot <- function(display_long, colors, sample_order,
                                    group_order = NULL, single = TRUE,
                                    rank = "taxonomic", subtitle = NULL,
                                    invalid_groups = character(0)) {
  if (!nrow(display_long)) stop("Cannot plot an empty composition table.", call. = FALSE)
  stack_levels <- display_long[order(display_long$StackOrder), "DisplayTaxon"]
  stack_levels <- unique(as.character(stack_levels))
  data <- display_long
  data$DisplayTaxon <- factor(data$DisplayTaxon, levels = stack_levels)
  data$SampleID <- factor(data$SampleID, levels = sample_order)
  if (!single) data$Group <- factor(data$Group, levels = group_order)

  plot_data <- data[is.finite(data$RelativeAbundance), , drop = FALSE]
  p <- ggplot(plot_data, aes(x = SampleID, y = RelativeAbundance, fill = DisplayTaxon)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = colors, drop = FALSE, name = "Taxon") +
    scale_y_continuous(
      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = function(x) sprintf("%d%%", round(100 * x)),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = sprintf("%s-level classified-read composition", tools::toTitleCase(rank)),
      subtitle = subtitle %||%
        "Relative abundance among classified reads; Other is the residual",
      x = NULL, y = "Relative abundance"
    ) +
    theme_amplicon() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      axis.text.x = element_text(angle = if (single) 0 else 90,
                                 hjust = if (single) 0.5 else 1,
                                 vjust = if (single) 0.5 else 0.5)
    )

  if (!single) {
    p <- p + facet_grid(cols = vars(Group), scales = "free_x", space = "free_x")
  } else {
    labels <- data[data$ValidDenominator & data$RelativeAbundance >= 0.03, , drop = FALSE]
    if (nrow(labels)) {
      p <- p + geom_text(
        data = labels,
        aes(label = sprintf("%.1f%%", 100 * RelativeAbundance)),
        position = position_stack(vjust = 0.5, reverse = TRUE),
        size = 3
      )
    }
  }
  if (length(invalid_groups)) {
    invalid_labels <- data.frame(
      SampleID = factor(as.character(invalid_groups), levels = sample_order),
      Label = "Undefined: 0 valid samples",
      stringsAsFactors = FALSE
    )
    p <- p + geom_text(
      data = invalid_labels,
      aes(x = SampleID, y = 0.02, label = Label),
      inherit.aes = FALSE, colour = "#666666", size = 3
    )
  }
  if (length(stack_levels) > 10L) p <- p + guides(fill = guide_legend(ncol = 2))
  p
}

prepare_heatmap_data <- function(rel_table, samples, metadata = NULL, top_n = 10L,
                                 include_other = TRUE, transform = "log10_relative",
                                 rank = attr(rel_table, "rank") %||% NA_character_) {
  values <- rank_relative_matrix(rel_table, samples)
  valid_samples <- samples
  if (!is.null(metadata) && "ValidDenominator" %in% names(metadata)) {
    valid_samples <- metadata$SampleID[as.logical(metadata$ValidDenominator)]
  }
  valid_samples <- intersect(samples, as.character(valid_samples))
  meta <- plot_metadata(samples, metadata, valid_samples = valid_samples)
  valid_by_sample <- stats::setNames(as.logical(meta$ValidDenominator), meta$SampleID)
  valid <- unname(valid_by_sample[samples])

  empty <- function(reason) list(
    skipped = TRUE,
    skip = data.frame(Rank = rank, Reason = reason, stringsAsFactors = FALSE)
  )
  if (!length(valid_samples) || !any(valid)) {
    return(empty("No sample has a positive classified-read denominator."))
  }

  valid_values <- values[, valid_samples, drop = FALSE]
  means <- rowMeans(valid_values)
  positive <- rowSums(valid_values > 0) > 0
  if (!any(positive)) return(empty("No positive classified relative abundance exists at this rank."))
  candidate_paths <- rownames(values)[positive]
  ordered <- candidate_paths[order(-means[candidate_paths], candidate_paths, method = "radix")]
  selected_paths <- head(ordered, as.integer(top_n))

  selected_values <- values[selected_paths, samples, drop = FALSE]
  if (any(!valid)) selected_values[, !valid] <- 0
  residual <- stats::setNames(pmax(0, 1 - colSums(selected_values)), samples)
  residual[!valid] <- 0
  display_paths <- selected_paths
  if (isTRUE(include_other) && any(residual > 0)) {
    display_paths <- c(display_paths, OTHER_TAXON_PATH)
  }
  raw <- selected_values
  if (length(display_paths) > length(selected_paths)) {
    raw <- rbind(raw, residual)
  }
  rownames(raw) <- display_paths

  pseudo_count <- NA_real_
  if (identical(transform, "log10_relative")) {
    transform_source <- raw[, valid, drop = FALSE]
    positive_values <- as.numeric(transform_source[transform_source > 0])
    if (!length(positive_values)) return(empty("No positive value is available for the log10 transform."))
    pseudo_count <- min(positive_values) / 2
    transformed <- log10(raw + pseudo_count)
  } else if (identical(transform, "none")) {
    transformed <- raw
  } else {
    stop("Unsupported heatmap transform.", call. = FALSE)
  }
  if (any(!is.finite(transformed))) stop("Heatmap transformed values must be finite.", call. = FALSE)

  cohort <- !is.null(metadata) && length(samples) >= 2L
  sample_order <- if (cohort) ordered_cohort_samples(samples, metadata)$sample_order else samples
  raw <- raw[, sample_order, drop = FALSE]
  transformed <- transformed[, sample_order, drop = FALSE]
  valid_ordered <- unname(valid_by_sample[sample_order])

  display_map <- make_taxon_display_map(rel_table)
  display_labels <- ifelse(
    display_paths == OTHER_TAXON_PATH, "Other",
    display_map$DisplayTaxon[match(display_paths, display_map$TaxonPath)]
  )
  rownames(raw) <- display_labels
  rownames(transformed) <- display_labels

  row_hclust <- NULL
  row_order <- seq_len(nrow(transformed))
  if (cohort && nrow(transformed) >= 2L) {
    cluster_values <- transformed[, valid_ordered, drop = FALSE]
    row_hclust <- stats::hclust(stats::dist(cluster_values, method = "euclidean"), method = "complete")
    row_order <- row_hclust$order
  } else {
    row_order <- order(-rowMeans(raw), display_paths, method = "radix")
    raw <- raw[row_order, , drop = FALSE]
    transformed <- transformed[row_order, , drop = FALSE]
    display_paths <- display_paths[row_order]
    display_labels <- display_labels[row_order]
  }

  matrix_plot <- transformed
  if (any(!valid_ordered)) matrix_plot[, !valid_ordered] <- NA_real_

  rows <- list()
  for (i in seq_along(display_paths)) {
    original_index <- if (cohort && nrow(transformed) >= 2L) row_order[i] else i
    path <- display_paths[original_index]
    label <- display_labels[original_index]
    for (j in seq_along(sample_order)) {
      sample_id <- sample_order[j]
      sample_meta <- meta[match(sample_id, meta$SampleID), , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        Rank = rank,
        TaxonPath = path,
        DisplayTaxon = label,
        SampleID = sample_id,
        Group = sample_meta$Group,
        RelativeAbundance = as.numeric(raw[original_index, sample_id]),
        Transform = transform,
        PseudoCount = as.numeric(pseudo_count),
        TransformedValue = as.numeric(transformed[original_index, sample_id]),
        IsOther = identical(path, OTHER_TAXON_PATH),
        RowOrder = i,
        ColumnOrder = j,
        ValidDenominator = sample_meta$ValidDenominator,
        stringsAsFactors = FALSE
      )
    }
  }
  sidecar <- do.call(rbind, rows)
  list(
    skipped = FALSE,
    matrix_raw = raw,
    matrix_transformed = transformed,
    matrix_plot = matrix_plot,
    sidecar = sidecar,
    row_hclust = row_hclust,
    sample_order = sample_order,
    group_order = if (cohort) unique(meta$Group[match(sample_order, meta$SampleID)]) else "Single",
    cohort = cohort,
    transform = transform,
    pseudo_count = pseudo_count
  )
}

draw_taxa_heatmap <- function(prepared, path, mode = if (isTRUE(prepared$cohort)) "cohort" else "single",
                              rank = unique(prepared$sidecar$Rank)) {
  if (isTRUE(prepared$skipped)) stop("Cannot draw a skipped heatmap.", call. = FALSE)
  matrix_values <- prepared$matrix_plot
  value_range <- range(matrix_values, finite = TRUE)
  if (!all(is.finite(value_range))) stop("Heatmap range must be finite.", call. = FALSE)
  if (diff(value_range) == 0) value_range <- value_range + c(-0.5, 0.5)
  breaks <- seq(value_range[1], value_range[2], length.out = 101L)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100)
  sample_valid <- prepared$sidecar$ValidDenominator[
    match(prepared$sample_order, prepared$sidecar$SampleID)
  ]
  annotation <- data.frame(
    Denominator = factor(
      ifelse(sample_valid, "valid", "undefined (0 classified reads)"),
      levels = c("valid", "undefined (0 classified reads)")
    ),
    row.names = prepared$sample_order,
    stringsAsFactors = FALSE
  )
  if (identical(mode, "cohort")) {
    annotation$Group <- factor(
      prepared$sidecar$Group[match(prepared$sample_order,
                                   prepared$sidecar$SampleID)],
      levels = prepared$group_order
    )
    annotation <- annotation[, c("Group", "Denominator"), drop = FALSE]
  }
  with_png_device(path, width = if (identical(mode, "cohort")) {
    max(8, min(24, 5 + 0.28 * length(prepared$sample_order)))
  } else 6, height = if (identical(mode, "cohort")) {
    max(6, min(12, 3.5 + 0.45 * nrow(matrix_values)))
  } else 7, dpi = 300, draw = function() {
    pheatmap::pheatmap(
      matrix_values,
      annotation_col = annotation,
      color = colors,
      breaks = breaks,
      na_col = "#E6E6E6",
      annotation_colors = list(Denominator = c(
        "valid" = "#FFFFFF",
        "undefined (0 classified reads)" = "#BDBDBD"
      )),
      border_color = "#D9D9D9",
      cluster_rows = if (identical(mode, "cohort") && !is.null(prepared$row_hclust)) {
        prepared$row_hclust
      } else FALSE,
      cluster_cols = FALSE,
      angle_col = if (identical(mode, "cohort")) 90 else 0,
      main = sprintf("%s abundance (%s; classified-read denominator; undefined = 0 classified reads)",
                     tools::toTitleCase(rank), prepared$transform),
      fontsize = 9
    )
  })
  invisible(path)
}

composition_stacked_columns <- c(
  "Rank", "TaxonPath", "DisplayTaxon", "SampleID", "Group",
  "RelativeAbundance", "MeanRelativeAbundance", "IsOther",
  "ValidDenominator", "StackOrder", "SampleOrder"
)

composition_group_columns <- c(
  "Rank", "Group", "TaxonPath", "DisplayTaxon", "MeanRelativeAbundance",
  "SamplesTotal", "SamplesUsed", "SamplesExcludedZeroClassified", "IsOther",
  "StackOrder", "GroupOrder"
)

composition_heatmap_columns <- c(
  "Rank", "TaxonPath", "DisplayTaxon", "SampleID", "Group",
  "RelativeAbundance", "Transform", "PseudoCount", "TransformedValue",
  "IsOther", "RowOrder", "ColumnOrder", "ValidDenominator"
)

validate_composition_display_labels <- function(data) {
  paths <- unique(as.character(data$TaxonPath[data$TaxonPath != OTHER_TAXON_PATH]))
  if (length(paths)) {
    map <- make_taxon_display_map(data.frame(TaxonPath = paths, stringsAsFactors = FALSE))
    expected <- map$DisplayTaxon[match(as.character(data$TaxonPath), map$TaxonPath)]
  } else {
    expected <- rep(NA_character_, nrow(data))
  }
  expected[data$TaxonPath == OTHER_TAXON_PATH] <- "Other"
  if (anyNA(expected) || !identical(as.character(data$DisplayTaxon), expected)) {
    stop("Composition DisplayTaxon values are inconsistent with TaxonPath identity.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_sidecar_sample_order <- function(data, expected_samples) {
  sample_order <- vapply(expected_samples, function(sample_id) {
    values <- unique(data$SampleOrder[data$SampleID == sample_id])
    if (length(values) != 1L) NA_real_ else as.numeric(values)
  }, numeric(1))
  if (anyNA(sample_order) || !identical(unname(sort(as.numeric(sample_order))),
                                         as.numeric(seq_along(expected_samples)))) {
    stop("Composition SampleOrder is not a contiguous deterministic sample order.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_sidecar_taxon_sets <- function(data, expected_samples, label) {
  path_sets <- lapply(expected_samples, function(sample_id) {
    sort(unique(as.character(data$TaxonPath[data$SampleID == sample_id])), method = "radix")
  })
  if (length(path_sets) > 1L && any(vapply(path_sets[-1L], function(paths) {
    !identical(paths, path_sets[[1L]])
  }, logical(1)))) {
    stop(sprintf("%s sidecar does not use one complete TaxonPath set per sample.", label),
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_stacked_sidecar <- function(x, expected_samples, expected_rank = NULL) {
  if (!is.data.frame(x) || !identical(names(x), composition_stacked_columns) || !nrow(x)) {
    stop("Invalid stacked composition sidecar schema.", call. = FALSE)
  }
  if (!is.null(expected_rank) && !all(as.character(x$Rank) == expected_rank)) {
    stop("Stacked composition sidecar has an unexpected rank.", call. = FALSE)
  }
  if (!setequal(as.character(unique(x$SampleID)), as.character(expected_samples)) ||
      anyDuplicated(paste(x$SampleID, x$TaxonPath, sep = "\r"))) {
    stop("Stacked composition sidecar has an invalid sample set or key.", call. = FALSE)
  }
  if (any(!is.finite(x$RelativeAbundance)) || any(x$RelativeAbundance < 0) ||
      any(x$RelativeAbundance > 1) || any(!is.finite(x$MeanRelativeAbundance)) ||
      any(x$MeanRelativeAbundance < 0) || any(x$MeanRelativeAbundance > 1) ||
      anyNA(x$ValidDenominator) || anyNA(x$IsOther)) {
    stop("Stacked composition sidecar contains invalid abundance or flag values.", call. = FALSE)
  }
  validate_composition_display_labels(x)
  validate_sidecar_sample_order(x, expected_samples)
  validate_sidecar_taxon_sets(x, expected_samples, "Stacked composition")
  for (sample_id in expected_samples) {
    sample <- x[x$SampleID == sample_id, , drop = FALSE]
    if (length(unique(sample$IsOther[sample$TaxonPath == OTHER_TAXON_PATH])) > 1L ||
        sum(sample$IsOther) > 1L ||
        any(sample$IsOther != (sample$TaxonPath == OTHER_TAXON_PATH))) {
      stop("Stacked composition sidecar has an invalid Other row.", call. = FALSE)
    }
    if (isTRUE(sample$ValidDenominator[1])) {
      if (abs(sum(sample$RelativeAbundance) - 1) > 1e-8) {
        stop("Valid stacked composition does not conserve to one.", call. = FALSE)
      }
    } else if (abs(sum(sample$RelativeAbundance)) > 1e-12) {
      stop("Invalid stacked composition must have zero abundance.", call. = FALSE)
    }
    stack_order <- sort(unique(as.integer(sample$StackOrder)))
    if (!identical(stack_order, seq_len(nrow(sample)))) {
      stop("Stacked composition StackOrder is not contiguous.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

validate_group_sidecar <- function(group, sample, expected_rank = NULL) {
  if (!is.data.frame(group) || !identical(names(group), composition_group_columns) || !nrow(group) ||
      !is.data.frame(sample)) {
    stop("Invalid group composition sidecar schema.", call. = FALSE)
  }
  if (!is.null(expected_rank) && !all(as.character(group$Rank) == expected_rank)) {
    stop("Group composition sidecar has an unexpected rank.", call. = FALSE)
  }
  if (anyDuplicated(paste(group$Group, group$TaxonPath, sep = "\r"))) {
    stop("Group composition sidecar has duplicate group/taxon keys.", call. = FALSE)
  }
  validate_composition_display_labels(group)
  sample_info <- unique(sample[, c("SampleID", "Group", "ValidDenominator"), drop = FALSE])
  if (anyDuplicated(sample_info$SampleID) || anyNA(sample_info$ValidDenominator)) {
    stop("Sample composition data cannot establish group denominator status.", call. = FALSE)
  }
  for (g in unique(as.character(group$Group))) {
    g_samples <- sample_info$SampleID[sample_info$Group == g]
    g_valid <- sample_info$ValidDenominator[sample_info$Group == g]
    g_rows <- group[group$Group == g, , drop = FALSE]
    if (any(g_rows$SamplesTotal != length(g_samples)) ||
        any(g_rows$SamplesUsed != sum(g_valid)) ||
        any(g_rows$SamplesExcludedZeroClassified != sum(!g_valid))) {
      stop("Group composition sample counts do not match sample denominator status.", call. = FALSE)
    }
    if (sum(g_valid) == 0L) {
      if (any(!is.na(g_rows$MeanRelativeAbundance))) {
        stop("A zero-valid group must have NA means.", call. = FALSE)
      }
    } else {
      if (any(!is.finite(g_rows$MeanRelativeAbundance)) ||
          any(g_rows$MeanRelativeAbundance < 0) || any(g_rows$MeanRelativeAbundance > 1)) {
        stop("A valid group has invalid mean abundance values.", call. = FALSE)
      }
      if (abs(sum(g_rows$MeanRelativeAbundance) - 1) > 1e-8) {
        stop("Valid group means do not conserve to one.", call. = FALSE)
      }
    }
    for (i in seq_len(nrow(g_rows))) {
      path <- as.character(g_rows$TaxonPath[i])
      values <- sample$RelativeAbundance[
        sample$Group == g & sample$TaxonPath == path & sample$ValidDenominator
      ]
      if (!length(values)) {
        if (!is.na(g_rows$MeanRelativeAbundance[i])) {
          stop("Group mean is defined despite zero valid sample values.", call. = FALSE)
        }
      } else if (!isTRUE(all.equal(as.numeric(g_rows$MeanRelativeAbundance[i]),
                                   mean(values), tolerance = 1e-12))) {
        stop("Group mean is not the arithmetic mean of valid sample proportions.", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

validate_heatmap_sidecar <- function(x, expected_samples, cfg, expected_rank = NULL) {
  if (!is.data.frame(x) || !identical(names(x), composition_heatmap_columns) || !nrow(x)) {
    stop("Invalid heatmap composition sidecar schema.", call. = FALSE)
  }
  if (!is.null(expected_rank) && !all(as.character(x$Rank) == expected_rank)) {
    stop("Heatmap composition sidecar has an unexpected rank.", call. = FALSE)
  }
  if (!setequal(as.character(unique(x$SampleID)), as.character(expected_samples)) ||
      anyDuplicated(paste(x$SampleID, x$TaxonPath, sep = "\r"))) {
    stop("Heatmap composition sidecar has an invalid sample set or key.", call. = FALSE)
  }
  if (any(!is.finite(x$RelativeAbundance)) || any(x$RelativeAbundance < 0) ||
      any(x$RelativeAbundance > 1) || any(!is.finite(x$TransformedValue)) ||
      anyNA(x$ValidDenominator) || anyNA(x$IsOther)) {
    stop("Heatmap composition sidecar contains invalid abundance or flag values.", call. = FALSE)
  }
  if (any(x$IsOther != (x$TaxonPath == OTHER_TAXON_PATH))) {
    stop("Heatmap composition sidecar has an invalid Other row.", call. = FALSE)
  }
  validate_composition_display_labels(x)
  validate_sidecar_taxon_sets(x, expected_samples, "Heatmap composition")
  transformed <- as.character(unique(x$Transform))
  if (length(transformed) != 1L || !transformed %in% c("none", "log10_relative")) {
    stop("Heatmap composition sidecar has an unsupported transform.", call. = FALSE)
  }
  include_other <- isTRUE(cfg$composition$heatmap_include_other)
  top_n <- as.integer(cfg$composition$heatmap_top_n_taxa)
  column_orders <- vapply(expected_samples, function(sample_id) {
    values <- unique(x$ColumnOrder[x$SampleID == sample_id])
    if (length(values) != 1L) NA_integer_ else as.integer(values)
  }, integer(1))
  if (anyNA(column_orders) || !identical(unname(sort(column_orders)), seq_along(expected_samples))) {
    stop("Heatmap ColumnOrder is not a contiguous deterministic sample order.", call. = FALSE)
  }
  for (sample_id in expected_samples) {
    sample <- x[x$SampleID == sample_id, , drop = FALSE]
    if (sum(sample$IsOther) > 1L || any(sample$IsOther != (sample$TaxonPath == OTHER_TAXON_PATH))) {
      stop("Heatmap composition sidecar has duplicate or inconsistent Other rows.", call. = FALSE)
    }
    if (isTRUE(sample$ValidDenominator[1])) {
      total <- sum(sample$RelativeAbundance)
      if (include_other && abs(total - 1) > 1e-8) {
        stop("Valid heatmap composition does not conserve to one with Other enabled.", call. = FALSE)
      }
      if (!include_other && (total < -1e-12 || total > 1 + 1e-8)) {
        stop("Valid heatmap composition exceeds one with Other disabled.", call. = FALSE)
      }
    } else if (abs(sum(sample$RelativeAbundance)) > 1e-12 ||
               any(sample$IsOther & sample$RelativeAbundance > 0)) {
      stop("Invalid heatmap composition must have zero abundance and no Other.", call. = FALSE)
    }
    if (sum(!sample$IsOther) > top_n) {
      stop("Heatmap composition exceeds the configured named taxon limit.", call. = FALSE)
    }
    if (!identical(sort(unique(as.integer(sample$RowOrder))), seq_len(nrow(sample)))) {
      stop("Heatmap row or column order is not contiguous.", call. = FALSE)
    }
  }
  if (identical(transformed, "none")) {
    if (any(!is.na(x$PseudoCount)) ||
        !isTRUE(all.equal(as.numeric(x$TransformedValue), as.numeric(x$RelativeAbundance),
                          tolerance = 0, check.attributes = FALSE))) {
      stop("A none heatmap transform must report no pseudocount and unchanged values.", call. = FALSE)
    }
  } else {
    if (any(!is.finite(x$PseudoCount)) || any(x$PseudoCount <= 0) ||
        length(unique(round(x$PseudoCount, 15))) != 1L) {
      stop("A log10 heatmap transform must report one positive pseudocount.", call. = FALSE)
    }
    positive <- x$RelativeAbundance[x$ValidDenominator & x$RelativeAbundance > 0]
    if (!length(positive) || abs(unique(x$PseudoCount)[1] - min(positive) / 2) > 1e-12) {
      stop("Heatmap pseudocount was not derived from valid positive abundance.", call. = FALSE)
    }
    expected_transformed <- log10(x$RelativeAbundance + unique(x$PseudoCount)[1])
    if (!isTRUE(all.equal(as.numeric(x$TransformedValue), expected_transformed,
                          tolerance = 1e-12, check.attributes = FALSE))) {
      stop("Heatmap transformed values do not match the log10 pseudocount contract.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

write_composition_table <- function(data, path, columns) {
  if (!is.data.frame(data) || !identical(names(data), columns)) {
    stop(sprintf("Composition sidecar schema mismatch for '%s'.", path), call. = FALSE)
  }
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE)
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    stop(sprintf("Composition sidecar was not written: '%s'.", path), call. = FALSE)
  }
  invisible(path)
}

run_taxa_composition <- function(context) {
  cfg <- context$config
  comp_dir <- cfg$output$dirs$composition
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  all_outputs <- character(0)

  top_n <- cfg$composition$top_n_taxa %||% 15L
  stacked_ranks <- cfg$composition$stacked_bar_ranks %||% c("phylum", "family", "genus")
  heatmap_ranks <- cfg$composition$heatmap_ranks %||% {
    legacy <- cfg$composition[["heatmap_rank"]]
    if (is.null(legacy)) c("phylum", "family", "genus") else legacy
  }
  heatmap_transform <- cfg$composition$heatmap_transform %||% "log10_relative"

  samples <- context$samples
  unclass_idx <- context$unclass_index
  count_matrix <- context$count_matrix
  taxonomy_df <- context$taxonomy
  sample_totals <- colSums(count_matrix)
  unclass_counts <- stats::setNames(context$sample_stats$UnclassifiedReads,
                                    context$sample_stats$SampleID)
  class_totals <- sample_totals - unclass_counts
  valid_samples <- samples[class_totals[samples] > 0]

  classification_df <- data.frame(
    SampleID = samples,
    TotalReads = as.numeric(sample_totals[samples]),
    ClassifiedReads = as.numeric(class_totals[samples]),
    UnclassifiedReads = as.numeric(unclass_counts[samples]),
    ClassifiedFraction = ifelse(sample_totals[samples] > 0,
                                as.numeric(class_totals[samples] / sample_totals[samples]), 0),
    stringsAsFactors = FALSE
  )
  classification_file <- file.path(comp_dir, "classification_fraction.tsv")
  write.table(classification_df, classification_file, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, classification_file)

  rank_tables <- list()
  class_counts <- count_matrix[-unclass_idx, , drop = FALSE]
  class_tax <- taxonomy_df[-unclass_idx, , drop = FALSE]
  for (rk in RANKS_TO_ANALYZE) {
    rk_idx <- which(colnames(taxonomy_df) == rk)
    rk_prefixes <- apply(class_tax[, 1:rk_idx, drop = FALSE], 1, paste, collapse = ";")
    leaf_names <- class_tax[[rk]]
    contextualized_names <- ifelse(
      leaf_names %in% c("Unknown", "unclassified", "uncultured"),
      sprintf("%s (%s)", leaf_names, class_tax[[rk_idx - 1]]), leaf_names
    )
    agg_df <- data.frame(
      TaxonPath = rk_prefixes, Taxon = contextualized_names, class_counts,
      check.names = FALSE, stringsAsFactors = FALSE
    ) %>%
      group_by(TaxonPath, Taxon) %>%
      summarise(across(all_of(samples), sum), .groups = "drop")
    duplicated_label <- duplicated(agg_df$Taxon) | duplicated(agg_df$Taxon, fromLast = TRUE)
    agg_df$Taxon[duplicated_label] <- sprintf(
      "%s [%s]", agg_df$Taxon[duplicated_label], agg_df$TaxonPath[duplicated_label]
    )
    rel_df <- agg_df
    for (s in samples) {
      denom <- class_totals[s]
      rel_df[[s]] <- if (denom > 0) rel_df[[s]] / denom else 0
    }
    attr(rel_df, "rank") <- rk
    count_file <- file.path(comp_dir, sprintf("count_%s.tsv", rk))
    rel_file <- file.path(comp_dir, sprintf("rel_abundance_%s.tsv", rk))
    write.table(agg_df, count_file, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(rel_df, rel_file, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, count_file, rel_file)
    rank_tables[[rk]] <- list(counts = agg_df, rel = rel_df)
  }

  # Preserve the historical horizontal single-sample plots and their names.
  if (length(samples) == 1L) {
    s_col <- samples[1]
    phylum_rel <- rank_tables[["phylum"]]$rel %>%
      select(Taxon, all_of(s_col)) %>% rename(rel = all_of(s_col)) %>% arrange(desc(rel))
    n_keep_phylum <- min(7L, nrow(phylum_rel))
    phylum_top <- phylum_rel %>%
      mutate(DisplayTaxon = if_else(row_number() <= n_keep_phylum, Taxon, "Other phyla")) %>%
      group_by(DisplayTaxon) %>% summarise(rel = sum(rel), .groups = "drop") %>% arrange(desc(rel))
    phylum_top$DisplayTaxon <- factor(phylum_top$DisplayTaxon, levels = rev(phylum_top$DisplayTaxon))
    p_phylum <- ggplot(phylum_top, aes(x = DisplayTaxon, y = rel, fill = DisplayTaxon)) +
      geom_col(width = 0.7) + geom_text(aes(label = sprintf("%.1f%%", 100 * rel)), hjust = -0.15, size = 3.3) +
      scale_fill_manual(values = get_phylum_colors(levels(phylum_top$DisplayTaxon))) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.18))) + coord_flip() +
      labs(title = sprintf("Phylum-Level Composition (%s)", s_col),
           subtitle = sprintf("Relative to %s classified reads", scales::comma(class_totals[s_col])),
           x = NULL, y = "Relative Abundance") + theme_amplicon() + theme(legend.position = "none")
    p_path <- file.path(comp_dir, "04a_phylum_composition.png")
    save_plot(p_path, p_phylum, width = 7.5, height = 5)
    all_outputs <- c(all_outputs, p_path)

    family_counts <- rank_tables[["family"]]$counts %>% select(Taxon, all_of(s_col)) %>%
      rename(count = all_of(s_col)) %>% arrange(desc(count)) %>% slice_head(n = top_n)
    family_counts$Taxon <- factor(family_counts$Taxon, levels = rev(family_counts$Taxon))
    p_family <- ggplot(family_counts, aes(x = Taxon, y = count)) + geom_col(width = 0.7, fill = "#1b9e77") +
      geom_text(aes(label = scales::comma(count)), hjust = -0.15, size = 3) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) + coord_flip() +
      labs(title = sprintf("Top %d Families by Read Count (%s)", top_n, s_col), x = NULL, y = "Read count") +
      theme_amplicon()
    p_path <- file.path(comp_dir, "04b_family_composition.png")
    save_plot(p_path, p_family, width = 8.5, height = 6)
    all_outputs <- c(all_outputs, p_path)

    genus_counts <- rank_tables[["genus"]]$counts %>% select(Taxon, all_of(s_col)) %>%
      rename(count = all_of(s_col)) %>% arrange(desc(count)) %>% slice_head(n = top_n)
    genus_counts$Taxon <- factor(genus_counts$Taxon, levels = rev(genus_counts$Taxon))
    p_genus <- ggplot(genus_counts, aes(x = Taxon, y = count)) + geom_col(width = 0.7, fill = "#d95f02") +
      geom_text(aes(label = scales::comma(count)), hjust = -0.15, size = 3) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) + coord_flip() +
      labs(title = sprintf("Top %d Genera by Read Count (%s)", top_n, s_col), x = NULL, y = "Read count") +
      theme_amplicon()
    p_path <- file.path(comp_dir, "04c_genus_composition.png")
    save_plot(p_path, p_genus, width = 8.5, height = 6)
    all_outputs <- c(all_outputs, p_path)

    species_rel <- rank_tables[["species"]]$rel %>% select(Taxon, all_of(s_col)) %>%
      rename(rel = all_of(s_col)) %>% arrange(desc(rel)) %>% slice_head(n = top_n)
    species_rel$Taxon <- factor(species_rel$Taxon, levels = rev(species_rel$Taxon))
    p_species <- ggplot(species_rel, aes(x = Taxon, y = rel)) + geom_col(width = 0.7, fill = "#7570b3") +
      geom_text(aes(label = sprintf("%.2f%%", 100 * rel)), hjust = -0.1, size = 3) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.25))) + coord_flip() +
      labs(title = sprintf("Top %d Species by Relative Abundance (%s)", top_n, s_col), x = NULL,
           y = "Relative Abundance (Classified reads)") + theme_amplicon() +
      theme(axis.text.y = element_text(face = "italic"))
    p_path <- file.path(comp_dir, "04d_species_composition.png")
    save_plot(p_path, p_species, width = 9, height = 6)
    all_outputs <- c(all_outputs, p_path)
  }

  plot_meta <- plot_metadata(samples, context$metadata, valid_samples = valid_samples)
  stacked_sidecar_columns <- composition_stacked_columns
  group_sidecar_columns <- composition_group_columns

  for (rk in stacked_ranks) {
    selection <- select_display_taxa(
      rank_tables[[rk]]$rel, samples, valid_samples,
      min_mean = cfg$composition$stacked_bar_min_mean_relative,
      min_n = cfg$composition$stacked_bar_min_taxa,
      max_n = cfg$composition$stacked_bar_max_taxa
    )
    display_long <- if (nrow(selection)) {
      collapse_rank_for_display(rank_tables[[rk]]$rel, selection$TaxonPath, samples, plot_meta)
    } else NULL
    positive_display <- !is.null(display_long) && any(display_long$RelativeAbundance > 0 &
      display_long$ValidDenominator)
    if (is.null(display_long) || !nrow(display_long) || !positive_display) {
      skipped <- data.frame(Rank = rk,
                            Reason = "No valid positive classified relative abundance exists at this rank.",
                            stringsAsFactors = FALSE)
      skip_path <- file.path(comp_dir, sprintf("04_%s_stacked_skipped.tsv", rk))
      write_composition_table(skipped, skip_path, c("Rank", "Reason"))
      all_outputs <- c(all_outputs, skip_path)
      if (length(samples) > 1L) {
        group_skip_path <- file.path(comp_dir, sprintf("04_%s_group_mean_stacked_skipped.tsv", rk))
        write_composition_table(skipped, group_skip_path, c("Rank", "Reason"))
        all_outputs <- c(all_outputs, group_skip_path)
      }
      next
    }

    display_long$Rank <- rk
    display_long <- display_long[, c("Rank", setdiff(names(display_long), "Rank")), drop = FALSE]
    display_long_full <- display_long
    display_long <- display_long[, stacked_sidecar_columns, drop = FALSE]
    validate_stacked_sidecar(display_long, samples, expected_rank = rk)
    stacked_path <- file.path(comp_dir, sprintf("04_%s_stacked.tsv", rk))
    write_composition_table(display_long, stacked_path, stacked_sidecar_columns)
    all_outputs <- c(all_outputs, stacked_path)

    tax_paths <- unique(display_long$TaxonPath[order(display_long$StackOrder)])
    tax_labels <- display_long$DisplayTaxon[match(tax_paths, display_long$TaxonPath)]
    colors <- composition_colors(tax_paths, tax_labels)
    sample_order <- if (length(samples) > 1L) ordered_cohort_samples(samples, context$metadata)$sample_order else samples
    group_order <- if (length(samples) > 1L) ordered_cohort_samples(samples, context$metadata)$group_order else NULL
    plot <- build_stacked_taxa_plot(display_long, colors, sample_order, group_order,
                                    single = length(samples) == 1L, rank = rk)
    stacked_png <- file.path(comp_dir, sprintf("04_%s_stacked.png", rk))
    save_plot(stacked_png, plot, width = if (length(samples) == 1L) 7 else {
      max(10, min(24, 6 + 0.32 * length(samples)))
    }, height = if (length(samples) == 1L) 7 else 8, dpi = 300)
    all_outputs <- c(all_outputs, stacked_png)

    if (length(samples) > 1L) {
      group_means <- summarize_group_means(display_long_full)
      group_means$Rank <- rk
      group_means <- group_means[, c("Rank", setdiff(names(group_means), "Rank")), drop = FALSE]
      group_means <- group_means[, group_sidecar_columns, drop = FALSE]
      validate_group_sidecar(group_means, display_long_full, expected_rank = rk)
      group_path <- file.path(comp_dir, sprintf("04_%s_group_mean_stacked.tsv", rk))
      write_composition_table(group_means, group_path, group_sidecar_columns)
      all_outputs <- c(all_outputs, group_path)
      group_plot_data <- group_means
      group_plot_data$SampleID <- group_plot_data$Group
      group_plot_data$RelativeAbundance <- group_plot_data$MeanRelativeAbundance
      group_plot_data$ValidDenominator <- group_plot_data$SamplesUsed > 0L
      group_plot_data$SampleOrder <- group_plot_data$GroupOrder
      group_plot <- build_stacked_taxa_plot(
        group_plot_data[, c("TaxonPath", "DisplayTaxon", "SampleID", "Group",
                            "RelativeAbundance", "IsOther", "ValidDenominator",
                            "StackOrder", "GroupOrder")],
        colors, group_order, group_order, single = TRUE, rank = paste(rk, "group mean"),
        subtitle = paste(
          "Arithmetic mean of per-sample classified-read relative abundance;",
          "samples are not pooled by read depth."
        ),
        invalid_groups = unique(group_means$Group[group_means$SamplesUsed == 0L])
      )
      group_png <- file.path(comp_dir, sprintf("04_%s_group_mean_stacked.png", rk))
      save_plot(group_png, group_plot,
                width = max(7, min(18, 3 + 1.25 * length(group_order))), height = 7, dpi = 300)
      all_outputs <- c(all_outputs, group_png)
    }
  }

  for (rk in heatmap_ranks) {
    prepared <- prepare_heatmap_data(
      rank_tables[[rk]]$rel, samples, plot_meta,
      top_n = cfg$composition$heatmap_top_n_taxa,
      include_other = cfg$composition$heatmap_include_other,
      transform = heatmap_transform, rank = rk
    )
    if (isTRUE(prepared$skipped)) {
      skip_path <- file.path(comp_dir, sprintf("04_heatmap_%s_skipped.tsv", rk))
      write_composition_table(prepared$skip, skip_path, c("Rank", "Reason"))
      all_outputs <- c(all_outputs, skip_path)
      next
    }
    validate_heatmap_sidecar(
      prepared$sidecar, samples, cfg, expected_rank = rk
    )
    heatmap_tsv <- file.path(comp_dir, sprintf("04_heatmap_%s.tsv", rk))
    write_composition_table(prepared$sidecar, heatmap_tsv,
                            c("Rank", "TaxonPath", "DisplayTaxon", "SampleID", "Group",
                              "RelativeAbundance", "Transform", "PseudoCount",
                              "TransformedValue", "IsOther", "RowOrder", "ColumnOrder",
                              "ValidDenominator"))
    all_outputs <- c(all_outputs, heatmap_tsv)
    heatmap_png <- file.path(comp_dir, sprintf("04_heatmap_%s.png", rk))
    draw_taxa_heatmap(prepared, heatmap_png,
                      mode = if (length(samples) > 1L) "cohort" else "single", rank = rk)
    all_outputs <- c(all_outputs, heatmap_png)
  }

  classification_long <- classification_df %>%
    select(SampleID, ClassifiedReads, UnclassifiedReads) %>%
    pivot_longer(-SampleID, names_to = "Status", values_to = "Reads") %>%
    mutate(Status = recode(Status, ClassifiedReads = "Classified", UnclassifiedReads = "Unclassified"))
  p_class <- ggplot(classification_long, aes(x = SampleID, y = Reads, fill = Status)) +
    geom_col(position = "fill") + scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(Classified = "#1b9e77", Unclassified = "#bdbdbd")) +
    labs(title = "All-read Classification Fraction", x = NULL, y = "Fraction of all reads") +
    theme_amplicon() + theme(axis.text.x = element_text(angle = if (length(samples) > 5) 45 else 0, hjust = 1))
  classification_plot <- file.path(comp_dir, "04_classification_fraction.png")
  save_plot(classification_plot, p_class, width = max(7, min(14, 0.45 * length(samples) + 5)), height = 5)
  all_outputs <- c(all_outputs, classification_plot)

  list(status = "completed", outputs = all_outputs)
}
