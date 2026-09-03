# =============================================================================
# Module 01: QC Diagnostics (Read Length, Classification Reconciliation)
# =============================================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

run_qc <- function(context) {
  cfg <- context$config
  qc_dir <- cfg$output$dirs$qc
  assignments_map <- context$assignments
  
  if (is.null(assignments_map) || length(assignments_map) == 0) {
    return(list(
      status = "skipped",
      reason = "No assignments mapping configured in input.assignments",
      outputs = character(0)
    ))
  }
  
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_outputs <- character(0)
  reconciliation_rows <- list()
  length_summary_rows <- list()
  
  # Target lengths from params.json or config fallbacks
  min_len_target <- if (!is.null(context$params$min_len)) {
    context$params$min_len
  } else {
    cfg$qc$target_min_length
  }
  max_len_target <- if (!is.null(context$params$max_len)) {
    context$params$max_len
  } else {
    cfg$qc$target_max_length
  }
  
  display_min <- cfg$qc$display_min_length %||% 1200L
  display_max <- cfg$qc$display_max_length %||% 1800L
  
  for (sample_id in context$samples) {
    asgn_path <- assignments_map[[sample_id]]
    if (is.null(asgn_path) || !file.exists(asgn_path)) {
      next
    }
    
    sample_out_dir <- file.path(qc_dir, sanitize_filename(sample_id))
    dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Expected counts from abundance context
    stat_row <- context$sample_stats[context$sample_stats$SampleID == sample_id, ]
    exp_total <- stat_row$TotalReads[1]
    exp_class <- stat_row$ClassifiedReads[1]
    exp_unclass <- stat_row$UnclassifiedReads[1]
    
    reads <- read_assignments_file(
      path = asgn_path,
      sample_id = sample_id,
      expected_total = exp_total,
      expected_classified = exp_class,
      expected_unclassified = exp_unclass
    )
    
    n_status_C <- sum(reads$status == "C")
    n_status_U <- sum(reads$status == "U")
    n_qc_reclass <- sum(reads$status == "C" & reads$taxid == 0)
    n_eff_class <- sum(reads$effective_classified)
    n_eff_unclass <- sum(!reads$effective_classified)
    
    # 1. Reconciliation Table Row
    reconciliation_rows[[sample_id]] <- data.frame(
      SampleID = sample_id,
      TotalReads = exp_total,
      AbundanceClassified = exp_class,
      AbundanceUnclassified = exp_unclass,
      AssignmentTotal = nrow(reads),
      AssignmentRawStatusC = n_status_C,
      AssignmentRawStatusU = n_status_U,
      AssignmentFailedFilter = n_qc_reclass,
      EffectiveClassified = n_eff_class,
      EffectiveUnclassified = n_eff_unclass,
      ReconciliationPass = (n_eff_class == exp_class) && (nrow(reads) == exp_total),
      stringsAsFactors = FALSE
    )
    
    # 2. Length summary statistics
    reads$status_category <- ifelse(
      reads$status == "C" & reads$taxid != 0, "Classified",
      ifelse(reads$status == "C" & reads$taxid == 0, "QC-filtered", "Never aligned")
    )
    reads$effective_status <- ifelse(reads$effective_classified, "Classified", "Unclassified")
    
    for (cat_name in c("All", "Classified", "Unclassified", "QC-filtered")) {
      sub_lens <- if (cat_name == "All") {
        reads$read_length
      } else if (cat_name %in% c("Classified", "Unclassified")) {
        reads$read_length[reads$effective_status == cat_name]
      } else {
        reads$read_length[reads$status_category == cat_name]
      }
      
      if (length(sub_lens) > 0) {
        length_summary_rows[[paste(sample_id, cat_name, sep = "_")]] <- data.frame(
          SampleID = sample_id,
          Category = cat_name,
          Count = length(sub_lens),
          Min = min(sub_lens),
          Q25 = as.numeric(quantile(sub_lens, 0.25)),
          Median = median(sub_lens),
          Mean = mean(sub_lens),
          Q75 = as.numeric(quantile(sub_lens, 0.75)),
          Max = max(sub_lens),
          SD = sd(sub_lens),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # 3. Figure 1a: Donut Plot
    donut_df <- reads %>%
      count(effective_status) %>%
      mutate(
        frac = n / sum(n),
        ymax = cumsum(frac),
        ymin = c(0, head(ymax, -1)),
        label = sprintf("%s\n%s reads\n(%.1f%%)", effective_status, scales::comma(n), 100 * frac)
      )
    
    p1a <- ggplot(donut_df, aes(ymin = ymin, ymax = ymax, xmin = 3, xmax = 4, fill = effective_status)) +
      geom_rect(color = "white", linewidth = 1.2) +
      coord_polar(theta = "y") +
      xlim(c(1, 4)) +
      scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
      annotate("text", x = 1, y = 0, label = sprintf("%s\nreads", scales::comma(sum(donut_df$n))),
               size = 4.2, fontface = "bold") +
      theme_void() +
      labs(title = sprintf("Read Classification Outcome: %s", sample_id), fill = NULL) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) +
      geom_text(aes(x = 3.5, y = (ymin + ymax) / 2, label = label), inherit.aes = FALSE,
                data = donut_df, color = "black", size = 3.3, fontface = "bold")
    
    p1a_path <- file.path(sample_out_dir, "01a_classification_donut.png")
    save_plot(p1a_path, p1a, width = 5.5, height = 5.5)
    all_outputs <- c(all_outputs, p1a_path)
    
    # 4. Figure 1b: Read Length Histogram
    p1b <- ggplot(reads, aes(x = read_length, fill = effective_status)) +
      geom_histogram(binwidth = 10, alpha = 0.85, position = "identity") +
      scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
      coord_cartesian(xlim = c(display_min, display_max)) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = sprintf("Read Length Distribution: %s", sample_id),
        subtitle = if (!is.null(min_len_target) && !is.null(max_len_target)) {
          sprintf("Expected full-length 16S ~1400-1550 bp | Upstream target [%d, %d] bp",
                  min_len_target, max_len_target)
        } else {
          "Expected full-length 16S ~1400-1550 bp"
        },
        x = "Read length (bp)", y = "Read count", fill = "Effective Status"
      ) +
      theme_amplicon() +
      theme(legend.position = "bottom")
    
    if (!is.null(min_len_target)) {
      p1b <- p1b + geom_vline(xintercept = min_len_target, linetype = "dotted", color = "grey30")
    }
    if (!is.null(max_len_target)) {
      p1b <- p1b + geom_vline(xintercept = max_len_target, linetype = "dotted", color = "grey30")
    }
    
    p1b_path <- file.path(sample_out_dir, "01b_read_length_distribution.png")
    save_plot(p1b_path, p1b, width = 7.5, height = 5)
    all_outputs <- c(all_outputs, p1b_path)
    
    # 5. Figure 1c: Diagnostic Bar (Why raw status != effective classification)
    qc_diag_df <- data.frame(
      Category = factor(
        c("True classified\n(aligned + passed QC)",
          "QC-filtered\n(aligned but failed\ncoverage/identity)",
          "Never aligned\n(status = U)"),
        levels = c("True classified\n(aligned + passed QC)",
                   "QC-filtered\n(aligned but failed\ncoverage/identity)",
                   "Never aligned\n(status = U)")
      ),
      Count = c(n_eff_class, n_qc_reclass, n_status_U)
    )
    
    p1c <- ggplot(qc_diag_df, aes(x = Category, y = Count, fill = Category)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = scales::comma(Count)), vjust = -0.4, size = 3.6, fontface = "bold") +
      scale_fill_manual(values = c(
        "True classified\n(aligned + passed QC)" = "#1b9e77",
        "QC-filtered\n(aligned but failed\ncoverage/identity)" = "#e6ab02",
        "Never aligned\n(status = U)" = "#bdbdbd"
      )) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = sprintf("Raw Alignment Status vs. Effective Classification (%s)", sample_id),
        subtitle = sprintf(
          "Alignment QC filter reclassifies %s reads to taxid=0 without altering raw status letter",
          scales::comma(n_qc_reclass)
        ),
        x = NULL, y = "Read count"
      ) +
      theme_amplicon() +
      theme(legend.position = "none")
    
    p1c_path <- file.path(sample_out_dir, "01c_qc_filter_diagnostic.png")
    save_plot(p1c_path, p1c, width = 7, height = 5.5)
    all_outputs <- c(all_outputs, p1c_path)
    
    # 6. Figure 1d: Violin & Boxplot of Read Length by Effective Class
    p1d <- ggplot(reads, aes(x = effective_status, y = read_length, fill = effective_status)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
      coord_cartesian(ylim = c(display_min, display_max)) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = sprintf("Read Length Distribution by Classification (%s)", sample_id),
        subtitle = "Full-length 16S amplicon read length comparison",
        x = "Effective Status", y = "Read length (bp)"
      ) +
      theme_amplicon() +
      theme(legend.position = "none")
    
    p1d_path <- file.path(sample_out_dir, "01d_read_length_violin.png")
    save_plot(p1d_path, p1d, width = 6, height = 5)
    all_outputs <- c(all_outputs, p1d_path)
  }
  
  # Export summary TSVs
  reconciliation_file <- file.path(qc_dir, "classification_reconciliation.tsv")
  if (length(reconciliation_rows) > 0) {
    write.table(do.call(rbind, reconciliation_rows), reconciliation_file,
                sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, reconciliation_file)
  }
  
  length_file <- file.path(qc_dir, "read_length_by_status.tsv")
  if (length(length_summary_rows) > 0) {
    write.table(do.call(rbind, length_summary_rows), length_file,
                sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, length_file)
  }
  
  list(
    status = "completed",
    outputs = all_outputs
  )
}
