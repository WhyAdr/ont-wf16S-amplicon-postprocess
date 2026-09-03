# =============================================================================
# Shared Data Layer & I/O Validation
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(digest)
})

RANKS_8 <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")

compute_file_hash <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

sanitize_filename <- function(s) {
  gsub("[^A-Za-z0-9_.-]", "_", s)
}

validate_sample_ids <- function(sample_ids) {
  if (length(sample_ids) == 0) {
    stop("Abundance table validation error: No sample columns detected.", call. = FALSE)
  }
  
  if (any(is.na(sample_ids)) || any(sample_ids == "")) {
    stop("Sample ID validation error: Empty or NA sample ID detected.", call. = FALSE)
  }
  
  if (any(sample_ids %in% c(".", ".."))) {
    stop("Sample ID validation error: Sample ID cannot be '.' or '..'.", call. = FALSE)
  }
  
  if (any(grepl("[/\\\\]", sample_ids))) {
    stop("Sample ID validation error: Sample ID cannot contain path separators ('/' or '\\').", call. = FALSE)
  }
  
  if (any(grepl("[[:cntrl:]]", sample_ids))) {
    stop("Sample ID validation error: Sample ID cannot contain control characters.", call. = FALSE)
  }
  
  sanitized <- vapply(sample_ids, sanitize_filename, character(1))
  if (any(duplicated(sanitized))) {
    stop("Sample ID validation error: Sample IDs collide after filename sanitization.", call. = FALSE)
  }
  
  invisible(TRUE)
}

read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("total"), include_samples = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("Abundance table not found: '%s'", path), call. = FALSE)
  }
  
  # Read header first
  raw_lines <- readLines(path, n = 5)
  if (length(raw_lines) == 0) {
    stop(sprintf("Abundance table is empty: '%s'", path), call. = FALSE)
  }
  
  header <- strsplit(raw_lines[1], "\t")[[1]]
  if (any(duplicated(header))) {
    stop(sprintf("Abundance table has duplicate column names: %s", 
                 paste(header[duplicated(header)], collapse = ", ")), call. = FALSE)
  }
  
  if (!tax_col %in% header) {
    stop(sprintf("Abundance table missing configured tax column '%s'. Columns found: %s",
                 tax_col, paste(header, collapse = ", ")), call. = FALSE)
  }
  
  raw_df <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # Check rows
  if (nrow(raw_df) == 0) {
    stop("Abundance table has 0 data rows.", call. = FALSE)
  }
  
  all_cols <- colnames(raw_df)
  all_sample_cols <- setdiff(all_cols, c(tax_col, aggregate_cols))
  
  validate_sample_ids(all_sample_cols)
  
  # Validate counts for all sample columns
  for (sc in all_sample_cols) {
    vals <- raw_df[[sc]]
    if (!is.numeric(vals)) {
      # Try converting to numeric
      num_vals <- suppressWarnings(as.numeric(vals))
      if (any(is.na(num_vals))) {
        stop(sprintf("Column '%s' contains non-numeric values.", sc), call. = FALSE)
      }
      vals <- num_vals
      raw_df[[sc]] <- vals
    }
    if (any(!is.finite(vals))) {
      stop(sprintf("Column '%s' contains non-finite values (NA, NaN, Inf).", sc), call. = FALSE)
    }
    if (any(vals < 0)) {
      stop(sprintf("Column '%s' contains negative counts.", sc), call. = FALSE)
    }
    if (any(abs(vals - round(vals)) > 1e-6)) {
      stop(sprintf("Column '%s' contains non-integer count values.", sc), call. = FALSE)
    }
  }
  
  # Validate aggregate columns (e.g. 'total') if present
  for (ac in aggregate_cols) {
    if (ac %in% all_cols) {
      actual_sum <- rowSums(as.matrix(raw_df[, all_sample_cols, drop = FALSE]))
      stated_total <- as.numeric(raw_df[[ac]])
      if (any(abs(actual_sum - stated_total) > 1e-4)) {
        diff_idx <- which(abs(actual_sum - stated_total) > 1e-4)[1]
        stop(sprintf(
          "Abundance table aggregate column '%s' does not equal sample row sums at row %d: stated %g != actual %g",
          ac, diff_idx + 1, stated_total[diff_idx], actual_sum[diff_idx]
        ), call. = FALSE)
      }
    }
  }
  
  # Select requested samples
  selected_samples <- if (!is.null(include_samples) && length(include_samples) > 0) {
    missing_sel <- setdiff(include_samples, all_sample_cols)
    if (length(missing_sel) > 0) {
      stop(sprintf("Requested sample(s) not found in abundance table: %s",
                   paste(missing_sel, collapse = ", ")), call. = FALSE)
    }
    include_samples
  } else {
    all_sample_cols
  }
  
  # Validate lineages
  lineages <- raw_df[[tax_col]]
  if (any(duplicated(lineages))) {
    dup <- lineages[duplicated(lineages)][1]
    stop(sprintf("Duplicate full lineage detected in abundance table: '%s'", dup), call. = FALSE)
  }
  
  parsed_lineages <- strsplit(lineages, ";")
  field_counts <- vapply(parsed_lineages, length, integer(1))
  if (any(field_counts != 8)) {
    bad_idx <- which(field_counts != 8)[1]
    stop(sprintf(
      "Lineage schema violation at row %d: expected 8 ranks, found %d ('%s')",
      bad_idx + 1, field_counts[bad_idx], lineages[bad_idx]
    ), call. = FALSE)
  }
  
  # Check unclassified rows
  unclass_indices <- which(vapply(parsed_lineages, function(x) x[1] == "Unclassified", logical(1)))
  if (length(unclass_indices) == 0) {
    stop("Abundance table validation error: No recognizable 'Unclassified' row found.", call. = FALSE)
  }
  if (length(unclass_indices) > 1) {
    stop(sprintf("Abundance table validation error: Multiple (%d) 'Unclassified' rows detected.",
                 length(unclass_indices)), call. = FALSE)
  }
  
  # Build count matrix (taxa x samples)
  count_mat <- as.matrix(raw_df[, selected_samples, drop = FALSE])
  rownames(count_mat) <- lineages
  mode(count_mat) <- "numeric"
  
  # Check sample read counts
  for (s in selected_samples) {
    tot_reads <- sum(count_mat[, s])
    if (tot_reads == 0) {
      stop(sprintf("Sample '%s' has 0 total reads in abundance table.", s), call. = FALSE)
    }
    uncl_reads <- sum(count_mat[unclass_indices, s])
    class_reads <- tot_reads - uncl_reads
    if (class_reads == 0) {
      stop(sprintf("Sample '%s' has 0 classified reads.", s), call. = FALSE)
    }
  }
  
  # Parse taxonomy data frame with 8 ranks
  tax_mat <- do.call(rbind, parsed_lineages)
  colnames(tax_mat) <- RANKS_8
  tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
  tax_df$TaxonPath <- lineages
  
  list(
    count_matrix = count_mat,
    taxonomy = tax_df,
    samples = selected_samples,
    unclass_index = unclass_indices[1],
    all_samples = all_sample_cols
  )
}

read_assignments_file <- function(path, sample_id, expected_total = NULL, expected_classified = NULL, expected_unclassified = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("Assignments file for sample '%s' not found: '%s'", sample_id, path), call. = FALSE)
  }
  
  # Schema: 5 tab-separated fields: status, read_id, taxid, len_field, lineage
  raw_reads <- read.delim(
    path,
    header = FALSE,
    sep = "\t",
    col.names = c("status", "read_id", "taxid", "len_field", "lineage"),
    colClasses = c("character", "character", "character", "character", "character"),
    stringsAsFactors = FALSE,
    quote = "",
    fill = FALSE
  )
  
  if (nrow(raw_reads) == 0) {
    stop(sprintf("Assignments file for sample '%s' is empty.", sample_id), call. = FALSE)
  }
  
  # Check unique read IDs
  if (any(duplicated(raw_reads$read_id))) {
    dup_id <- raw_reads$read_id[duplicated(raw_reads$read_id)][1]
    stop(sprintf("Assignments file for sample '%s' contains duplicate read ID: '%s'", sample_id, dup_id), call. = FALSE)
  }
  
  # Validate status
  valid_statuses <- raw_reads$status %in% c("C", "U")
  if (!all(valid_statuses)) {
    bad_idx <- which(!valid_statuses)[1]
    stop(sprintf("Assignments file for sample '%s' has invalid status '%s' at line %d",
                 sample_id, raw_reads$status[bad_idx], bad_idx), call. = FALSE)
  }
  
  # Parse TaxID
  taxid_num <- suppressWarnings(as.integer(raw_reads$taxid))
  if (any(is.na(taxid_num))) {
    bad_idx <- which(is.na(taxid_num))[1]
    stop(sprintf("Assignments file for sample '%s' has non-integer TaxID '%s' at line %d",
                 sample_id, raw_reads$taxid[bad_idx], bad_idx), call. = FALSE)
  }
  raw_reads$taxid <- taxid_num
  
  # Parse read length defensively: single integer or last numeric part of pipe-delimited string
  # Examples: "0|1481" -> 1481; "1500" -> 1500
  parsed_len <- vapply(raw_reads$len_field, function(lf) {
    if (is.null(lf) || is.na(lf) || lf == "") return(NA_integer_)
    parts <- strsplit(lf, "\\|")[[1]]
    last_p <- parts[length(parts)]
    suppressWarnings(as.integer(last_p))
  }, integer(1), USE.NAMES = FALSE)
  
  if (any(is.na(parsed_len))) {
    bad_idx <- which(is.na(parsed_len))[1]
    stop(sprintf("Assignments file for sample '%s' has malformed length field '%s' at line %d",
                 sample_id, raw_reads$len_field[bad_idx], bad_idx), call. = FALSE)
  }
  raw_reads$read_length <- parsed_len
  
  # Effective classification: taxid > 0
  raw_reads$effective_classified <- (raw_reads$taxid > 0)
  
  n_total <- nrow(raw_reads)
  n_eff_class <- sum(raw_reads$effective_classified)
  n_eff_unclass <- n_total - n_eff_class
  
  # Reconcile against abundance expectations
  if (!is.null(expected_total) && n_total != expected_total) {
    stop(sprintf(
      "Reconciliation error for sample '%s': assignment rows (%d) != abundance total reads (%d)",
      sample_id, n_total, expected_total
    ), call. = FALSE)
  }
  if (!is.null(expected_classified) && n_eff_class != expected_classified) {
    stop(sprintf(
      "Reconciliation error for sample '%s': effective classified reads (%d) != abundance classified reads (%d)",
      sample_id, n_eff_class, expected_classified
    ), call. = FALSE)
  }
  if (!is.null(expected_unclassified) && n_eff_unclass != expected_unclassified) {
    stop(sprintf(
      "Reconciliation error for sample '%s': effective unclassified reads (%d) != abundance unclassified reads (%d)",
      sample_id, n_eff_unclass, expected_unclassified
    ), call. = FALSE)
  }
  
  raw_reads
}

read_metadata_table <- function(path, selected_samples) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    stop(sprintf("Metadata file not found: '%s'", path), call. = FALSE)
  }
  
  meta <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!"SampleID" %in% colnames(meta)) {
    stop("Metadata table must contain a 'SampleID' column.", call. = FALSE)
  }
  if (!"Group" %in% colnames(meta)) {
    stop("Metadata table must contain a 'Group' column.", call. = FALSE)
  }
  
  if (any(duplicated(meta$SampleID))) {
    dup_ids <- meta$SampleID[duplicated(meta$SampleID)]
    stop(sprintf("Metadata contains duplicate SampleID values: %s", paste(unique(dup_ids), collapse = ", ")), call. = FALSE)
  }
  
  meta_ids <- meta$SampleID
  missing_in_meta <- setdiff(selected_samples, meta_ids)
  extra_in_meta <- setdiff(meta_ids, selected_samples)
  
  if (length(missing_in_meta) > 0 || length(extra_in_meta) > 0) {
    msg <- "Metadata SampleID mismatch with selected abundance samples:"
    if (length(missing_in_meta) > 0) {
      msg <- paste0(msg, sprintf("\n  Missing from metadata: %s", paste(missing_in_meta, collapse = ", ")))
    }
    if (length(extra_in_meta) > 0) {
      msg <- paste0(msg, sprintf("\n  Extra in metadata: %s", paste(extra_in_meta, collapse = ", ")))
    }
    stop(msg, call. = FALSE)
  }
  
  # Align metadata to exact order of selected_samples
  meta_aligned <- meta[match(selected_samples, meta$SampleID), , drop = FALSE]
  rownames(meta_aligned) <- selected_samples
  meta_aligned
}

build_context <- function(cfg) {
  # 1. Read abundance table
  ab_res <- read_abundance_table(
    path = cfg$input$abundance_table,
    tax_col = cfg$input$tax_column,
    aggregate_cols = cfg$input$aggregate_columns,
    include_samples = cfg$input$include_samples
  )
  
  selected_samples <- ab_res$samples
  count_matrix <- ab_res$count_matrix
  taxonomy_df <- ab_res$taxonomy
  unclass_idx <- ab_res$unclass_index
  
  # Calculate per-sample read stats
  sample_stats <- data.frame(
    SampleID = selected_samples,
    TotalReads = colSums(count_matrix),
    UnclassifiedReads = count_matrix[unclass_idx, selected_samples],
    ClassifiedReads = colSums(count_matrix[-unclass_idx, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  
  # 2. Mode resolution
  configured_mode <- cfg$mode
  resolved_mode <- if (configured_mode == "auto") {
    if (length(selected_samples) == 1) "single" else "cohort"
  } else if (configured_mode %in% c("single", "cohort")) {
    configured_mode
  } else {
    stop(sprintf("Invalid mode '%s' in configuration. Must be 'auto', 'single', or 'cohort'.", configured_mode), call. = FALSE)
  }
  
  if (resolved_mode == "single" && length(selected_samples) != 1) {
    stop(sprintf("Mode is 'single' but %d samples are selected.", length(selected_samples)), call. = FALSE)
  }
  
  # 3. Read metadata
  metadata <- read_metadata_table(cfg$input$metadata, selected_samples)
  
  if (resolved_mode == "cohort" && is.null(metadata)) {
    stop("Cohort mode requires a metadata table mapping SampleID to Group.", call. = FALSE)
  }
  
  # 4. Assignments mapping
  assignments_map <- cfg$input$assignments
  if (!is.null(assignments_map) && !is.list(assignments_map)) {
    stop("Config 'input.assignments' must be a mapping of SampleID -> path or null.", call. = FALSE)
  }
  
  # 5. Read params.json if available
  params <- NULL
  if (!is.null(cfg$input$params_json) && file.exists(cfg$input$params_json)) {
    params <- suppressWarnings(tryCatch(jsonlite::fromJSON(cfg$input$params_json), error = function(e) NULL))
  }
  
  # 6. File hashes
  file_hashes <- list(
    abundance_table = compute_file_hash(cfg$input$abundance_table),
    metadata = compute_file_hash(cfg$input$metadata),
    params_json = compute_file_hash(cfg$input$params_json),
    taxonomy_cache = compute_file_hash(cfg$taxonomy$cache)
  )
  if (!is.null(assignments_map)) {
    for (s in names(assignments_map)) {
      file_hashes[[paste0("assignment_", s)]] <- compute_file_hash(assignments_map[[s]])
    }
  }
  
  list(
    config = cfg,
    mode = resolved_mode,
    samples = selected_samples,
    count_matrix = count_matrix,
    taxonomy = taxonomy_df,
    unclass_index = unclass_idx,
    sample_stats = sample_stats,
    metadata = metadata,
    assignments = assignments_map,
    params = params,
    file_hashes = file_hashes,
    warnings = character(0)
  )
}
