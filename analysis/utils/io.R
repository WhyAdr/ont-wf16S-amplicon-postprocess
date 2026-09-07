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

SUPPORTED_NCBI_DATABASE_SETS <- c("ncbi_16s_18s", "ncbi_16s_18s_28s_ITS")
BAMSTATS_REQUIRED_COLUMNS <- c("name", "sample_name", "iden", "ref_coverage")

read_upstream_params <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("input.params_json is required and must identify an existing wf-16s params.json.",
         call. = FALSE)
  }
  params <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) stop(sprintf("Could not parse params.json '%s': %s", path, e$message),
                             call. = FALSE)
  )
  required <- c("classifier", "database_set", "taxonomic_rank", "min_len", "max_len",
                "min_read_qual", "min_percent_identity", "min_ref_coverage",
                "abundance_threshold")
  missing <- required[vapply(required, function(x) is.null(params[[x]]), logical(1))]
  if (length(missing)) {
    stop(sprintf("params.json is missing required contract field(s): %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  scalar_string <- function(field) {
    value <- params[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop(sprintf("params.json field '%s' must be one non-empty string.", field), call. = FALSE)
    }
    value
  }
  classifier <- scalar_string("classifier")
  database_set <- scalar_string("database_set")
  taxonomic_rank <- scalar_string("taxonomic_rank")
  if (!identical(classifier, "minimap2")) {
    stop(sprintf(
      "Unsupported wf-16s classifier '%s'. This pipeline supports minimap2 only; Kraken2/Bracken requires a classifier-specific denominator model.",
      classifier
    ), call. = FALSE)
  }
  if (!database_set %in% SUPPORTED_NCBI_DATABASE_SETS) {
    stop(sprintf(
      "Unsupported wf-16s database_set '%s'. This pipeline supports the bundled NCBI database sets only.",
      database_set
    ), call. = FALSE)
  }
  if (!identical(taxonomic_rank, "S")) {
    stop(sprintf("Unsupported wf-16s taxonomic_rank '%s'; expected species rank 'S'.",
                 taxonomic_rank), call. = FALSE)
  }
  override_fields <- c("taxonomy", "reference", "ref2taxid", "database")
  active_overrides <- override_fields[vapply(override_fields, function(field) {
    value <- params[[field]]
    !is.null(value) && length(value) > 0L && !all(is.na(value)) && any(nzchar(as.character(value)))
  }, logical(1))]
  if (length(active_overrides)) {
    stop(sprintf("Custom wf-16s reference/taxonomy overrides are unsupported: %s",
                 paste(active_overrides, collapse = ", ")), call. = FALSE)
  }
  numeric_fields <- c("min_len", "max_len", "min_read_qual", "min_percent_identity",
                      "min_ref_coverage", "abundance_threshold")
  for (field in numeric_fields) {
    value <- params[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
      stop(sprintf("params.json field '%s' must be one finite number.", field), call. = FALSE)
    }
  }
  integer_fields <- c("min_len", "max_len", "abundance_threshold")
  for (field in integer_fields) {
    value <- params[[field]]
    if (abs(value - round(value)) > sqrt(.Machine$double.eps) ||
        value > .Machine$integer.max) {
      stop(sprintf(
        "params.json field '%s' must be a whole number no greater than %d.",
        field, .Machine$integer.max
      ), call. = FALSE)
    }
  }
  if (params$min_len <= 0 || params$max_len <= params$min_len || params$min_read_qual < 0 ||
      params$min_percent_identity < 0 || params$min_percent_identity > 100 ||
      params$min_ref_coverage < 0 || params$min_ref_coverage > 100 ||
      params$abundance_threshold < 0) {
    stop("params.json contains an invalid length, quality, identity, coverage, or abundance threshold.",
         call. = FALSE)
  }
  params
}

extract_upstream_contract <- function(params) {
  database_meta <- params$database_sets[[params$database_set]]
  list(
    workflow_name = "epi2me-labs/wf-16s",
    workflow_version = NULL,
    workflow_revision = NULL,
    wf_agent = params$wf$agent %||% NULL,
    classifier = params$classifier,
    database_set = params$database_set,
    taxonomy_namespace = "NCBI",
    database_taxonomy_source = database_meta$taxonomy %||% NULL,
    taxonomic_rank = params$taxonomic_rank,
    min_len = params$min_len,
    max_len = params$max_len,
    min_read_qual = params$min_read_qual,
    min_percent_identity = params$min_percent_identity,
    min_ref_coverage = params$min_ref_coverage,
    abundance_threshold = params$abundance_threshold,
    output_unclassified = params$output_unclassified %||% NULL,
    include_read_assignments = params$include_read_assignments %||% NULL
  )
}

read_bamstats_table <- function(path) {
  header_con <- gzfile(path, open = "rt")
  on.exit(close(header_con), add = TRUE)
  header_line <- readLines(header_con, n = 1L, warn = FALSE)
  if (length(header_line) == 0L) {
    stop(sprintf("Invalid bamstats schema or empty file: '%s'", path), call. = FALSE)
  }
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  if (anyDuplicated(header)) {
    stop(sprintf("Bamstats file '%s' contains duplicate column names.", path), call. = FALSE)
  }
  missing <- setdiff(BAMSTATS_REQUIRED_COLUMNS, header)
  if (length(missing) > 0L) {
    stop(sprintf("Bamstats file '%s' lacks required column(s): %s.", path,
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  col_classes <- stats::setNames(rep("NULL", length(header)), header)
  col_classes[BAMSTATS_REQUIRED_COLUMNS] <- "character"
  stats <- read.delim(gzfile(path), header = TRUE, sep = "\t", colClasses = col_classes,
                      check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(stats) == 0L) {
    stop(sprintf("Invalid bamstats schema or empty file: '%s'", path), call. = FALSE)
  }
  stats
}

discover_bamstats <- function(root, sample_ids) {
  mapped <- stats::setNames(rep(NA_character_, length(sample_ids)), sample_ids)
  if (is.null(root)) return(mapped)
  if (!dir.exists(root)) {
    stop(sprintf("Configured input.wf16s_output_root does not exist: '%s'", root), call. = FALSE)
  }
  candidates <- sort(list.files(
    root, pattern = "^bamstats[.]readstats[.]tsv[.]gz$",
    recursive = TRUE, full.names = TRUE
  ))
  if (!length(candidates)) {
    message(sprintf("[INFO] No bamstats.readstats.tsv.gz found under '%s'.", root))
    return(mapped)
  }
  for (path in candidates) {
    stats <- read_bamstats_table(path)
    samples_in_file <- unique(stats$sample_name)
    if (length(samples_in_file) == 0L || anyNA(samples_in_file) || any(!nzchar(trimws(samples_in_file)))) {
      stop(sprintf("Bamstats file '%s' contains missing or empty sample_name values.", path), call. = FALSE)
    }
    if (any(samples_in_file != trimws(samples_in_file))) {
      stop(sprintf("Bamstats file '%s' contains sample_name values with leading/trailing whitespace.",
                   path), call. = FALSE)
    }
    if (length(samples_in_file) > 1L) {
      stop(sprintf("Bamstats file '%s' contains inconsistent sample names: %s",
                   path, paste(samples_in_file, collapse = ", ")), call. = FALSE)
    }
    file_sample <- samples_in_file[[1]]
    if (!file_sample %in% sample_ids) next
    if (!is.na(mapped[[file_sample]])) {
      stop(sprintf("Multiple bamstats files discovered for sample '%s'.", file_sample), call. = FALSE)
    }
    mapped[[file_sample]] <- normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  missing <- names(mapped)[is.na(mapped)]
  if (length(missing)) {
    message(sprintf("[INFO] No bamstats file mapped for sample(s): %s",
                    paste(missing, collapse = ", ")))
  }
  mapped
}

partition_minimap2_failures <- function(reads, bamstats_path, params, sample_id) {
  stats <- read_bamstats_table(bamstats_path)
  if (anyNA(stats$name) || any(!nzchar(stats$name)) || anyDuplicated(stats$name)) {
    stop(sprintf("Bamstats read names for '%s' must be non-empty and unique.", sample_id),
         call. = FALSE)
  }
  sample_col <- as.character(stats$sample_name)
  if (anyNA(sample_col) || any(!nzchar(trimws(sample_col))) || any(sample_col != sample_id)) {
    stop(sprintf("Bamstats sample_name does not consistently equal '%s'.", sample_id),
         call. = FALSE)
  }
  c_reads <- reads[reads$status == "C", , drop = FALSE]
  matched <- match(c_reads$read_id, stats$name)
  if (anyNA(matched)) {
    stop(sprintf("Bamstats is missing %d status-C read(s) for '%s'.",
                 sum(is.na(matched)), sample_id), call. = FALSE)
  }
  aligned <- stats[matched, , drop = FALSE]
  aligned$iden <- suppressWarnings(as.numeric(aligned$iden))
  aligned$ref_coverage <- suppressWarnings(as.numeric(aligned$ref_coverage))
  if (any(!is.finite(aligned$iden)) || any(!is.finite(aligned$ref_coverage))) {
    stop(sprintf("Bamstats identity/coverage values for '%s' must be finite numbers.", sample_id),
         call. = FALSE)
  }
  identity_failed <- aligned$iden < params$min_percent_identity
  coverage_failed <- aligned$ref_coverage < params$min_ref_coverage
  positive <- c_reads$taxid > 0
  c0 <- !positive
  if (any(identity_failed[positive] | coverage_failed[positive])) {
    stop(sprintf("At least one TaxID>0 read for '%s' fails the recorded thresholds.", sample_id),
         call. = FALSE)
  }
  if (any(!identity_failed[c0] & !coverage_failed[c0])) {
    stop(sprintf("At least one C+TaxID0 read for '%s' passes both recorded thresholds.", sample_id),
         call. = FALSE)
  }
  list(
    matched = sum(c0),
    identity_only = sum(c0 & identity_failed & !coverage_failed),
    coverage_only = sum(c0 & !identity_failed & coverage_failed),
    both = sum(c0 & identity_failed & coverage_failed)
  )
}

validate_sample_ids <- function(sample_ids) {
  if (length(sample_ids) == 0) {
    stop("Abundance table validation error: No sample columns detected.", call. = FALSE)
  }

  if (any(is.na(sample_ids)) || any(sample_ids == "")) {
    stop("Sample ID validation error: Empty or NA sample ID detected.", call. = FALSE)
  }

  if (any(sample_ids != trimws(sample_ids))) {
    stop("Sample ID validation error: Sample IDs must not have leading or trailing whitespace.",
         call. = FALSE)
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
  if (any(grepl("[.]$", sanitized))) {
    stop("Sample ID validation error: Portable output basenames must not end in a dot.",
         call. = FALSE)
  }
  if (anyDuplicated(tolower(sanitized))) {
    stop("Sample ID validation error: Sample IDs collide after portable filename normalization.",
         call. = FALSE)
  }
  windows_base <- toupper(sub("[.].*$", "", sanitized))
  reserved <- windows_base %in% c("CON", "PRN", "AUX", "NUL",
                                  paste0("COM", 1:9), paste0("LPT", 1:9))
  if (any(reserved)) {
    stop(sprintf("Sample ID validation error: '%s' is a reserved Windows device basename.",
                 sample_ids[which(reserved)[1]]), call. = FALSE)
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
      stated_total <- suppressWarnings(as.numeric(raw_df[[ac]]))
      if (anyNA(stated_total) || any(!is.finite(stated_total)) || any(stated_total < 0) ||
          any(abs(stated_total - round(stated_total)) > 1e-6)) {
        stop(sprintf("Aggregate column '%s' must contain finite, non-negative integer counts.", ac),
             call. = FALSE)
      }
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
  if (anyNA(lineages) || any(!nzchar(lineages)) || any(lineages != trimws(lineages))) {
    stop("Taxonomy lineages must be non-empty and must not have leading/trailing whitespace.", call. = FALSE)
  }
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
  for (row_index in seq_along(parsed_lineages)) {
    fields <- parsed_lineages[[row_index]]
    bad_rank <- which(!nzchar(trimws(fields)) | fields != trimws(fields))
    if (length(bad_rank) > 0L) {
      rank_index <- bad_rank[1]
      stop(sprintf(
        "Lineage schema violation at row %d: rank '%s' must be non-empty and have no leading/trailing whitespace ('%s').",
        row_index + 1L, RANKS_8[rank_index], lineages[row_index]
      ), call. = FALSE)
    }
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
  expected_unclassified <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
  if (!identical(lineages[unclass_indices], expected_unclassified)) {
    stop(sprintf("Unclassified lineage must be exactly '%s'.", expected_unclassified), call. = FALSE)
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

open_assignments_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

split_assignment_fields <- function(line) {
  # Appending a sentinel tab makes strsplit() retain physical trailing empty
  # fields, which is required for exact five-column validation.
  strsplit(paste0(line, "\t"), "\t", fixed = TRUE)[[1]]
}

parse_assignment_length <- function(length_field) {
  if (is.na(length_field) || length_field == "") return(NA_integer_)
  parts <- strsplit(length_field, "\\|")[[1]]
  last_part <- parts[length(parts)]
  if (!grepl("^[0-9]+$", last_part)) return(NA_integer_)
  suppressWarnings(as.integer(last_part))
}

validate_assignment_chunk <- function(lines, sample_id, line_numbers, seen_read_ids) {
  fields <- strsplit(paste0(lines, "\t"), "\t", fixed = TRUE)
  field_counts <- lengths(fields)
  if (any(field_counts != 5L)) {
    bad_index <- which(field_counts != 5L)[1]
    stop(sprintf(
      "Assignments file for sample '%s' has %d fields at line %d; expected exactly 5.",
      sample_id, field_counts[bad_index], line_numbers[bad_index]
    ), call. = FALSE)
  }

  fields <- do.call(rbind, fields)
  status <- fields[, 1]
  read_id <- fields[, 2]
  taxid_text <- fields[, 3]
  length_field <- fields[, 4]

  if (any(!nzchar(read_id))) {
    bad_index <- which(!nzchar(read_id))[1]
    stop(sprintf("Assignments file for sample '%s' has an empty read ID at line %d",
                 sample_id, line_numbers[bad_index]), call. = FALSE)
  }
  already_seen <- vapply(read_id, exists, logical(1), envir = seen_read_ids, inherits = FALSE)
  duplicate <- duplicated(read_id) | already_seen
  if (any(duplicate)) {
    duplicate_id <- read_id[which(duplicate)[1]]
    stop(sprintf("Assignments file for sample '%s' contains duplicate read ID: '%s'",
                 sample_id, duplicate_id), call. = FALSE)
  }
  list2env(stats::setNames(as.list(rep(TRUE, length(read_id))), read_id), envir = seen_read_ids)

  valid_status <- status %in% c("C", "U")
  if (any(!valid_status)) {
    bad_index <- which(!valid_status)[1]
    stop(sprintf("Assignments file for sample '%s' has invalid status '%s' at line %d",
                 sample_id, status[bad_index], line_numbers[bad_index]), call. = FALSE)
  }

  valid_taxid <- grepl("^[0-9]+$", taxid_text)
  taxid <- suppressWarnings(as.numeric(taxid_text))
  invalid_taxid <- !valid_taxid | is.na(taxid) | !is.finite(taxid)
  if (any(invalid_taxid)) {
    bad_index <- which(invalid_taxid)[1]
    stop(sprintf("Assignments file for sample '%s' has non-integer TaxID '%s' at line %d",
                 sample_id, taxid_text[bad_index], line_numbers[bad_index]), call. = FALSE)
  }
  inconsistent <- status == "U" & taxid > 0
  if (any(inconsistent)) {
    bad_index <- which(inconsistent)[1]
    stop(sprintf("Assignments file for sample '%s' has status U with positive TaxID at line %d",
                 sample_id, line_numbers[bad_index]), call. = FALSE)
  }

  read_length <- vapply(length_field, parse_assignment_length, integer(1), USE.NAMES = FALSE)
  invalid_length <- is.na(read_length) | !is.finite(read_length) | read_length <= 0
  if (any(invalid_length)) {
    bad_index <- which(invalid_length)[1]
    stop(sprintf("Assignments file for sample '%s' has malformed length field '%s' at line %d",
                 sample_id, length_field[bad_index], line_numbers[bad_index]), call. = FALSE)
  }

  invisible(NULL)
}

scan_assignments_file <- function(path, sample_id, chunk_size) {
  connection <- open_assignments_connection(path)
  on.exit(close(connection), add = TRUE)
  seen_read_ids <- new.env(hash = TRUE, parent = emptyenv())
  row_count <- 0L
  physical_line <- 0L

  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) break
    line_numbers <- physical_line + seq_along(lines)
    physical_line <- physical_line + length(lines)
    nonempty <- nzchar(lines)
    if (!any(nonempty)) next
    validate_assignment_chunk(lines[nonempty], sample_id, line_numbers[nonempty], seen_read_ids)
    row_count <- row_count + sum(nonempty)
  }

  if (row_count == 0L) {
    stop(sprintf("Assignments file for sample '%s' is empty.", sample_id), call. = FALSE)
  }
  row_count
}

read_assignments_file <- function(path, sample_id, expected_total = NULL, expected_classified = NULL,
                                  expected_unclassified = NULL, chunk_size = 100000L) {
  if (!file.exists(path)) {
    stop(sprintf("Assignments file for sample '%s' not found: '%s'", sample_id, path), call. = FALSE)
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) ||
      chunk_size < 1 || chunk_size != as.integer(chunk_size)) {
    stop("'chunk_size' must be a positive integer.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)
  row_count <- scan_assignments_file(path, sample_id, chunk_size)

  # The validation pass avoids retaining raw text; this pass only builds the
  # read-level data frame still required by QC downstream.
  raw_reads <- data.frame(
    status = character(row_count),
    read_id = character(row_count),
    taxid = numeric(row_count),
    len_field = character(row_count),
    lineage = character(row_count),
    read_length = integer(row_count),
    stringsAsFactors = FALSE
  )
  connection <- open_assignments_connection(path)
  on.exit(close(connection), add = TRUE)
  row_index <- 0L
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) break
    nonempty <- nzchar(lines)
    if (!any(nonempty)) next
    fields <- do.call(rbind, lapply(lines[nonempty], split_assignment_fields))
    indices <- seq.int(row_index + 1L, length.out = nrow(fields))
    raw_reads$status[indices] <- fields[, 1]
    raw_reads$read_id[indices] <- fields[, 2]
    raw_reads$taxid[indices] <- suppressWarnings(as.numeric(fields[, 3]))
    raw_reads$len_field[indices] <- fields[, 4]
    raw_reads$lineage[indices] <- fields[, 5]
    raw_reads$read_length[indices] <- vapply(fields[, 4], parse_assignment_length,
                                              integer(1), USE.NAMES = FALSE)
    row_index <- row_index + nrow(fields)
  }
  raw_reads$effective_classified <- raw_reads$taxid > 0

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

  header_line <- readLines(path, n = 1L, warn = FALSE)
  if (length(header_line) == 0L) {
    stop("Metadata table is empty.", call. = FALSE)
  }
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  if (anyDuplicated(header)) {
    stop("Metadata table contains duplicate column names.", call. = FALSE)
  }
  required_identity <- c("SampleID", "Group")
  missing_identity <- setdiff(required_identity, header)
  if (length(missing_identity) > 0L) {
    stop(sprintf("Metadata table must contain column(s): %s.",
                 paste(missing_identity, collapse = ", ")), call. = FALSE)
  }
  col_classes <- rep(NA_character_, length(header))
  col_classes[match(required_identity, header)] <- "character"
  meta <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     check.names = FALSE, colClasses = col_classes)

  invalid_sample <- is.na(meta$SampleID) | !nzchar(trimws(meta$SampleID))
  invalid_group <- is.na(meta$Group) | !nzchar(trimws(meta$Group))
  if (any(invalid_sample)) {
    stop(sprintf("Metadata contains an empty SampleID at data row %d.", which(invalid_sample)[1]), call. = FALSE)
  }
  if (any(invalid_group)) {
    stop(sprintf("Metadata contains an empty Group at data row %d.", which(invalid_group)[1]), call. = FALSE)
  }
  if (any(meta$SampleID != trimws(meta$SampleID)) || any(meta$Group != trimws(meta$Group))) {
    stop("Metadata SampleID and Group values must not have leading or trailing whitespace.", call. = FALSE)
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
  # 1. Fail closed on the producer contract before parsing classifier-specific files.
  params <- read_upstream_params(cfg$input$params_json)
  upstream_contract <- extract_upstream_contract(params)

  # 2. Read abundance table
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

  # 3. Mode resolution
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
  if (resolved_mode == "cohort" && length(selected_samples) < 2L) {
    stop(sprintf("Mode is 'cohort' but only %d sample is selected; at least 2 are required.",
                 length(selected_samples)), call. = FALSE)
  }

  # 4. Read metadata
  metadata <- read_metadata_table(cfg$input$metadata, selected_samples)

  if (resolved_mode == "cohort" && is.null(metadata)) {
    stop("Cohort mode requires a metadata table mapping SampleID to Group.", call. = FALSE)
  }

  # 5. Assignments mapping
  assignments_map <- cfg$input$assignments
  if (!is.null(assignments_map) && !is.list(assignments_map)) {
    stop("Config 'input.assignments' must be a mapping of SampleID -> path or null.", call. = FALSE)
  }
  assignment_data <- list()
  retain_assignment_rows <- !is.null(cfg$cli) && "qc" %in% cfg$cli$modules && !isTRUE(cfg$cli$validate_only)
  if (!is.null(assignments_map)) {
    extra_assignment_ids <- setdiff(names(assignments_map), selected_samples)
    if (length(extra_assignment_ids) > 0L) {
      stop(sprintf("Assignments configured for unselected/unknown sample(s): %s",
                   paste(extra_assignment_ids, collapse = ", ")), call. = FALSE)
    }
    for (s in names(assignments_map)) {
      stat_row <- sample_stats[sample_stats$SampleID == s, , drop = FALSE]
      parsed_assignments <- read_assignments_file(
        assignments_map[[s]], s,
        expected_total = stat_row$TotalReads[1],
        expected_classified = stat_row$ClassifiedReads[1],
        expected_unclassified = stat_row$UnclassifiedReads[1]
      )
      if (retain_assignment_rows) assignment_data[[s]] <- parsed_assignments
      rm(parsed_assignments)
    }
  }

  if (resolved_mode == "cohort" && !is.null(cfg$beta$strata_column)) {
    strata_col <- cfg$beta$strata_column
    if (!strata_col %in% colnames(metadata)) {
      stop(sprintf("Configured beta.strata_column '%s' is absent from metadata.", strata_col), call. = FALSE)
    }
    strata <- metadata[[strata_col]]
    if (anyNA(strata) || any(!nzchar(trimws(as.character(strata))))) {
      stop(sprintf("Metadata strata column '%s' contains missing/empty values.", strata_col), call. = FALSE)
    }
    if (length(unique(strata)) < 2L) {
      stop(sprintf("Metadata strata column '%s' must contain at least two strata.", strata_col), call. = FALSE)
    }
  }

  # 6. Optional per-sample bamstats discovery
  bamstats <- discover_bamstats(cfg$input$wf16s_output_root, selected_samples)

  # 7. File hashes
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
  for (s in names(bamstats)[!is.na(bamstats)]) {
    file_hashes[[paste0("bamstats_", s)]] <- compute_file_hash(bamstats[[s]])
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
    assignment_data = assignment_data,
    bamstats = bamstats,
    params = params,
    upstream_contract = upstream_contract,
    file_hashes = file_hashes,
    warnings = character(0)
  )
}
