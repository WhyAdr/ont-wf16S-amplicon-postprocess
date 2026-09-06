# =============================================================================
# Module 08: FAPROTAX 1.2.12 Taxon-Based Functional Inference
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(scales)
})

FAPROTAX_EXPECTED_VERSION <- "1.2.12"

get_faprotax_database <- function() {
  env <- new.env(parent = emptyenv())
  utils::data("prok_func_FAPROTAX", package = "microeco", envir = env)
  db <- env[["prok_func_FAPROTAX"]]
  if (!is.list(db) || is.null(db[["ver"]]) || length(db[["ver"]]) != 1L) {
    stop("microeco did not expose embedded FAPROTAX version metadata.", call. = FALSE)
  }
  db
}

validate_faprotax_runtime <- function() {
  if (!require("microeco", quietly = TRUE, character.only = TRUE)) {
    stop("Module 'faprotax' requires the optional package 'microeco'.", call. = FALSE)
  }

  min_version <- if (exists("MICROECO_MIN_VERSION", inherits = TRUE)) {
    MICROECO_MIN_VERSION
  } else {
    "2.3.0"
  }
  installed_version <- utils::packageVersion("microeco")
  if (utils::compareVersion(as.character(installed_version), min_version) < 0) {
    stop(sprintf("Module 'faprotax' requires microeco >= %s.", min_version), call. = FALSE)
  }

  db <- get_faprotax_database()
  observed_version <- as.character(db[["ver"]])
  if (!identical(observed_version, FAPROTAX_EXPECTED_VERSION)) {
    stop(sprintf("Expected FAPROTAX %s, but microeco embeds %s.",
                 FAPROTAX_EXPECTED_VERSION, observed_version), call. = FALSE)
  }

  list(
    microeco_version = as.character(installed_version),
    faprotax_version = observed_version
  )
}

prepare_faprotax_input <- function(context) {
  count_matrix <- context$count_matrix[, context$samples, drop = FALSE]
  classified_indices <- setdiff(seq_len(nrow(count_matrix)), context$unclass_index)
  if (length(classified_indices) == 0L) return(NULL)

  source_taxonomy <- context$taxonomy[classified_indices, , drop = FALSE]
  source_counts <- count_matrix[classified_indices, , drop = FALSE]
  positive_prokaryotic <- !is.na(source_taxonomy$superkingdom) &
    source_taxonomy$superkingdom %in% c("Bacteria", "Archaea") &
    rowSums(source_counts) > 0
  if (!any(positive_prokaryotic)) return(NULL)

  source_taxonomy <- source_taxonomy[positive_prokaryotic, , drop = FALSE]
  source_counts <- source_counts[positive_prokaryotic, , drop = FALSE]
  source_indices <- classified_indices[positive_prokaryotic]
  feature_ids <- sprintf("Taxon_%06d", source_indices)
  if (anyDuplicated(feature_ids)) {
    stop("FAPROTAX feature identifiers are not unique.", call. = FALSE)
  }
  rownames(source_counts) <- feature_ids

  # microeco expects Kingdom to contain Bacteria/Archaea. Deliberately omit
  # the intervening NCBI kingdom rank (Bacillati, Pseudomonadati, etc.).
  tax_table <- data.frame(
    Kingdom = as.character(source_taxonomy$superkingdom),
    Phylum = as.character(source_taxonomy$phylum),
    Class = as.character(source_taxonomy$class),
    Order = as.character(source_taxonomy$order),
    Family = as.character(source_taxonomy$family),
    Genus = as.character(source_taxonomy$genus),
    Species = as.character(source_taxonomy$species),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  unknown_cells <- is.na(tax_table) |
    tolower(trimws(as.matrix(tax_table))) %in% c("unknown", "unclassified")
  tax_table[unknown_cells] <- ""
  rownames(tax_table) <- feature_ids

  sample_table <- context$metadata
  if (is.null(sample_table)) {
    sample_table <- data.frame(
      SampleID = context$samples,
      Group = context$samples,
      row.names = context$samples,
      stringsAsFactors = FALSE
    )
  } else {
    sample_table <- sample_table[context$samples, , drop = FALSE]
    rownames(sample_table) <- context$samples
  }

  list(
    counts = source_counts,
    tax_table = tax_table,
    source_taxonomy = source_taxonomy,
    sample_table = sample_table,
    feature_ids = feature_ids,
    source_indices = source_indices,
    samples = context$samples
  )
}

as_integer_valued_numeric <- function(x, field) {
  x <- as.numeric(x)
  if (anyNA(x) || any(!is.finite(x)) || any(x < 0) ||
      any(abs(x - round(x)) > 1e-6)) {
    stop(sprintf("FAPROTAX output field '%s' contains invalid read counts.", field),
         call. = FALSE)
  }
  round(x)
}

validate_faprotax_outputs <- function(function_abundance, coverage) {
  required_abundance <- c(
    "SampleID", "Function", "FunctionReadCount",
    "PctEligibleProkaryoticReads", "PctClassifiedReads", "PctTotalReads"
  )
  if (!identical(names(function_abundance), required_abundance)) {
    stop("FAPROTAX function abundance has an invalid column contract.", call. = FALSE)
  }
  required_coverage <- c(
    "SampleID", "TotalReads", "UnclassifiedReads", "ClassifiedReads",
    "EligibleProkaryoticReads", "ExcludedNonProkaryoticClassifiedReads",
    "FunctionMappedReads", "FunctionUnmappedEligibleReads",
    "EligibleProkaryoticTaxa", "FunctionMappedTaxa"
  )
  if (!identical(names(coverage), required_coverage)) {
    stop("FAPROTAX mapping coverage has an invalid column contract.", call. = FALSE)
  }

  count_columns <- setdiff(required_coverage, "SampleID")
  for (field in count_columns) {
    coverage[[field]] <- as_integer_valued_numeric(coverage[[field]], field)
  }
  if (anyDuplicated(coverage$SampleID) || anyNA(coverage$SampleID)) {
    stop("FAPROTAX mapping coverage must contain unique sample IDs.", call. = FALSE)
  }

  if (any(coverage$ClassifiedReads + coverage$UnclassifiedReads != coverage$TotalReads)) {
    stop("FAPROTAX read accounting failed: classified plus unclassified does not equal total.",
         call. = FALSE)
  }
  if (any(coverage$EligibleProkaryoticReads +
          coverage$ExcludedNonProkaryoticClassifiedReads != coverage$ClassifiedReads)) {
    stop("FAPROTAX read accounting failed: eligible plus excluded does not equal classified.",
         call. = FALSE)
  }
  if (any(coverage$FunctionMappedReads + coverage$FunctionUnmappedEligibleReads !=
          coverage$EligibleProkaryoticReads)) {
    stop("FAPROTAX read accounting failed: mapped plus unmapped does not equal eligible.",
         call. = FALSE)
  }

  if (nrow(function_abundance) > 0L) {
    if (!all(function_abundance$SampleID %in% coverage$SampleID)) {
      stop("FAPROTAX function abundance contains an unknown sample ID.", call. = FALSE)
    }
    eligible_by_sample <- stats::setNames(
      coverage$EligibleProkaryoticReads, coverage$SampleID
    )
    if (any(function_abundance$FunctionReadCount < 0) ||
        any(!is.finite(function_abundance$FunctionReadCount)) ||
        any(abs(function_abundance$FunctionReadCount -
               round(function_abundance$FunctionReadCount)) > 1e-6)) {
      stop("FAPROTAX function abundance contains invalid counts.", call. = FALSE)
    }
    if (any(function_abundance$FunctionReadCount >
            eligible_by_sample[function_abundance$SampleID] + 1e-6)) {
      stop("FAPROTAX function abundance exceeds eligible prokaryotic reads.", call. = FALSE)
    }
    pct_columns <- c(
      "PctEligibleProkaryoticReads", "PctClassifiedReads", "PctTotalReads"
    )
    for (field in pct_columns) {
      values <- function_abundance[[field]]
      if (anyNA(values) || any(!is.finite(values)) || any(values < -1e-8) ||
          any(values > 100 + 1e-8)) {
        stop(sprintf("FAPROTAX percentage field '%s' is outside [0, 100].", field),
             call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

build_faprotax_function_abundance <- function(function_counts, samples, coverage) {
  positive_functions <- rownames(function_counts)[rowSums(function_counts) > 0]
  result <- data.frame(
    SampleID = character(0),
    Function = character(0),
    FunctionReadCount = numeric(0),
    PctEligibleProkaryoticReads = numeric(0),
    PctClassifiedReads = numeric(0),
    PctTotalReads = numeric(0),
    stringsAsFactors = FALSE
  )
  if (length(positive_functions) == 0L) return(result)

  result <- expand.grid(
    Function = positive_functions,
    SampleID = samples,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  result$FunctionReadCount <- as.numeric(function_counts[
    cbind(match(result$Function, rownames(function_counts)),
          match(result$SampleID, colnames(function_counts)))
  ])

  coverage_by_sample <- coverage[match(result$SampleID, coverage$SampleID), , drop = FALSE]
  safe_pct <- function(numerator, denominator) {
    ifelse(denominator > 0, 100 * numerator / denominator, 0)
  }
  result$PctEligibleProkaryoticReads <- safe_pct(
    result$FunctionReadCount, coverage_by_sample$EligibleProkaryoticReads
  )
  result$PctClassifiedReads <- safe_pct(
    result$FunctionReadCount, coverage_by_sample$ClassifiedReads
  )
  result$PctTotalReads <- safe_pct(
    result$FunctionReadCount, coverage_by_sample$TotalReads
  )

  result <- result[, c(
    "SampleID", "Function", "FunctionReadCount",
    "PctEligibleProkaryoticReads", "PctClassifiedReads", "PctTotalReads"
  )]
  rownames(result) <- NULL
  result
}

build_faprotax_taxon_function_assignments <- function(binary, prepared) {
  mapped_taxa <- which(rowSums(binary) > 0)
  if (length(mapped_taxa) == 0L) {
    return(data.frame(
      FeatureID = character(0),
      Function = character(0),
      TaxonPath = character(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(mapped_taxa, function(i) {
    functions <- colnames(binary)[which(binary[i, ] == 1)]
    if (length(functions) == 0L) return(NULL)
    data.frame(
      FeatureID = rep(prepared$feature_ids[i], length(functions)),
      Function = functions,
      TaxonPath = rep(prepared$source_taxonomy$TaxonPath[i], length(functions)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      FeatureID = character(0),
      Function = character(0),
      TaxonPath = character(0),
      stringsAsFactors = FALSE
    ))
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

plot_faprotax_functions <- function(function_abundance, top_n, path) {
  if (nrow(function_abundance) == 0L) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No positive FAPROTAX assignments") +
      labs(title = "FAPROTAX functional inference") +
      theme_void()
    save_plot(path, p, width = 8, height = 5)
    return(invisible(path))
  }

  totals <- aggregate(
    FunctionReadCount ~ Function,
    data = function_abundance,
    FUN = sum
  )
  totals <- totals[order(-totals$FunctionReadCount, totals$Function), , drop = FALSE]
  selected <- head(totals$Function, top_n)
  plot_data <- function_abundance[function_abundance$Function %in% selected, , drop = FALSE]
  plot_data$Function <- factor(plot_data$Function, levels = rev(selected))

  p <- ggplot(plot_data, aes(x = Function, y = FunctionReadCount)) +
    geom_col(width = 0.72, fill = "#1b9e77") +
    facet_wrap(~SampleID, scales = "free_y") +
    coord_flip() +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = sprintf("Top %d FAPROTAX Functions", length(selected)),
      subtitle = paste(
        "Independent bars; function counts can overlap across functions.",
        "Counts are sums of classified taxon counts."
      ),
      x = NULL,
      y = "Function read count"
    ) +
    theme_amplicon()
  save_plot(path, p, width = 9, height = max(5, 4 + 0.18 * length(selected)))
  invisible(path)
}

run_faprotax <- function(context) {
  runtime <- validate_faprotax_runtime()
  out_dir <- context$config$output$dirs$faprotax
  if (is.null(out_dir) || !nzchar(out_dir)) {
    stop("FAPROTAX output directory is not configured.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  prepared <- prepare_faprotax_input(context)
  if (is.null(prepared)) {
    skip_file <- file.path(out_dir, "faprotax_skipped.tsv")
    write.table(
      data.frame(
        Status = "skipped",
        Reason = "No positive-count Bacteria or Archaea rows",
        stringsAsFactors = FALSE
      ),
      skip_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    return(list(
      status = "skipped",
      reason = "No eligible prokaryotic taxa",
      outputs = skip_file
    ))
  }

  dataset <- getElement(microeco::microtable, "new")(
    otu_table = as.data.frame(prepared$counts, check.names = FALSE),
    sample_table = prepared$sample_table,
    tax_table = prepared$tax_table,
    auto_tidy = FALSE
  )
  predictor <- getElement(microeco::trans_func, "new")(dataset)
  if (!identical(predictor$for_what, "prok")) {
    stop("microeco did not recognize the prepared taxonomy as prokaryotic.", call. = FALSE)
  }
  predictor$cal_func(prok_database = "FAPROTAX")
  raw_binary <- predictor$res_func
  if (is.null(raw_binary) || is.null(dim(raw_binary)) ||
      !all(prepared$feature_ids %in% rownames(raw_binary))) {
    stop("microeco returned a function matrix without the prepared feature IDs.", call. = FALSE)
  }
  binary <- as.matrix(raw_binary[prepared$feature_ids, , drop = FALSE])
  storage.mode(binary) <- "numeric"
  if (anyNA(binary) || any(!is.finite(binary)) || any(!binary %in% c(0, 1))) {
    stop("microeco returned a non-binary or missing taxon-function matrix.", call. = FALSE)
  }
  if (is.null(colnames(binary)) || anyNA(colnames(binary)) ||
      any(!nzchar(colnames(binary))) || anyDuplicated(colnames(binary))) {
    stop("microeco returned an invalid function-name contract.", call. = FALSE)
  }

  sample_stats <- context$sample_stats[match(prepared$samples, context$sample_stats$SampleID), , drop = FALSE]
  if (anyNA(sample_stats$SampleID)) {
    stop("FAPROTAX could not align sample statistics to the input samples.", call. = FALSE)
  }
  total_reads <- as.numeric(sample_stats$TotalReads)
  unclassified_reads <- as.numeric(sample_stats$UnclassifiedReads)
  classified_reads <- as.numeric(sample_stats$ClassifiedReads)
  eligible_reads <- as.numeric(colSums(prepared$counts))
  names(eligible_reads) <- prepared$samples
  excluded_nonprok <- classified_reads - eligible_reads

  mapped_taxa <- rowSums(binary) > 0
  function_mapped_reads <- if (any(mapped_taxa)) {
    as.numeric(colSums(prepared$counts[mapped_taxa, , drop = FALSE]))
  } else {
    rep(0, length(prepared$samples))
  }
  names(function_mapped_reads) <- prepared$samples
  function_unmapped_reads <- eligible_reads - function_mapped_reads
  eligible_taxa <- as.numeric(colSums(prepared$counts > 0))
  function_mapped_taxa <- if (any(mapped_taxa)) {
    as.numeric(colSums(prepared$counts[mapped_taxa, , drop = FALSE] > 0))
  } else {
    rep(0, length(prepared$samples))
  }

  coverage <- data.frame(
    SampleID = prepared$samples,
    TotalReads = total_reads,
    UnclassifiedReads = unclassified_reads,
    ClassifiedReads = classified_reads,
    EligibleProkaryoticReads = eligible_reads[prepared$samples],
    ExcludedNonProkaryoticClassifiedReads = excluded_nonprok,
    FunctionMappedReads = function_mapped_reads[prepared$samples],
    FunctionUnmappedEligibleReads = function_unmapped_reads[prepared$samples],
    EligibleProkaryoticTaxa = eligible_taxa,
    FunctionMappedTaxa = function_mapped_taxa,
    stringsAsFactors = FALSE
  )

  function_counts <- t(binary) %*% as.matrix(prepared$counts)
  function_counts <- as.matrix(function_counts)
  rownames(function_counts) <- colnames(binary)
  colnames(function_counts) <- prepared$samples
  function_abundance <- build_faprotax_function_abundance(
    function_counts, prepared$samples, coverage
  )
  taxon_function_assignments <- build_faprotax_taxon_function_assignments(binary, prepared)
  validate_faprotax_outputs(function_abundance, coverage)

  function_file <- file.path(out_dir, "faprotax_function_abundance.tsv")
  coverage_file <- file.path(out_dir, "faprotax_mapping_coverage.tsv")
  assignments_file <- file.path(out_dir, "faprotax_taxon_function_assignments.tsv")
  plot_file <- file.path(out_dir, "faprotax_top_functions.png")
  provenance_file <- file.path(out_dir, "faprotax_provenance.json")

  write.table(function_abundance, function_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(coverage, coverage_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(taxon_function_assignments, assignments_file,
              sep = "\t", row.names = FALSE, quote = FALSE)
  plot_faprotax_functions(
    function_abundance,
    top_n = context$config$faprotax$top_n_functions,
    path = plot_file
  )

  provenance <- c(runtime, list(
    method = "microeco::trans_func$cal_func(prok_database = \"FAPROTAX\")",
    database = "FAPROTAX",
    database_source = "microeco embedded data object prok_func_FAPROTAX",
    abundance_model = "sum of classified taxon counts carrying each function",
    input_filter = list(
      included_superkingdoms = c("Bacteria", "Archaea"),
      excluded_rows = c("canonical unclassified", "classified non-prokaryotic"),
      positive_count_rows_only = TRUE
    ),
    taxonomy_mapping = list(
      microeco_Kingdom = "pipeline superkingdom",
      omitted_pipeline_rank = "kingdom"
    ),
    denominator_definitions = list(
      TotalReads = "sum of all abundance-table rows, including canonical unclassified",
      UnclassifiedReads = "canonical unclassified abundance-table row",
      ClassifiedReads = "TotalReads minus UnclassifiedReads",
      EligibleProkaryoticReads = "positive-count Bacteria and Archaea rows",
      ExcludedNonProkaryoticClassifiedReads = "ClassifiedReads minus EligibleProkaryoticReads",
      FunctionMappedReads = "each eligible read counted once if its taxon has at least one function",
      FunctionUnmappedEligibleReads = "EligibleProkaryoticReads minus FunctionMappedReads",
      EligibleProkaryoticTaxa = "positive-count eligible taxon rows per sample",
      FunctionMappedTaxa = "positive-count eligible taxon rows with at least one function per sample"
    ),
    functions_are_nonexclusive = TRUE,
    overlap_warning = paste(
      "Functions overlap; percentages are independent fractions of their",
      "explicit read denominators and must not be stacked or expected to sum to 100%."
    ),
    interpretation = paste(
      "Taxon-based functional assignment/inference; this does not establish",
      "gene presence, pathway completeness, expression, or activity."
    ),
    concordance = list(
      official_collapse_table_review = "tracked fixture differences documented",
      review_document = "faprotax1.2.12-concordance-review.md",
      note = paste(
        "microeco output is not described as byte-for-byte equivalent to the",
        "official FAPROTAX collapse_table.py workflow; see the tracked-fixture",
        "review for the observed engine differences."
      )
    )
  ))
  jsonlite::write_json(provenance, provenance_file,
                       pretty = TRUE, auto_unbox = TRUE, null = "null")

  list(
    status = "completed",
    outputs = c(function_file, coverage_file, assignments_file, plot_file, provenance_file),
    microeco_version = runtime$microeco_version,
    faprotax_version = runtime$faprotax_version
  )
}
