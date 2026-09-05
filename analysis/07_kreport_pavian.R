# =============================================================================
# Module 07: Kraken Report (.kreport) Generation for Pavian Sankey
# =============================================================================

suppressMessages({
  library(jsonlite)
  library(dplyr)
  library(processx)
})

run_kreport <- function(context) {
  cfg <- context$config
  kreport_dir <- cfg$output$dirs$kreport
  dir.create(kreport_dir, recursive = TRUE, showWarnings = FALSE)

  all_outputs <- character(0)
  cache_file <- cfg$taxonomy$cache
  network_mode <- cfg$taxonomy$network_mode %||% "cache_only"
  unresolved_policy <- cfg$taxonomy$unresolved_policy %||% "warn"

  unresolved_tsv <- file.path(kreport_dir, "unresolved_taxids.tsv")
  conflicts_tsv <- file.path(kreport_dir, "taxonomy_conflicts.tsv")
  resolution_sources_tsv <- file.path(kreport_dir, "taxonomy_resolution_sources.tsv")
  prov_json <- file.path(kreport_dir, "taxonomy_provenance.json")
  resolved_cache <- file.path(kreport_dir, "resolved_taxonomy_cache.json")

  py_script <- file.path(cfg$pipeline_root, "analysis", "utils", "ncbi_taxonomy.py")
  if (!file.exists(py_script)) {
    stop(sprintf("Taxonomy resolver script not found: '%s'", py_script), call. = FALSE)
  }

  python_cmd <- find_python()

  cmd_args <- c(
    py_script,
    "--abundance", cfg$input$abundance_table,
    "--tax-column", cfg$input$tax_column,
    "--cache", cache_file,
    "--resolved-cache", resolved_cache,
    "--mode", network_mode,
    "--email-env", cfg$taxonomy$email_env,
    "--api-key-env", cfg$taxonomy$api_key_env,
    "--unresolved-policy", unresolved_policy,
    "--unresolved-tsv", unresolved_tsv,
    "--conflicts-tsv", conflicts_tsv,
    "--resolution-sources-tsv", resolution_sources_tsv,
    "--provenance", prov_json
  )
  assignment_paths <- unname(unlist(context$assignments, use.names = FALSE))
  if (length(assignment_paths) > 0L) {
    cmd_args <- c(cmd_args, as.vector(rbind("--assignments", assignment_paths)))
  }

  # processx passes a true argument vector on Windows and Unix; do not shell-quote.
  resolver <- processx::run(
    command = python_cmd,
    args = cmd_args,
    echo = TRUE,
    error_on_status = FALSE
  )
  if (!identical(resolver$status, 0L)) {
    stop(sprintf("NCBI taxonomy resolver failed with exit status %d", resolver$status),
         call. = FALSE)
  }

  required_resolver_outputs <- c(
    resolved_cache, unresolved_tsv, conflicts_tsv, resolution_sources_tsv, prov_json
  )
  missing_resolver_outputs <- required_resolver_outputs[!file.exists(required_resolver_outputs)]
  if (length(missing_resolver_outputs) > 0L) {
    stop(sprintf("Taxonomy resolver omitted expected output(s): %s",
                 paste(missing_resolver_outputs, collapse = ", ")), call. = FALSE)
  }
  all_outputs <- c(all_outputs, required_resolver_outputs)

  # Load the run-local cache so assignment-derived TaxIDs are available without
  # mutating the configured source cache in cache_only mode.
  tax_cache <- jsonlite::fromJSON(resolved_cache, simplifyVector = FALSE)
  source_df <- read.delim(resolution_sources_tsv, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(source_df), c("TaxonPath", "TaxID", "ResolutionSource")) ||
      anyDuplicated(source_df$TaxonPath)) {
    stop("Invalid taxonomy resolution-source output from resolver.", call. = FALSE)
  }
  resolution_sources <- setNames(source_df$ResolutionSource, source_df$TaxonPath)

  # Generate .kreport for each sample
  samples <- context$samples
  unclass_idx <- context$unclass_index
  count_matrix <- context$count_matrix
  lineages <- context$taxonomy$TaxonPath

  resolution_rows <- list()

  for (s in samples) {
    counts_s <- count_matrix[, s]
    total_reads <- as.numeric(sum(counts_s))
    uncl_reads <- as.numeric(counts_s[unclass_idx])

    # Build DFS abundance-sorted tree
    nodes_sorted <- build_kreport_tree(lineages, counts_s)

    # Validate tree invariants
    validate_kreport_tree(nodes_sorted, total_reads, uncl_reads)

    # Format 6-column lines
    kreport_lines <- format_kreport_lines(nodes_sorted, total_reads, uncl_reads, tax_cache)

    # Write .kreport file
    out_file <- file.path(kreport_dir, sprintf("%s.kreport", sanitize_filename(s)))
    writeLines(kreport_lines, out_file)
    all_outputs <- c(all_outputs, out_file)

    # Collect resolution info
    for (i in seq_len(nrow(nodes_sorted))) {
      p <- nodes_sorted$path[i]
      tid <- tax_cache[[p]]
      if (is.null(tid)) tid <- 0L

      resolution_rows[[length(resolution_rows) + 1L]] <- data.frame(
        SampleID = s,
        Depth = nodes_sorted$depth[i],
        RankCode = nodes_sorted$rank_code[i],
        NodeName = nodes_sorted$name[i],
        TaxonPath = p,
        TaxID = as.integer(tid),
        Status = if (tid > 0) "Resolved" else "Unresolved",
        ResolutionSource = unname(resolution_sources[[p]] %||% "unresolved"),
        stringsAsFactors = FALSE
      )
    }
  }

  # Export taxonomy resolution summary
  res_summary_file <- file.path(kreport_dir, "taxonomy_resolution.tsv")
  if (length(resolution_rows) > 0) {
    res_df <- do.call(rbind, resolution_rows)
    write.table(res_df, res_summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, res_summary_file)
  }

  list(
    status = "completed",
    outputs = all_outputs
  )
}
