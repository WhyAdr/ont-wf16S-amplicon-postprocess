# =============================================================================
# Module 07: Kraken Report (.kreport) and optional Krona export
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
  if (identical(network_mode, "refresh")) {
    cmd_args <- c(cmd_args, "--defer-cache-commit", "--cache-lock-held",
                  "--cache-lock-owner-pid", as.character(Sys.getpid()),
                  "--transaction-id", context$transaction_id)
  }
  assignment_paths <- unname(unlist(context$assignments, use.names = FALSE))
  if (length(assignment_paths) > 0L) {
    cmd_args <- c(cmd_args, as.vector(rbind("--assignments", assignment_paths)))
  }
  expected_records <- c(
    list(context$input_inventory$abundance_table, context$input_inventory$taxonomy_cache),
    context$input_inventory$assignments
  )
  expected_records <- Filter(function(record) !is.null(record) && !is.null(record$path) &&
    !is.null(record$sha256), expected_records)
  if (length(expected_records)) {
    expected_specs <- vapply(expected_records, function(record) paste(record$path, record$sha256, sep = "\t"),
                             character(1))
    cmd_args <- c(cmd_args, as.vector(rbind("--expected-input", expected_specs)))
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
  conflict_df <- read.delim(conflicts_tsv, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(conflict_df), c("Lineage", "WinnerTaxID", "TaxIDCountsJSON")) ||
      anyDuplicated(conflict_df$Lineage)) {
    stop("Invalid taxonomy conflict output from resolver.", call. = FALSE)
  }
  if (nrow(conflict_df) > 0L) {
    warning(sprintf(
      "%d lineage-to-TaxID conflict(s) used the documented modal-count/minimum-TaxID tie-break; review taxonomy_conflicts.tsv.",
      nrow(conflict_df)), call. = FALSE)
  }

  krona_cfg <- cfg$krona %||% list(
    enabled = FALSE,
    render_html = FALSE,
    html_renderer = "builtin",
    executable = "ktImportText"
  )
  krona_enabled <- isTRUE(krona_cfg$enabled)
  render_html <- krona_enabled && isTRUE(krona_cfg$render_html)
  krona_dir <- file.path(kreport_dir, "krona")
  krona_provenance_file <- file.path(krona_dir, "krona_provenance.json")
  krona_executable <- NA_character_
  krona_requested_executable <- krona_cfg$executable %||% "ktImportText"
  krona_renderer <- NULL
  krona_vendor <- NULL
  html_status <- if (render_html) "pending" else "not_requested"
  krona_records <- list()

  if (krona_enabled) {
    dir.create(krona_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(krona_dir)) {
      stop(sprintf("Could not create Krona output directory '%s'.", krona_dir), call. = FALSE)
    }

    if (render_html) {
      krona_renderer <- resolve_krona_renderer(krona_cfg, cfg$pipeline_root)
      krona_executable <- krona_renderer$resolved_executable %||% NA_character_
      if (identical(krona_renderer$provider, "builtin")) {
        krona_vendor <- krona_vendor_manifest(krona_renderer$vendor_dir)
      }
    }
  }

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

    if (krona_enabled) {
      sample_filename <- sanitize_filename(s)
      krona_tsv <- file.path(krona_dir, sprintf("%s.krona.tsv", sample_filename))
      write_krona_input(krona_tsv, nodes_sorted, total_reads, uncl_reads)
      all_outputs <- c(all_outputs, krona_tsv)

      emitted_magnitudes <- nodes_sorted$reads_taxon[nodes_sorted$reads_taxon > 0]
      emitted_sum <- uncl_reads + sum(emitted_magnitudes)
      sample_record <- list(
        sample_id = s,
        total_reads = as.numeric(total_reads),
        classified_reads = as.numeric(total_reads - uncl_reads),
        unclassified_reads = as.numeric(uncl_reads),
        emitted_magnitude_sum = as.numeric(emitted_sum),
        tsv_path = krona_artifact_relpath(krona_tsv, cfg$output$base_dir),
        tsv_sha256 = compute_file_hash(krona_tsv)
      )

      if (render_html) {
        krona_html <- file.path(krona_dir, sprintf("%s.krona.html", sample_filename))
        if (identical(krona_renderer$provider, "builtin")) {
          render_builtin_krona_html(
            python_cmd = python_cmd,
            builder_path = krona_renderer$builder,
            vendor_dir = krona_renderer$vendor_dir,
            output_path = krona_html,
            sample_id = s,
            input_path = krona_tsv,
            expected_total = total_reads
          )
        } else {
          render_kronatools_html(krona_executable, krona_html, s, krona_tsv)
        }
        all_outputs <- c(all_outputs, krona_html)
        sample_record$html_path <- krona_artifact_relpath(krona_html, cfg$output$base_dir)
        sample_record$html_sha256 <- compute_file_hash(krona_html)
      }
      krona_records[[length(krona_records) + 1L]] <- sample_record
    }

    # Collect resolution info
    for (i in seq_len(nrow(nodes_sorted))) {
      p <- nodes_sorted$path[i]
      tid <- as.character(tax_cache[[p]] %||% "0")
      if (!grepl("^(0|[1-9][0-9]{0,15})$", tid) ||
          (nchar(tid) == 16L && tid > "9007199254740991")) {
        stop(sprintf("Resolver returned a noncanonical TaxID for '%s'.", p), call. = FALSE)
      }

      resolution_source <- unname(resolution_sources[[p]] %||% "unresolved")
      resolution_rows[[length(resolution_rows) + 1L]] <- data.frame(
        SampleID = s,
        Depth = nodes_sorted$depth[i],
        RankCode = nodes_sorted$rank_code[i],
        NodeName = nodes_sorted$name[i],
        TaxonPath = p,
        TaxID = tid,
        Status = if (identical(resolution_source, "assignment_conflict")) {
          "Conflicted"
        } else if (tid != "0") {
          "Resolved"
        } else {
          "Unresolved"
        },
        ResolutionSource = resolution_source,
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

  if (krona_enabled) {
    html_rendered <- isTRUE(render_html)
    if (html_rendered) {
      html_status <- "rendered"
    }
    renderer_policy <- tolower(as.character(krona_cfg$html_renderer %||% "builtin"))
    renderer_manifest <- if (identical(krona_renderer$provider %||% "", "builtin")) {
      krona_vendor
    } else {
      NULL
    }
    krona_provenance <- list(
      schema_version = 1L,
      path_basis = "run_dir",
      format = "Krona 2.0-compatible direct-count tab-delimited lineage format",
      renderer_policy = renderer_policy,
      renderer = if (html_rendered) krona_renderer$renderer else NULL,
      renderer_version = if (html_rendered) krona_renderer$renderer_version else NULL,
      krona_version = if (html_rendered) krona_renderer$krona_version else NULL,
      vendor_source = if (!is.null(renderer_manifest)) renderer_manifest$source_url else NULL,
      vendor_manifest_sha256 = if (!is.null(renderer_manifest)) {
        renderer_manifest$manifest_sha256
      } else {
        NULL
      },
      vendor_sha256_manifest = if (!is.null(renderer_manifest)) {
        json_array(renderer_manifest$entries)
      } else {
        json_array(character(0))
      },
      html_status = html_status,
      standalone_html = isTRUE(html_rendered),
      count_model = "direct abundance-table taxon counts plus canonical unclassified count",
      denominator = "TotalReads",
      classified_definition = "sum of direct positive-count classified taxonomy rows",
      requested_executable = krona_requested_executable,
      resolved_executable = if (html_rendered) {
        krona_renderer$resolved_executable %||% krona_renderer$builder
      } else {
        NULL
      },
      render_html = render_html,
      samples = json_array(krona_records)
    )
    jsonlite::write_json(
      krona_provenance,
      krona_provenance_file,
      pretty = TRUE,
      auto_unbox = TRUE,
      null = "null"
    )
    if (!file.exists(krona_provenance_file) ||
        !isTRUE(file.info(krona_provenance_file)$size > 0)) {
      stop(sprintf("Krona provenance was not written: '%s'.", krona_provenance_file),
           call. = FALSE)
    }
    all_outputs <- c(all_outputs, krona_provenance_file)
  }

  list(
    status = "completed",
    outputs = all_outputs
  )
}
