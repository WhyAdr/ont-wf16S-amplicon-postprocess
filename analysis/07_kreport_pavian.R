# =============================================================================
# Module 07: Kraken Report (.kreport) Generation for Pavian Sankey
# =============================================================================

suppressMessages({
  library(jsonlite)
  library(dplyr)
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
  prov_json <- file.path(kreport_dir, "taxonomy_provenance.json")
  
  # Invoke Python resolver if cache does not exist or refresh is requested
  need_resolver <- !file.exists(cache_file) || (network_mode == "refresh")
  
  py_script <- if (!is.null(context$config$config_dir) && nzchar(context$config$config_dir)) {
    file.path(context$config$config_dir, "analysis", "utils", "ncbi_taxonomy.py")
  } else {
    file.path("analysis", "utils", "ncbi_taxonomy.py")
  }
  
  if (!file.exists(py_script)) {
    candidates <- c(
      file.path("analysis", "utils", "ncbi_taxonomy.py"),
      file.path("..", "analysis", "utils", "ncbi_taxonomy.py"),
      file.path("..", "..", "analysis", "utils", "ncbi_taxonomy.py")
    )
    for (cand in candidates) {
      if (file.exists(cand)) {
        py_script <- cand
        break
      }
    }
  }
  
  first_sample <- context$samples[1]
  first_asgn <- if (!is.null(context$assignments)) context$assignments[[first_sample]] else NULL
  
  cmd_args <- c(
    py_script,
    "--abundance", cfg$input$abundance_table,
    "--cache", cache_file,
    "--mode", network_mode,
    "--unresolved-policy", unresolved_policy,
    "--unresolved-tsv", unresolved_tsv,
    "--provenance", prov_json
  )
  if (!is.null(first_asgn) && file.exists(first_asgn)) {
    cmd_args <- c(cmd_args, "--assignments", first_asgn)
  }
  
  # Run Python script
  res_code <- system2("python", args = cmd_args)
  if (res_code != 0) {
    stop(sprintf("NCBI taxonomy resolver failed with exit status %d", res_code), call. = FALSE)
  }
  
  if (file.exists(unresolved_tsv)) all_outputs <- c(all_outputs, unresolved_tsv)
  if (file.exists(prov_json)) all_outputs <- c(all_outputs, prov_json)
  
  # Load taxonomy cache
  tax_cache <- if (file.exists(cache_file)) {
    jsonlite::fromJSON(cache_file)
  } else {
    list()
  }
  
  # Generate .kreport for each sample
  samples <- context$samples
  unclass_idx <- context$unclass_index
  count_matrix <- context$count_matrix
  lineages <- context$taxonomy$TaxonPath
  
  resolution_rows <- list()
  
  for (s in samples) {
    counts_s <- count_matrix[, s]
    total_reads <- as.integer(sum(counts_s))
    uncl_reads <- as.integer(counts_s[unclass_idx])
    
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
      
      resolution_rows[[p]] <- data.frame(
        SampleID = s,
        Depth = nodes_sorted$depth[i],
        RankCode = nodes_sorted$rank_code[i],
        NodeName = nodes_sorted$name[i],
        TaxonPath = p,
        TaxID = as.integer(tid),
        Status = if (tid > 0) "Resolved" else "Unresolved",
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
