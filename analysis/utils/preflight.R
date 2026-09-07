# =============================================================================
# Side-effect-free preflight and immutable-input provenance
# =============================================================================

preflight_error <- function(code, message) {
  stop(sprintf("%s: %s", code, message), call. = FALSE)
}

utc_timestamp <- function(value) {
  format(as.POSIXct(value, tz = "UTC"), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
}

canonical_existing_path <- function(path, label) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !file.exists(path) || dir.exists(path)) {
    preflight_error("E_INPUT_MISSING", sprintf("%s is not an existing regular file: '%s'", label, path))
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

fingerprint_file <- function(path, label = "input") {
  canonical <- canonical_existing_path(path, label)
  info <- file.info(canonical)
  list(
    path = canonical,
    size_bytes = as.numeric(info$size),
    mtime_utc = utc_timestamp(info$mtime),
    sha256 = compute_file_hash(canonical)
  )
}

inventory_inputs <- function(cfg, bamstats = NULL) {
  records <- list(
    abundance_table = fingerprint_file(cfg$input$abundance_table, "input.abundance_table"),
    params_json = fingerprint_file(cfg$input$params_json, "input.params_json")
  )
  if (!is.null(cfg$input$metadata)) {
    records$metadata <- fingerprint_file(cfg$input$metadata, "input.metadata")
  }
  if (!is.null(cfg$taxonomy$cache)) {
    records$taxonomy_cache <- fingerprint_file(cfg$taxonomy$cache, "taxonomy.cache")
  }
  assignments <- cfg$input$assignments %||% list()
  assignment_records <- lapply(names(assignments), function(sample_id) {
    c(list(sample_id = sample_id), fingerprint_file(assignments[[sample_id]],
      sprintf("input.assignments.%s", sample_id)))
  })
  canonical_assignments <- vapply(assignment_records, `[[`, character(1), "path")
  if (anyDuplicated(tolower(canonical_assignments))) {
    preflight_error("E_DUPLICATE_INPUT", "one physical assignment file is configured for multiple samples")
  }
  records$assignments <- assignment_records

  bamstats <- bamstats %||% character(0)
  present <- names(bamstats)[!is.na(bamstats)]
  records$bamstats <- lapply(present, function(sample_id) {
    c(list(sample_id = sample_id), fingerprint_file(bamstats[[sample_id]],
      sprintf("bamstats.%s", sample_id)))
  })
  records
}

assert_inputs_unchanged <- function(inventory, allow_taxonomy_cache_change = FALSE,
                                     taxonomy_cache_expected_sha256 = NULL) {
  flatten <- c(
    inventory[c("abundance_table", "params_json", "metadata", "taxonomy_cache")],
    inventory$assignments,
    inventory$bamstats
  )
  for (expected in Filter(Negate(is.null), flatten)) {
    actual <- fingerprint_file(expected$path, expected$path)
    fields <- c("path", "size_bytes", "mtime_utc", "sha256")
    taxonomy_record <- !is.null(inventory$taxonomy_cache) &&
      identical(expected$path, inventory$taxonomy_cache$path)
    if (taxonomy_record && isTRUE(allow_taxonomy_cache_change)) {
      expected_hash <- taxonomy_cache_expected_sha256 %||% expected$sha256
      if (!identical(actual$sha256, expected_hash)) {
        preflight_error("E_INPUT_CHANGED", sprintf("taxonomy cache did not match the expected transition: '%s'",
                                                     expected$path))
      }
      next
    }
    if (!identical(expected[fields], actual[fields])) {
      preflight_error("E_INPUT_CHANGED", sprintf("input changed after inventory: '%s'", expected$path))
    }
  }
  invisible(TRUE)
}

read_lock_status <- function(repo_root, packages) {
  lockfile <- file.path(repo_root, "renv.lock")
  library_paths <- normalizePath(.libPaths(), winslash = "/", mustWork = FALSE)
  project_library <- tryCatch(
    normalizePath(renv::paths$library(project = repo_root), winslash = "/", mustWork = FALSE),
    error = function(e) NULL
  )
  if (length(project_library) && is.na(project_library)) project_library <- NULL
  if (!file.exists(lockfile)) {
    return(list(lock_status = "missing", locked = FALSE, lockfile = NULL,
                lockfile_sha256 = NULL, expected_r = NULL, actual_r = as.character(getRversion()),
                r_discrepancies = json_array("renv.lock is missing"),
                package_discrepancies = json_array(character(0)),
                library_discrepancies = json_array(character(0)),
                library_paths = json_array(library_paths), project_library = project_library,
                package_locations = json_array(list())))
  }
  lock <- jsonlite::fromJSON(lockfile, simplifyVector = FALSE)
  expected_r <- as.character(lock$R$Version %||% "")
  actual_r <- as.character(getRversion())
  r_diff <- if (identical(expected_r, actual_r)) character(0) else
    sprintf("expected R %s, found %s", expected_r, actual_r)
  package_diff <- character(0)
  lock_packages <- sort(unique(c(names(lock$Packages %||% list()), packages)))
  package_locations <- vector("list", length(lock_packages))
  for (index in seq_along(lock_packages)) {
    package <- lock_packages[[index]]
    lock_record <- lock$Packages[[package]]
    expected <- as.character(lock_record$Version %||% "<missing from lock>")
    location <- tryCatch(find.package(package, quiet = TRUE)[1], error = function(e) NA_character_)
    if (!length(location) || is.na(location) || !nzchar(location)) location <- NULL
    actual <- if (!is.null(location)) as.character(utils::packageVersion(package)) else "<not installed>"
    versions_match <- if (grepl("^[0-9]+([.][0-9]+|[-][0-9]+)+$", expected) &&
                          grepl("^[0-9]+([.][0-9]+|[-][0-9]+)+$", actual)) {
      tryCatch(utils::compareVersion(expected, actual) == 0L, error = function(e) FALSE)
    } else {
      identical(expected, actual)
    }
    if (!versions_match) {
      package_diff <- c(package_diff, sprintf("%s: expected %s, found %s", package, expected, actual))
    }
    package_locations[[index]] <- list(package = package, version = actual, path = location)
  }
  library_diff <- character(0)
  if (is.null(project_library) || !any(tolower(normalizePath(library_paths, winslash = "/", mustWork = FALSE)) ==
                                       tolower(project_library))) {
    library_diff <- sprintf("active .libPaths() does not include project library '%s'",
                            project_library %||% "<unknown>")
  }
  synchronized <- length(r_diff) == 0L && length(package_diff) == 0L && length(library_diff) == 0L
  list(
    lock_status = if (synchronized) "synchronized" else "mismatch",
    locked = synchronized,
    lockfile = "renv.lock",
    lockfile_sha256 = compute_file_hash(lockfile),
    expected_r = expected_r,
    actual_r = actual_r,
    r_discrepancies = json_array(r_diff),
    package_discrepancies = json_array(package_diff),
    library_discrepancies = json_array(library_diff),
    library_paths = json_array(library_paths),
    project_library = project_library,
    package_locations = json_array(package_locations)
  )
}

source_file_allowed <- function(path) {
  basename(path) %in% c("VERSION", "renv.lock") |
    grepl("[.](R|r|py|json|ya?ml)$", path, perl = TRUE)
}

maintained_source_files <- function(repo_root) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  tracked <- tryCatch({
    result <- processx::run("git", c("-c", paste0("safe.directory=", root), "-C", root,
      "ls-files", "--", "analysis", "config.example.yml", "VERSION", "renv.lock"),
      error_on_status = FALSE)
    if (result$status != 0L) character(0) else {
      relative <- strsplit(trimws(result$stdout), "\\r?\\n", perl = TRUE)[[1]]
      file.path(root, relative[nzchar(relative)])
    }
  }, error = function(e) character(0))
  candidates <- if (length(tracked)) tracked else c(
    list.files(file.path(root, "analysis"), recursive = TRUE, full.names = TRUE),
    file.path(root, c("config.example.yml", "VERSION", "renv.lock"))
  )
  candidates <- candidates[file.exists(candidates) & !dir.exists(candidates) & source_file_allowed(candidates)]
  sort(normalizePath(candidates, winslash = "/", mustWork = TRUE))
}

source_provenance <- function(repo_root) {
  files <- maintained_source_files(repo_root)
  relative <- substring(files, nchar(normalizePath(repo_root, winslash = "/")) + 2L)
  entries <- paste(relative, vapply(files, compute_file_hash, character(1)), sep = "\t")
  digest_value <- digest::digest(paste(entries, collapse = "\n"), algo = "sha256", serialize = FALSE)
  git_commit <- NULL
  git_dirty <- NULL
  try({
    commit <- processx::run("git", c("-c", paste0("safe.directory=", repo_root), "-C", repo_root,
      "rev-parse", "HEAD"), error_on_status = FALSE)
    if (commit$status == 0L) git_commit <- trimws(commit$stdout)
    dirty <- processx::run("git", c("-c", paste0("safe.directory=", repo_root), "-C", repo_root,
      "diff", "--quiet", "--", "analysis", "config.example.yml", "VERSION", "renv.lock"),
      error_on_status = FALSE)
    status <- processx::run("git", c("-c", paste0("safe.directory=", repo_root), "-C", repo_root,
      "status", "--porcelain", "--untracked-files=all", "--", "analysis", "config.example.yml", "VERSION", "renv.lock"),
      error_on_status = FALSE)
    git_dirty <- !identical(dirty$status, 0L) || nzchar(trimws(status$stdout))
  }, silent = TRUE)
  list(git_commit = git_commit, git_dirty = git_dirty, source_digest_sha256 = digest_value,
       source_files = json_array(relative))
}

validate_output_root <- function(cfg, repo_root, extra_paths = character(0)) {
  output <- normalizePath(cfg$output$base_dir, winslash = "/", mustWork = FALSE)
  homes <- c(Sys.getenv("USERPROFILE"), Sys.getenv("HOME"))
  forbidden <- unique(normalizePath(c(repo_root, homes[nzchar(homes)]),
                                    winslash = "/", mustWork = FALSE))
  volume <- normalizePath(dirname(output), winslash = "/", mustWork = FALSE)
  while (!identical(dirname(volume), volume)) volume <- dirname(volume)
  forbidden <- c(forbidden, volume)
  if (tolower(output) %in% tolower(forbidden)) {
    preflight_error("E_OUTPUT_UNSAFE", sprintf("unsafe output root '%s'", output))
  }
  inputs <- c(cfg$input$abundance_table, cfg$input$metadata, cfg$input$params_json,
              cfg$taxonomy$cache, unlist(cfg$input$assignments, use.names = FALSE), extra_paths)
  input_paths <- normalizePath(inputs[!is.na(inputs) & nzchar(inputs)], winslash = "/", mustWork = FALSE)
  overlaps <- function(left, right) {
    tolower(left) == tolower(right) ||
      startsWith(tolower(left), paste0(tolower(right), "/")) ||
      startsWith(tolower(right), paste0(tolower(left), "/"))
  }
  if (any(vapply(input_paths, overlaps, logical(1), right = output))) {
    preflight_error("E_OUTPUT_UNSAFE", "output root overlaps an input path")
  }
  upstream <- cfg$input$wf16s_output_root
  if (!is.null(upstream) && any(vapply(normalizePath(upstream, winslash = "/", mustWork = FALSE),
                                     overlaps, logical(1), right = output))) {
    preflight_error("E_OUTPUT_UNSAFE", "output root is nested within the upstream input root")
  }
  invisible(output)
}

run_module_preflight <- function(context, modules) {
  cfg <- context$config
  if ("qc" %in% modules) {
    for (sample_id in names(context$assignment_data)) {
      bamstats <- context$bamstats[[sample_id]]
      if (!is.null(bamstats) && !is.na(bamstats)) {
        partition_minimap2_failures(context$assignment_data[[sample_id]], bamstats,
                                    context$params, sample_id)
      }
    }
  }
  if ("faprotax" %in% modules) invisible(validate_faprotax_runtime())
  if ("kreport" %in% modules) {
    python <- tryCatch(find_python(), error = function(e) preflight_error("E_KREPORT_PREFLIGHT", e$message))
    script <- file.path(cfg$pipeline_root, "analysis", "utils", "ncbi_taxonomy.py")
    compile <- processx::run(python, c("-m", "py_compile", script), error_on_status = FALSE)
    if (compile$status != 0L) preflight_error("E_KREPORT_PREFLIGHT", trimws(compile$stderr))
    args <- c(script, "--validate-only", "--abundance", cfg$input$abundance_table,
      "--tax-column", cfg$input$tax_column, "--cache", cfg$taxonomy$cache,
      "--mode", cfg$taxonomy$network_mode, "--email-env", cfg$taxonomy$email_env,
      "--api-key-env", cfg$taxonomy$api_key_env, "--unresolved-policy", cfg$taxonomy$unresolved_policy)
    if (isTRUE(cfg$cli$online_preflight)) args <- c(args, "--online-preflight")
    assignments <- unname(unlist(context$assignments, use.names = FALSE))
    if (length(assignments)) args <- c(args, as.vector(rbind("--assignments", assignments)))
    expected_records <- c(
      list(context$input_inventory$abundance_table, context$input_inventory$taxonomy_cache),
      context$input_inventory$assignments
    )
    expected_records <- Filter(function(record) !is.null(record) && !is.null(record$path) &&
      !is.null(record$sha256), expected_records)
    if (length(expected_records)) {
      expected_specs <- vapply(expected_records,
        function(record) paste(record$path, record$sha256, sep = "\t"), character(1))
      args <- c(args, as.vector(rbind("--expected-input", expected_specs)))
    }
    probe <- processx::run(python, args, error_on_status = FALSE)
    if (probe$status != 0L) preflight_error("E_KREPORT_PREFLIGHT", trimws(paste(probe$stderr, probe$stdout)))
  }
  if (isTRUE(cfg$krona$enabled) && isTRUE(cfg$krona$render_html) &&
      is.na(find_krona_executable(cfg$krona$executable))) {
    preflight_error("E_KRONA_PREFLIGHT", sprintf("Krona executable '%s' was not found", cfg$krona$executable))
  }
  invisible(TRUE)
}
