# =============================================================================
# Side-effect-free preflight and immutable-input provenance
# =============================================================================

preflight_error <- function(code, message) {
  stop(sprintf("%s: %s", code, message), call. = FALSE)
}

utc_timestamp <- function(value) {
  format(as.POSIXct(value, tz = "UTC"), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
}

preflight_path_key <- function(path, os_type = .Platform$OS.type) {
  normalized <- gsub("\\\\", "/", path)
  if (identical(os_type, "windows")) tolower(normalized) else normalized
}

preflight_paths_are_same <- function(left, right, os_type = .Platform$OS.type) {
  identical(preflight_path_key(left, os_type), preflight_path_key(right, os_type))
}

preflight_path_is_same_or_descendant <- function(path, root, os_type = .Platform$OS.type) {
  path_key <- preflight_path_key(path, os_type)
  root_key <- preflight_path_key(root, os_type)
  identical(path_key, root_key) || startsWith(path_key, paste0(root_key, "/"))
}

package_location_is_project_bound <- function(location, project_library, renv_cache = NULL,
                                               resolved_project_entry = NULL,
                                               os_type = .Platform$OS.type) {
  in_project_tree <- preflight_path_is_same_or_descendant(
    location, project_library, os_type
  )
  via_project_cache_entry <- !is.null(renv_cache) &&
    preflight_path_is_same_or_descendant(location, renv_cache, os_type) &&
    !is.null(resolved_project_entry) &&
    preflight_paths_are_same(location, resolved_project_entry, os_type)
  in_project_tree || via_project_cache_entry
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
    config_file = fingerprint_file(cfg$config_file, "config_file"),
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
                                     taxonomy_cache_expected_sha256 = NULL,
                                     allow_taxonomy_mtime_change = FALSE) {
  flatten <- c(
    inventory[c("config_file", "abundance_table", "params_json", "metadata", "taxonomy_cache")],
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
    if (taxonomy_record && isTRUE(allow_taxonomy_mtime_change)) {
      if (!identical(expected[c("path", "size_bytes", "sha256")], actual[c("path", "size_bytes", "sha256")])) {
        preflight_error("E_INPUT_CHANGED", sprintf("taxonomy cache changed after restoration: '%s'", expected$path))
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
  base_pkgs <- c("base", "compiler", "datasets", "graphics", "grDevices", "grid",
                 "methods", "parallel", "splines", "stats", "stats4", "tcltk",
                 "tools", "utils", "boot", "class", "cluster", "codetools",
                 "foreign", "KernSmooth", "lattice", "MASS", "Matrix", "mgcv",
                 "nlme", "nnet", "rpart", "spatial", "survival")
  norm_proj_lib <- if (!is.null(project_library)) normalizePath(project_library, winslash = "/", mustWork = FALSE) else NULL
  renv_cache <- tryCatch(renv::paths$cache(), error = function(e) NULL)
  norm_cache <- if (!is.null(renv_cache)) normalizePath(renv_cache, winslash = "/", mustWork = FALSE) else NULL
  library_diff <- character(0)
  if (is.null(project_library) || !any(vapply(
    normalizePath(library_paths, winslash = "/", mustWork = FALSE),
    preflight_paths_are_same, logical(1), right = norm_proj_lib
  ))) {
    library_diff <- c(library_diff, sprintf("active .libPaths() does not include project library '%s'",
                                            project_library %||% "<unknown>"))
  }

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
    is_base <- package %in% base_pkgs || identical(lock_record$Priority, "recommended") || identical(lock_record$Priority, "base")
    if (!is_base && !is.null(location) && !is.null(norm_proj_lib)) {
      norm_loc <- normalizePath(location, winslash = "/", mustWork = FALSE)
      project_entry <- file.path(project_library, package)
      resolved_project_entry <- if (dir.exists(project_entry)) {
        normalizePath(project_entry, winslash = "/", mustWork = TRUE)
      } else {
        NULL
      }
      location_is_bound <- package_location_is_project_bound(
        norm_loc, norm_proj_lib, norm_cache, resolved_project_entry
      )
      if (!location_is_bound) {
        library_diff <- c(library_diff, sprintf("%s: package loaded from outside project library: '%s'",
                                                package, location))
      }
    }
    description_path <- if (!is.null(location)) file.path(location, "DESCRIPTION") else NULL
    description_sha256 <- if (!is.null(description_path) && file.exists(description_path) &&
                              !dir.exists(description_path)) {
      compute_file_hash(description_path)
    } else {
      NULL
    }
    # These fields record lock metadata only.  Version/library checks above do not
    # prove that installed package contents equal the lockfile source artifact.
    lock_identity <- lapply(c("Source", "Repository", "RemoteType", "RemoteHost", "RemoteRepo",
                              "RemoteRef", "RemoteSha", "Hash"), function(field) {
      value <- lock_record[[field]] %||% NULL
      if (is.null(value) || length(value) != 1L || is.na(value)) NULL else as.character(value)
    })
    names(lock_identity) <- c("lock_source", "lock_repository", "lock_remote_type", "lock_remote_host",
                              "lock_remote_repo", "lock_remote_ref", "lock_remote_sha", "lock_hash")
    package_locations[[index]] <- c(
      list(package = package, version = actual, path = location,
           description_sha256 = description_sha256),
      lock_identity
    )
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
  basename(path) %in% c("VERSION", "renv.lock", ".Rprofile", "activate.R", "settings.json") |
    grepl("[.](R|r|py|json|ya?ml)$", path, perl = TRUE)
}

is_krona_vendor_file <- function(path, repo_root) {
  vendor_root <- normalizePath(file.path(repo_root, "analysis", "vendor", "krona-2.8.1"),
                               winslash = "/", mustWork = FALSE)
  candidate <- normalizePath(path, winslash = "/", mustWork = FALSE)
  preflight_path_is_same_or_descendant(candidate, vendor_root)
}

maintained_source_files <- function(repo_root) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  git_scoped_paths <- c("analysis", "config.example.yml", "VERSION", "renv.lock",
                        ".Rprofile", "renv/activate.R", "renv/settings.json")
  tracked <- tryCatch({
    result <- processx::run("git", c("-c", paste0("safe.directory=", root), "-C", root,
      "ls-files", "--", git_scoped_paths),
      error_on_status = FALSE)
    if (result$status != 0L) character(0) else {
      relative <- strsplit(trimws(result$stdout), "\\r?\\n", perl = TRUE)[[1]]
      file.path(root, relative[nzchar(relative)])
    }
  }, error = function(e) character(0))
  candidates <- if (length(tracked)) tracked else c(
    list.files(file.path(root, "analysis"), recursive = TRUE, full.names = TRUE),
    file.path(root, c("config.example.yml", "VERSION", "renv.lock",
                      ".Rprofile", "renv/activate.R", "renv/settings.json"))
  )
  candidates <- candidates[file.exists(candidates) & !dir.exists(candidates) &
    (source_file_allowed(candidates) | vapply(candidates, is_krona_vendor_file,
                                               logical(1), repo_root = root))]
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
    git_scoped_paths <- c("analysis", "config.example.yml", "VERSION", "renv.lock",
                          ".Rprofile", "renv/activate.R", "renv/settings.json")
    commit <- processx::run("git", c("-c", paste0("safe.directory=", repo_root), "-C", repo_root,
      "rev-parse", "HEAD"), error_on_status = FALSE)
    if (commit$status == 0L) git_commit <- trimws(commit$stdout)
    dirty <- processx::run("git", c("-c", paste0("safe.directory=", repo_root), "-C", repo_root,
      "diff", "--quiet", "--", git_scoped_paths),
      error_on_status = FALSE)
    status <- processx::run("git", c("-c", paste0("safe.directory=", repo_root), "-C", repo_root,
      "status", "--porcelain", "--untracked-files=all", "--", git_scoped_paths),
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
  warnings <- character(0)
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
    compile <- processx::run(python, c("-c", "import ast, sys; p = sys.argv[1]; ast.parse(open(p, 'rb').read(), filename=p)", script), error_on_status = FALSE)
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
    if (probe$status != 0L) {
      err_text <- trimws(paste(probe$stderr, probe$stdout))
      code <- if (grepl("E_ONLINE_PREFLIGHT_REQUIRED", err_text)) "E_ONLINE_PREFLIGHT_REQUIRED" else "E_KREPORT_PREFLIGHT"
      preflight_error(code, err_text)
    }
  }
  if (isTRUE(cfg$krona$enabled) && isTRUE(cfg$krona$render_html)) {
    krona_exe <- cfg$krona$executable %||% "ktImportText"
    if (nzchar(krona_exe) && dir.exists(krona_exe)) {
      preflight_error("E_KRONA_PREFLIGHT", sprintf("Krona executable '%s' is a directory, not a file.", krona_exe))
    }
    renderer <- tryCatch(
      resolve_krona_renderer(cfg$krona, cfg$pipeline_root),
      error = function(e) preflight_error("E_KRONA_PREFLIGHT", e$message)
    )
    if (identical(renderer$provider, "builtin")) {
      python <- tryCatch(find_python(), error = function(e) {
        preflight_error("E_KRONA_PREFLIGHT", e$message)
      })
      builder <- renderer$builder
      compile <- processx::run(
        python,
        c("-c", "import ast, sys; p = sys.argv[1]; ast.parse(open(p, 'rb').read(), filename=p)", builder),
        error_on_status = FALSE
      )
      if (!identical(compile$status, 0L)) {
        preflight_error("E_KRONA_PREFLIGHT", trimws(paste(compile$stderr, compile$stdout)))
      }
      validate <- processx::run(
        python,
        c(builder, "--validate-only", "--vendor-dir", renderer$vendor_dir),
        error_on_status = FALSE
      )
      if (!identical(validate$status, 0L)) {
        preflight_error("E_KRONA_PREFLIGHT", trimws(paste(validate$stderr, validate$stdout)))
      }
    }
  }
  invisible(warnings)
}
