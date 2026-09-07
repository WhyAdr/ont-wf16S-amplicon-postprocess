# =============================================================================
# Atomic output helpers and whole-run publication
# =============================================================================

atomic_replace <- function(path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
  writer(temp)
  if (!file.exists(temp) || dir.exists(temp)) stop(sprintf("Atomic writer did not create '%s'.", path), call. = FALSE)
  if (file.exists(path) && !unlink(path, force = TRUE)) stop(sprintf("Could not replace '%s'.", path), call. = FALSE)
  if (!file.rename(temp, path)) stop(sprintf("Could not atomically publish '%s'.", path), call. = FALSE)
  invisible(path)
}

atomic_write_lines <- function(text, path) atomic_replace(path, function(temp) writeLines(text, temp))
atomic_write_yaml <- function(value, path) atomic_replace(path, function(temp) yaml::write_yaml(value, temp))
atomic_write_json <- function(value, path, auto_unbox = TRUE) {
  atomic_replace(path, function(temp) jsonlite::write_json(value, temp, pretty = TRUE,
    auto_unbox = auto_unbox, null = "null"))
}

prepare_run_staging <- function(final_root) {
  parent <- dirname(final_root)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(pattern = paste0(".", basename(final_root), ".staging-"), tmpdir = parent)
  if (!dir.create(stage, recursive = FALSE, showWarnings = FALSE)) {
    stop(sprintf("Could not create run staging directory '%s'.", stage), call. = FALSE)
  }
  normalizePath(stage, winslash = "/", mustWork = TRUE)
}

snapshot_staging_directory <- function(stage) {
  parent <- dirname(stage)
  snapshot <- tempfile(pattern = paste0(".", basename(stage), ".snapshot-"), tmpdir = parent)
  if (!dir.create(snapshot, recursive = FALSE, showWarnings = FALSE)) {
    stop(sprintf("Could not create module staging snapshot '%s'.", snapshot), call. = FALSE)
  }
  entries <- list.files(stage, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                        full.names = FALSE)
  directories <- entries[dir.exists(file.path(stage, entries))]
  for (relative in directories) {
    dir.create(file.path(snapshot, relative), recursive = TRUE, showWarnings = FALSE)
  }
  files <- setdiff(entries, directories)
  for (relative in files) {
    target <- file.path(snapshot, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(stage, relative), target, overwrite = FALSE)) {
      unlink(snapshot, recursive = TRUE, force = TRUE)
      stop(sprintf("Could not snapshot staged path '%s'.", relative), call. = FALSE)
    }
  }
  normalizePath(snapshot, winslash = "/", mustWork = TRUE)
}

restore_staging_directory <- function(stage, snapshot) {
  if (!dir.exists(snapshot)) stop("Module staging snapshot is missing.", call. = FALSE)
  if (dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
  if (!dir.create(stage, recursive = FALSE, showWarnings = FALSE)) {
    stop(sprintf("Could not recreate module staging directory '%s'.", stage), call. = FALSE)
  }
  entries <- list.files(snapshot, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                        full.names = FALSE)
  directories <- entries[dir.exists(file.path(snapshot, entries))]
  for (relative in directories) {
    dir.create(file.path(stage, relative), recursive = TRUE, showWarnings = FALSE)
  }
  files <- setdiff(entries, directories)
  for (relative in files) {
    target <- file.path(stage, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(snapshot, relative), target, overwrite = FALSE)) {
      stop(sprintf("Could not restore staged path '%s'.", relative), call. = FALSE)
    }
  }
  invisible(stage)
}

validate_prior_output <- function(final_root, overwrite) {
  if (!dir.exists(final_root)) return(NULL)
  entries <- list.files(final_root, all.files = TRUE, no.. = TRUE)
  if (!length(entries)) return(NULL)
  manifest_path <- file.path(final_root, "run_manifest.json")
  if (!isTRUE(overwrite)) {
    stop(sprintf("Output directory '%s' is non-empty; pass --overwrite for an owned prior run.", final_root), call. = FALSE)
  }
  if (!file.exists(manifest_path)) {
    stop(sprintf("E_OUTPUT_UNOWNED: non-empty output directory lacks a valid prior manifest: '%s'", final_root), call. = FALSE)
  }
  prior <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE), error = function(e) NULL)
  root <- normalizePath(final_root, winslash = "/", mustWork = FALSE)
  prior_root <- if (!is.null(prior$output_root)) {
    normalizePath(as.character(prior$output_root), winslash = "/", mustWork = FALSE)
  } else {
    NA_character_
  }
  if (is.null(prior) || !identical(prior$pipeline, "ont-wf16s-postprocess") ||
      is.na(prior_root) || !identical(tolower(prior_root), tolower(root))) {
    stop(sprintf("E_OUTPUT_UNOWNED: prior manifest does not prove ownership of '%s'", final_root), call. = FALSE)
  }
  if (!identical(prior$run_status, "completed")) {
    stop(sprintf("E_OUTPUT_UNOWNED: prior output '%s' is not a completed run.", final_root), call. = FALSE)
  }

  path_is_within <- function(path) {
    normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
    identical(tolower(normalized), tolower(root)) ||
      startsWith(tolower(normalized), paste0(tolower(root), "/"))
  }
  relative_owned_path <- function(path) {
    candidate <- if (is_absolute_path(path)) path else file.path(root, path)
    if (!path_is_within(candidate)) {
      stop(sprintf("E_OUTPUT_UNOWNED: prior manifest declares path outside output root: '%s'", path),
           call. = FALSE)
    }
    normalized <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    substring(normalized, nchar(root) + 2L)
  }
  prior_files <- list.files(final_root, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                            full.names = FALSE)
  prior_files <- prior_files[!dir.exists(file.path(final_root, prior_files))]
  current_contract <- identical(prior$schema_version, 2L) &&
    identical(prior$schema_revision, 1L) && !is.null(prior$owned_outputs) &&
    is.list(prior$environment) && !is.null(prior$environment$library_paths) &&
    !is.null(prior$environment$package_locations)
  if (current_contract) {
    tryCatch(validate_manifest_v2(prior, physical_root = final_root),
             error = function(e) stop(sprintf("E_OUTPUT_UNOWNED: invalid prior manifest: %s", e$message),
                                       call. = FALSE))
    owned <- as.character(unlist(prior$owned_outputs, use.names = FALSE))
    owned <- unique(vapply(owned, relative_owned_path, character(1)))
  } else {
    # v0.4.0 and earlier manifests had no ownership inventory. Derive one only
    # when every physical file is explicitly listed by a legacy module record.
    if (is.null(prior$modules) || !is.list(prior$modules)) {
      stop("E_OUTPUT_MIGRATION_REQUIRED: prior manifest predates ownership tracking and cannot be migrated safely.",
           call. = FALSE)
    }
    declared <- c(unlist(prior$owned_outputs %||% list(), use.names = FALSE),
                  unlist(lapply(prior$modules, function(record) record$outputs %||% character(0)),
                         use.names = FALSE))
    owned <- unique(c("run_manifest.json", "resolved_config.yml", "session_info.txt",
                      vapply(as.character(declared), relative_owned_path, character(1))))
    missing_declarations <- setdiff(prior_files, owned)
    if (length(missing_declarations)) {
      stop(sprintf("E_OUTPUT_MIGRATION_REQUIRED: legacy output contains undeclared file(s): %s",
                   paste(missing_declarations, collapse = ", ")), call. = FALSE)
    }
    prior$ownership_migrated <- TRUE
  }
  if (length(setdiff(owned, prior_files))) {
    stop(sprintf("E_OUTPUT_UNOWNED: prior manifest declares missing file(s): %s",
                 paste(setdiff(owned, prior_files), collapse = ", ")), call. = FALSE)
  }
  prior$owned_outputs <- if (exists("json_array", mode = "function")) json_array(owned) else owned
  prior
}

preserve_unowned_outputs <- function(final_root, stage, prior_manifest) {
  if (is.null(prior_manifest) || !dir.exists(final_root)) return(invisible(TRUE))
  owned <- unique(c("run_manifest.json", "resolved_config.yml", "session_info.txt",
                    unlist(prior_manifest$owned_outputs %||% list(), use.names = FALSE)))
  files <- list.files(final_root, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  files <- files[!dir.exists(file.path(final_root, files))]
  unowned <- setdiff(files, owned)
  for (relative in unowned) {
    target <- file.path(stage, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(final_root, relative), target, overwrite = FALSE)) {
      stop(sprintf("Could not preserve unowned output '%s'.", relative), call. = FALSE)
    }
  }
  invisible(TRUE)
}

publish_staged_run <- function(stage, final_root) {
  parent <- dirname(final_root)
  backup <- tempfile(pattern = paste0(".", basename(final_root), ".previous-"), tmpdir = parent)
  had_prior <- dir.exists(final_root)
  if (had_prior && !file.rename(final_root, backup)) stop("Could not preserve the prior completed run.", call. = FALSE)
  committed <- FALSE
  on.exit({
    if (!committed && had_prior && dir.exists(backup) && !dir.exists(final_root)) file.rename(backup, final_root)
  }, add = TRUE)
  if (!file.rename(stage, final_root)) stop("Could not publish staged run.", call. = FALSE)
  committed <- TRUE
  if (had_prior && dir.exists(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  invisible(final_root)
}
