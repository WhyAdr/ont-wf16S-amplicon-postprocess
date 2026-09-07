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
  if (is.null(prior) || !identical(prior$pipeline, "ont-wf16s-postprocess") ||
      !identical(prior$output_root, normalizePath(final_root, winslash = "/", mustWork = FALSE))) {
    stop(sprintf("E_OUTPUT_UNOWNED: prior manifest does not prove ownership of '%s'", final_root), call. = FALSE)
  }
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
