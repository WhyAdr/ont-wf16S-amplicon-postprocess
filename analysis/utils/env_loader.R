# =============================================================================
# Self-sufficient environment loader for external CWD execution (P2-1)
# =============================================================================

find_pipeline_repo_root <- function(file_arg = NULL) {
  if (is.null(file_arg)) {
    args <- commandArgs(trailingOnly = FALSE)
    match_arg <- grep("^--file=", args, value = TRUE)
    file_arg <- if (length(match_arg)) sub("^--file=", "", match_arg[1]) else NULL
  }
  if (!is.null(file_arg) && nzchar(file_arg)) {
    norm_file <- normalizePath(file_arg, winslash = "/", mustWork = FALSE)
    script_dir <- dirname(norm_file)
    if (basename(script_dir) == "analysis") {
      return(normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE))
    }
    if (file.exists(file.path(script_dir, "VERSION"))) {
      return(script_dir)
    }
    if (file.exists(file.path(dirname(script_dir), "VERSION"))) {
      return(normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

load_pipeline_environment <- function(repo_root = NULL, fatal_fn = stop) {
  if (is.null(repo_root)) repo_root <- find_pipeline_repo_root()
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)

  activate_script <- file.path(repo_root, "renv", "activate.R")
  if (!file.exists(activate_script)) {
    fatal_fn("Environment activation", simpleError(sprintf(
      "E_RENV_NOT_RESTORED: renv bootstrap script '%s' is missing.", activate_script
    )))
  }

  tryCatch({
    source(activate_script, local = FALSE)
  }, error = function(e) {
    fatal_fn("Environment activation", simpleError(sprintf(
      "E_RENV_NOT_RESTORED: could not bootstrap renv: %s", conditionMessage(e)
    )))
  })

  if (!requireNamespace("renv", quietly = TRUE)) {
    fatal_fn("Environment activation", simpleError(
      "E_RENV_NOT_RESTORED: package 'renv' is not available after bootstrapping."
    ))
  }

  tryCatch({
    renv::load(project = repo_root)
  }, error = function(e) {
    fatal_fn("Environment activation", simpleError(sprintf(
      "E_RENV_NOT_RESTORED: failed to load project environment at '%s': %s",
      repo_root, conditionMessage(e)
    )))
  })

  project_library <- tryCatch(
    normalizePath(renv::paths$library(project = repo_root), winslash = "/", mustWork = FALSE),
    error = function(e) NA_character_
  )
  lib_desc <- if (is.na(project_library) || is.null(project_library)) "<unknown>" else project_library
  if (is.na(project_library) || !dir.exists(project_library)) {
    fatal_fn("Environment activation", simpleError(sprintf(
      "E_RENV_NOT_RESTORED: project library '%s' does not exist; run 'Rscript analysis/install_packages.R --restore'",
      lib_desc
    )))
  }

  .libPaths(unique(c(project_library, .libPaths())))
  invisible(repo_root)
}
