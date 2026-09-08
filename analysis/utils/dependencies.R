# =============================================================================
# Dependency Inspector Utility
# =============================================================================

RUNTIME_PACKAGES <- c(
  "yaml",
  "optparse",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "scales",
  "vegan",
  "permute",
  "RColorBrewer",
  "jsonlite",
  "pheatmap",
  "UpSetR",
  "digest",
  "processx",
  "filelock"
)

MODULE_PACKAGES <- list(
  faprotax = "microeco"
)
MICROECO_MIN_VERSION <- "2.3.0"

TEST_PACKAGES <- c("testthat")
REQUIRED_PACKAGES <- unique(c(RUNTIME_PACKAGES, TEST_PACKAGES))

get_module_packages <- function(modules) {
  requested <- intersect(unique(modules), names(MODULE_PACKAGES))
  unique(unlist(MODULE_PACKAGES[requested], use.names = FALSE))
}

check_module_dependencies <- function(modules) {
  pkgs <- get_module_packages(modules)
  check_dependencies(pkgs)

  if ("faprotax" %in% modules &&
      utils::compareVersion(as.character(utils::packageVersion("microeco")),
                            MICROECO_MIN_VERSION) < 0) {
    stop(sprintf("Module 'faprotax' requires microeco >= %s.", MICROECO_MIN_VERSION),
         call. = FALSE)
  }
  invisible(TRUE)
}

get_required_packages <- function(include_tests = FALSE, include_modules = FALSE) {
  ans <- if (isTRUE(include_tests)) REQUIRED_PACKAGES else RUNTIME_PACKAGES
  if (isTRUE(include_modules)) {
    ans <- unique(c(ans, unlist(MODULE_PACKAGES, use.names = FALSE)))
  }
  ans
}

check_dependencies <- function(pkgs = RUNTIME_PACKAGES) {
  installed <- rownames(installed.packages())
  missing_pkgs <- setdiff(pkgs, installed)
  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "Missing required R packages: %s\nPlease run: Rscript analysis/install_packages.R --install",
      paste(missing_pkgs, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

get_dependency_versions <- function(pkgs = RUNTIME_PACKAGES) {
  installed <- rownames(installed.packages())
  vapply(pkgs, function(pkg) {
    if (pkg %in% installed) as.character(packageVersion(pkg)) else NA_character_
  }, FUN.VALUE = character(1))
}

find_python <- function() {
  candidates <- unname(c(Sys.which("python3"), Sys.which("python")))
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    probe <- tryCatch(
      processx::run(candidate, "--version", error_on_status = FALSE),
      error = function(e) NULL
    )
    if (!is.null(probe) && identical(probe$status, 0L)) return(candidate)
  }
  stop("Neither 'python3' nor 'python' was found on PATH.", call. = FALSE)
}

find_krona_executable <- function(executable = "ktImportText") {
  if (is.null(executable) || length(executable) != 1L ||
      is.na(executable) || !nzchar(trimws(executable))) {
    return(NA_character_)
  }

  executable <- trimws(executable)
  if (file.exists(executable) || grepl("[/\\\\]", executable)) {
    if (!file.exists(executable) || dir.exists(executable)) return(NA_character_)
    return(normalizePath(executable, winslash = "/", mustWork = TRUE))
  }

  found <- Sys.which(executable)
  if (!nzchar(found)) {
    NA_character_
  } else {
    normalizePath(found, winslash = "/", mustWork = TRUE)
  }
}

get_krona_version <- function(executable) {
  if (is.null(executable) || length(executable) != 1L || is.na(executable) ||
      !nzchar(executable)) {
    return(NULL)
  }

  probe <- tryCatch(
    processx::run(executable, args = character(0), error_on_status = FALSE),
    error = function(e) NULL
  )
  if (is.null(probe)) return(NULL)

  output <- paste(c(probe$stdout, probe$stderr), collapse = "\n")
  match <- regexec("KronaTools[[:space:]]+([0-9]+([.][0-9]+){1,3})",
                   output, perl = TRUE)
  captures <- regmatches(output, match)[[1]]
  if (length(captures) >= 2L) captures[2] else NULL
}
