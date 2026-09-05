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
  "RColorBrewer",
  "jsonlite",
  "pheatmap",
  "UpSetR",
  "ggrepel",
  "digest",
  "processx"
)

TEST_PACKAGES <- c("testthat")
REQUIRED_PACKAGES <- unique(c(RUNTIME_PACKAGES, TEST_PACKAGES))

get_required_packages <- function(include_tests = FALSE) {
  if (isTRUE(include_tests)) REQUIRED_PACKAGES else RUNTIME_PACKAGES
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
