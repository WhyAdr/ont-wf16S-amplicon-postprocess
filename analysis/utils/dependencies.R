# =============================================================================
# Dependency Inspector Utility
# =============================================================================

REQUIRED_PACKAGES <- c(
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
  "testthat"
)

get_required_packages <- function() {
  REQUIRED_PACKAGES
}

check_dependencies <- function(pkgs = REQUIRED_PACKAGES) {
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

get_dependency_versions <- function(pkgs = REQUIRED_PACKAGES) {
  installed <- rownames(installed.packages())
  vapply(pkgs, function(pkg) {
    if (pkg %in% installed) as.character(packageVersion(pkg)) else NA_character_
  }, FUN.VALUE = character(1))
}
