#!/usr/bin/env Rscript
# =============================================================================
# Dependency Checker and Installer for ONT wf-16s Post-Processing Pipeline
#
# Usage:
#   Rscript analysis/install_packages.R            # Check only (exits non-zero if missing)
#   Rscript analysis/install_packages.R --install  # Install missing packages from CRAN
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
do_install <- "--install" %in% args

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

installed <- rownames(installed.packages())
missing_pkgs <- setdiff(REQUIRED_PACKAGES, installed)

cat("=== ONT wf-16s Pipeline Dependency Check ===\n")
cat(sprintf("Total required packages: %d\n", length(REQUIRED_PACKAGES)))
cat(sprintf("Installed: %d\n", length(REQUIRED_PACKAGES) - length(missing_pkgs)))
cat(sprintf("Missing:   %d\n\n", length(missing_pkgs)))

for (pkg in REQUIRED_PACKAGES) {
  is_inst <- pkg %in% installed
  ver <- if (is_inst) as.character(packageVersion(pkg)) else "NOT INSTALLED"
  status <- if (is_inst) "[OK]" else "[MISSING]"
  cat(sprintf("  %-15s %-10s %s\n", pkg, status, ver))
}

if (length(missing_pkgs) > 0) {
  if (do_install) {
    cat(sprintf("\nAttempting installation of %d missing packages...\n", length(missing_pkgs)))
    repos <- "https://cloud.r-project.org"
    install.packages(missing_pkgs, repos = repos)

    # Re-verify
    installed_now <- rownames(installed.packages())
    still_missing <- setdiff(REQUIRED_PACKAGES, installed_now)
    if (length(still_missing) > 0) {
      cat(sprintf("\nERROR: Failed to install: %s\n", paste(still_missing, collapse = ", ")))
      quit(status = 1)
    } else {
      cat("\nAll missing packages installed successfully.\n")
      quit(status = 0)
    }
  } else {
    cat("\nERROR: Missing dependencies detected. Run with --install to install them.\n")
    quit(status = 1)
  }
} else {
  cat("\nAll dependencies are installed and available.\n")
  quit(status = 0)
}
