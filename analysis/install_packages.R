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

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath("analysis", winslash = "/", mustWork = TRUE)
}
source(file.path(script_dir, "utils", "dependencies.R"))
include_modules <- any(c("--all", "--modules", "--faprotax") %in% args)
REQUIRED_PACKAGES <- get_required_packages(include_tests = TRUE, include_modules = include_modules)

installed <- rownames(installed.packages())
missing_pkgs <- setdiff(REQUIRED_PACKAGES, installed)

cat("=== ONT wf-16s Pipeline Dependency Check ===\n")
cat(sprintf("Total required packages: %d\n", length(REQUIRED_PACKAGES)))
cat(sprintf("Installed: %d\n", length(REQUIRED_PACKAGES) - length(missing_pkgs)))
cat(sprintf("Missing:   %d\n\n", length(missing_pkgs)))

incompatible_pkgs <- character(0)

for (pkg in REQUIRED_PACKAGES) {
  is_inst <- pkg %in% installed
  ver <- if (is_inst) as.character(packageVersion(pkg)) else "NOT INSTALLED"
  status <- if (is_inst) "[OK]" else "[MISSING]"
  if (is_inst && pkg == "microeco") {
    if (utils::compareVersion(ver, "2.3.0") < 0) {
      status <- "[INCOMPATIBLE]"
      incompatible_pkgs <- c(incompatible_pkgs,
                             sprintf("microeco >= 2.3.0 required, found %s", ver))
    } else {
      env <- new.env(parent = emptyenv())
      utils::data("prok_func_FAPROTAX", package = "microeco", envir = env)
      db <- env[["prok_func_FAPROTAX"]]
      db_ver <- if (is.list(db) && !is.null(db[["ver"]])) as.character(db[["ver"]]) else "unknown"
      if (!identical(db_ver, "1.2.12")) {
        status <- "[INCOMPATIBLE]"
        incompatible_pkgs <- c(incompatible_pkgs,
                               sprintf("embedded FAPROTAX 1.2.12 required, found %s", db_ver))
      } else {
        ver <- sprintf("%s (FAPROTAX %s)", ver, db_ver)
      }
    }
  }
  cat(sprintf("  %-15s %-15s %s\n", pkg, status, ver))
}

has_missing <- length(missing_pkgs) > 0
has_incompatible <- length(incompatible_pkgs) > 0

if (has_missing || has_incompatible) {
  if (has_incompatible) {
    cat(sprintf("\nERROR: Incompatible dependencies detected:\n  %s\n",
                paste(incompatible_pkgs, collapse = "\n  ")))
  }
  if (has_missing) {
    if (do_install) {
      cat(sprintf("\nAttempting installation of %d missing packages...\n", length(missing_pkgs)))
      repos <- "https://cloud.r-project.org"
      install.packages(missing_pkgs, repos = repos)

      # Re-verify
      installed_now <- rownames(installed.packages())
      still_missing <- setdiff(REQUIRED_PACKAGES, installed_now)
      if (length(still_missing) > 0 || has_incompatible) {
        if (length(still_missing) > 0) {
          cat(sprintf("\nERROR: Failed to install: %s\n", paste(still_missing, collapse = ", ")))
        }
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
    quit(status = 1)
  }
} else {
  cat("\nAll dependencies are installed and available.\n")
  quit(status = 0)
}
