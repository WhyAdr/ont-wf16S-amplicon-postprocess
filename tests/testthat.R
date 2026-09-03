# =============================================================================
# Testthat Suite Runner for ONT wf-16s Post-Processing Pipeline
# =============================================================================

suppressMessages({
  library(testthat)
})

# Locate test directory relative to this script
test_dir_path <- file.path("tests", "testthat")
if (!dir.exists(test_dir_path)) {
  # If called from tests/
  test_dir_path <- "testthat"
}

cat("Running testthat suite in:", normalizePath(test_dir_path, winslash = "/"), "\n")
res <- testthat::test_dir(test_dir_path, reporter = testthat::ProgressReporter$new())

# Check results and exit with appropriate code
df_res <- as.data.frame(res)
n_failed <- sum(df_res$failed)
n_errors <- sum(df_res$error)

if (n_failed > 0 || n_errors > 0) {
  cat(sprintf("\nTEST SUITE FAILED: %d failures, %d errors\n", n_failed, n_errors))
  quit(status = 1)
} else {
  cat(sprintf("\nTEST SUITE PASSED: %d tests completed successfully.\n", nrow(df_res)))
  quit(status = 0)
}
