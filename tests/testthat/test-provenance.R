# =============================================================================
# Source and lock provenance tests
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "manifest.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "dependencies.R"))
source(file.path("..", "..", "analysis", "utils", "preflight.R"))

test_that("source digest ignores generated Python bytecode", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  pycache <- file.path(repo_root, "analysis", "utils", "__pycache__")
  dir.create(pycache, recursive = TRUE, showWarnings = FALSE)
  generated <- file.path(pycache, "digest-stability-test.pyc")
  on.exit(unlink(generated, force = TRUE), add = TRUE)
  before <- source_provenance(repo_root)
  writeBin(as.raw(c(0xCA, 0xFE, 0xBA, 0xBE)), generated)
  after <- source_provenance(repo_root)
  expect_identical(after$source_digest_sha256, before$source_digest_sha256)
  expect_false(generated %in% as.character(before$source_files))
})

test_that("lock status reports the full lock closure and active library paths", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(requireNamespace("renv", quietly = TRUE))
  renv::activate(project = repo_root)
  project_library <- normalizePath(renv::paths$library(project = repo_root), winslash = "/")
  .libPaths(unique(c(project_library, .libPaths())))
  status <- read_lock_status(repo_root, RUNTIME_PACKAGES)
  expect_true(length(status$package_locations) >= length(RUNTIME_PACKAGES))
  expect_true(length(status$library_paths) >= 1L)
  expect_true(project_library %in% unlist(status$library_paths, use.names = FALSE))
})
