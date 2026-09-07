# =============================================================================
# Unit Tests: Current Release Metadata Consistency
# =============================================================================

test_that("Current release metadata agrees with VERSION", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  version <- trimws(readLines(file.path(repo_root, "VERSION"), n = 1L, warn = FALSE))
  citation <- readLines(file.path(repo_root, "CITATION.cff"), warn = FALSE)
  readme <- readLines(file.path(repo_root, "README.md"), warn = FALSE)
  changelog <- readLines(file.path(repo_root, "CHANGELOG.md"), warn = FALSE)

  citation_version <- sub("^version:[[:space:]]*", "", grep("^version:", citation, value = TRUE)[1])
  readme_version <- sub("^Current pipeline version: \\*\\*([^*]+)\\*\\*[.]$", "\\1",
                        grep("^Current pipeline version:", readme, value = TRUE)[1])
  changelog_version <- sub("^## \\[([^]]+)\\].*$", "\\1",
                           grep("^## \\[", changelog, value = TRUE)[1])

  expect_identical(citation_version, version)
  expect_identical(readme_version, version)
  expect_identical(changelog_version, version)
})

test_that("renv lock records the complete declared R environment", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  lock_path <- file.path(repo_root, "renv.lock")
  expect_true(file.exists(lock_path))
  lock <- jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
  expect_identical(lock$R$Version, "4.5.3")
  expect_true(all(c("renv", "testthat", "microeco") %in% names(lock$Packages)))
  expect_identical(lock$Packages$microeco$Version, "2.3.0")
})
