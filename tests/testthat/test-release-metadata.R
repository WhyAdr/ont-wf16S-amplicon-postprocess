# =============================================================================
# Unit Tests: Current Release Metadata Consistency
# =============================================================================

test_that("Current release metadata agrees with VERSION", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  version_path <- file.path(repo_root, "VERSION")
  version_lines <- readLines(version_path, warn = FALSE)
  version <- if (length(version_lines) == 1L) trimws(version_lines[[1]]) else ""
  version_bytes <- readBin(version_path, what = "raw", n = file.info(version_path)$size)
  source(file.path(repo_root, "analysis", "utils", "version.R"))
  expect_length(version_lines, 1L)
  expect_match(version, "^[0-9]+[.][0-9]+[.][0-9]+$")
  expect_true(length(version_bytes) > 0L && identical(tail(version_bytes, 1L), as.raw(0x0a)))
  expect_identical(read_pipeline_version(version_path), version)
  citation <- readLines(file.path(repo_root, "CITATION.cff"), warn = FALSE)
  readme <- readLines(file.path(repo_root, "README.md"), warn = FALSE)
  changelog <- readLines(file.path(repo_root, "CHANGELOG.md"), warn = FALSE)

  citation_version <- sub("^version:[[:space:]]*", "", grep("^version:", citation, value = TRUE)[1])
  current_line <- grep("^Current pipeline version:", readme, value = TRUE)
  readme_header_version <- sub("^Current pipeline version: \\*\\*([^*]+)\\*\\*[.]$", "\\1",
                               current_line)
  readme_mentions <- unlist(regmatches(
    readme, gregexpr("Version [0-9]+[.][0-9]+[.][0-9]+", readme, perl = TRUE)
  ))
  readme_versions <- sub("^Version ", "", readme_mentions)
  changelog_versions <- sub("^## \\[([^]]+)\\].*$", "\\1",
                            grep("^## \\[", changelog, value = TRUE))
  release_versions <- changelog_versions[changelog_versions != "Unreleased"]

  expect_identical(citation_version, version)
  expect_identical(readme_header_version, version)
  expect_length(readme_versions, 4L)
  expect_true(all(readme_versions == version))
  expect_true(changelog_versions[1] %in% c("Unreleased", version))
  expect_lte(sum(changelog_versions == "Unreleased"), 1L)
  expect_identical(release_versions[1], version)
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
