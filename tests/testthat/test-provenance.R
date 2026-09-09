# =============================================================================
# Source and lock provenance tests
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "manifest.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "dependencies.R"))
source(file.path("..", "..", "analysis", "utils", "preflight.R"))

test_that("provenance path identity preserves Unix case distinctions", {
  expect_false(preflight_paths_are_same("/cache/pkg", "/cache/Pkg", os_type = "unix"))
  expect_true(preflight_paths_are_same("C:/cache/pkg", "C:/cache/Pkg", os_type = "windows"))
  expect_false(preflight_path_is_same_or_descendant("/cache/Pkg", "/cache/pkg", os_type = "unix"))
  expect_true(preflight_path_is_same_or_descendant("C:/cache/Pkg", "C:/cache/pkg", os_type = "windows"))
})

test_that("renv-cache locations must resolve from the matching project entry", {
  cache <- "/renv/cache"
  project <- "/project/renv/library"
  loaded <- "/renv/cache/yaml/2.3.12/hash/yaml"
  expect_true(package_location_is_project_bound(
    loaded, project, cache, resolved_project_entry = loaded, os_type = "unix"
  ))
  expect_false(package_location_is_project_bound(
    loaded, project, cache,
    resolved_project_entry = "/renv/cache/yaml/2.3.11/other/yaml",
    os_type = "unix"
  ))
  expect_false(package_location_is_project_bound(
    "/external/yaml", project, cache,
    resolved_project_entry = "/external/yaml", os_type = "unix"
  ))
})

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

test_that("pinned Krona renderer assets are included in source provenance", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  relative <- as.character(source_provenance(repo_root)$source_files)
  expected <- file.path(
    "analysis", "vendor", "krona-2.8.1",
    c("SOURCE.json", "LICENSE.txt", "src/krona-2.0.js", "img/favicon.ico",
      "img/hidden.png", "img/loading.gif", "img/logo-med.png")
  )
  expect_true(all(gsub("\\\\", "/", expected) %in% relative))
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
  location <- status$package_locations[[1]]
  expect_true(all(c("description_sha256", "lock_source", "lock_repository", "lock_remote_type",
                    "lock_remote_host", "lock_remote_repo", "lock_remote_ref", "lock_remote_sha",
                    "lock_hash") %in% names(location)))
  if (!is.null(location$description_sha256)) expect_match(location$description_sha256, "^[0-9a-f]{64}$")
})

test_that("input inventory fingerprints the supplied configuration file", {
  config_file <- tempfile(fileext = ".yml")
  writeLines("project_name: provenance_fixture", config_file)
  cfg <- get_default_config()
  cfg$config_file <- normalizePath(config_file, winslash = "/", mustWork = TRUE)
  cfg$input$abundance_table <- config_file
  cfg$input$params_json <- config_file
  cfg$taxonomy$cache <- config_file
  cfg$input$assignments <- list()

  inventory <- inventory_inputs(cfg)
  expect_identical(inventory$config_file$path, cfg$config_file)
  expect_match(inventory$config_file$sha256, "^[0-9a-f]{64}$")
  writeLines("project_name: changed_fixture", config_file)
  expect_error(assert_inputs_unchanged(inventory), "E_INPUT_CHANGED")
})

test_that("startup files are included in maintained source files", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  files <- maintained_source_files(repo_root)
  norm_files <- tolower(normalizePath(files, winslash = "/", mustWork = FALSE))
  rprofile <- tolower(normalizePath(file.path(repo_root, ".Rprofile"), winslash = "/", mustWork = FALSE))
  activate <- tolower(normalizePath(file.path(repo_root, "renv", "activate.R"), winslash = "/", mustWork = FALSE))
  settings <- tolower(normalizePath(file.path(repo_root, "renv", "settings.json"), winslash = "/", mustWork = FALSE))
  expect_true(rprofile %in% norm_files)
  expect_true(activate %in% norm_files)
  expect_true(settings %in% norm_files)
})

test_that("package resolved outside project library is recorded as library discrepancy", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  temp_lib <- tempfile("mock_lib_")
  dir.create(temp_lib, recursive = TRUE)
  fake_repo <- tempfile("mock_repo_")
  dir.create(fake_repo, recursive = TRUE)
  lock_file <- file.path(fake_repo, "renv.lock")
  writeLines(jsonlite::toJSON(list(
    R = list(Version = as.character(getRversion())),
    Packages = list(
      yaml = list(Package = "yaml", Version = as.character(utils::packageVersion("yaml")))
    )
  ), auto_unbox = TRUE), lock_file)

  status <- read_lock_status(fake_repo, c("yaml"))
  expect_identical(status$lock_status, "mismatch")
  expect_true(length(status$library_discrepancies) >= 1L)
  expect_true(any(grepl("package loaded from outside project library", unlist(status$library_discrepancies))))
})
