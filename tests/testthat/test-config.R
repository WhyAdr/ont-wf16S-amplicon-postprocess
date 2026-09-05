# =============================================================================
# Unit Tests: Configuration & Path Resolution
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))

test_that("load_config loads default config.yml and resolves relative paths to config dir", {
  config_path <- file.path("..", "..", "config.yml")
  expect_true(file.exists(config_path))

  cfg <- load_config(config_path)

  expect_equal(cfg$schema_version, 1L)
  expect_equal(cfg$mode, "auto")
  expect_equal(cfg$seed, 42L)

  # Base output directory must be derived
  expect_true(is.character(cfg$output$base_dir))
  expect_true(nzchar(cfg$output$base_dir))

  # All 7 module directories must be present and derived from base_dir
  expect_true("qc" %in% names(cfg$output$dirs))
  expect_true("alpha" %in% names(cfg$output$dirs))
  expect_true("beta" %in% names(cfg$output$dirs))
  expect_true("composition" %in% names(cfg$output$dirs))
  expect_true("ordination" %in% names(cfg$output$dirs))
  expect_true("shared_taxa" %in% names(cfg$output$dirs))
  expect_true("kreport" %in% names(cfg$output$dirs))

  expect_equal(cfg$output$dirs$qc, file.path(cfg$output$base_dir, "01_QC"))
  expect_equal(cfg$output$dirs$alpha, file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
  expect_equal(cfg$output$dirs$beta, file.path(cfg$output$base_dir, "03_Beta_Diversity"))
  expect_equal(cfg$output$dirs$composition, file.path(cfg$output$base_dir, "04_Taxa_Composition"))
  expect_equal(cfg$output$dirs$ordination, file.path(cfg$output$base_dir, "05_Ordination"))
  expect_equal(cfg$output$dirs$shared_taxa, file.path(cfg$output$base_dir, "06_Shared_Taxa"))
  expect_equal(cfg$output$dirs$kreport, file.path(cfg$output$base_dir, "07_Kreport"))
})

test_that("CLI --output-dir overrides base_dir and all derived paths", {
  config_path <- file.path("..", "..", "config.yml")
  override_dir <- file.path(tempdir(), "test_override_output")

  cfg <- load_config(config_path, cli_opts = list(output_dir = override_dir))

  expect_equal(normalizePath(cfg$output$base_dir, winslash = "/", mustWork = FALSE),
               normalizePath(override_dir, winslash = "/", mustWork = FALSE))
  expect_equal(cfg$output$dirs$qc, file.path(cfg$output$base_dir, "01_QC"))
  expect_equal(cfg$output$dirs$alpha, file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
})

test_that("load_config errors on missing config file", {
  expect_error(load_config("non_existent_config.yml"), "Configuration file not found")
})

test_that("unknown config keys fail closed", {
  bad <- get_default_config()
  expect_error(merge_config(bad, list(alhpa = list())), "Unknown configuration key.*alhpa")
})

test_that("taxonomy refresh requires explicit CLI opt-in", {
  tmp <- tempfile("refresh_config_")
  dir.create(tmp)
  cfg <- get_default_config()
  cfg$taxonomy$network_mode <- "refresh"
  path <- file.path(tmp, "config.yml")
  yaml::write_yaml(cfg, path)

  expect_error(load_config(path), "YAML cannot enable taxonomy refresh")
  resolved <- load_config(path, cli_opts = list(refresh_taxonomy = TRUE))
  expect_equal(resolved$taxonomy$network_mode, "refresh")
  expect_true(resolved$cli$refresh_taxonomy)
})
