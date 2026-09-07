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

  # All module directories must be present and derived from base_dir
  expect_true("qc" %in% names(cfg$output$dirs))
  expect_true("alpha" %in% names(cfg$output$dirs))
  expect_true("beta" %in% names(cfg$output$dirs))
  expect_true("composition" %in% names(cfg$output$dirs))
  expect_true("ordination" %in% names(cfg$output$dirs))
  expect_true("shared_taxa" %in% names(cfg$output$dirs))
  expect_true("kreport" %in% names(cfg$output$dirs))
  expect_true("faprotax" %in% names(cfg$output$dirs))

  expect_equal(cfg$output$dirs$qc, file.path(cfg$output$base_dir, "01_QC"))
  expect_equal(cfg$output$dirs$alpha, file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
  expect_equal(cfg$output$dirs$beta, file.path(cfg$output$base_dir, "03_Beta_Diversity"))
  expect_equal(cfg$output$dirs$composition, file.path(cfg$output$base_dir, "04_Taxa_Composition"))
  expect_equal(cfg$output$dirs$ordination, file.path(cfg$output$base_dir, "05_Ordination"))
  expect_equal(cfg$output$dirs$shared_taxa, file.path(cfg$output$base_dir, "06_Shared_Taxa"))
  expect_equal(cfg$output$dirs$kreport, file.path(cfg$output$base_dir, "07_Kreport"))
  expect_equal(cfg$output$dirs$faprotax, file.path(cfg$output$base_dir, "08_FAPROTAX"))
  expect_equal(cfg$faprotax$top_n_functions, 20L)
  expect_false(cfg$krona$enabled)
  expect_true(is.logical(cfg$krona$render_html))
  expect_true(cfg$krona$render_html)
  expect_equal(cfg$krona$executable, "ktImportText")
  expect_false(cfg$cli$krona)
  expect_false("faprotax" %in% cfg$cli$modules)
})

test_that("Krona CLI opt-in is recorded in config and manifest settings", {
  config_path <- file.path("..", "..", "config.yml")
  cfg <- load_config(config_path, cli_opts = list(krona = TRUE))

  expect_true(cfg$krona$enabled)
  expect_true(cfg$cli$krona)
})

test_that("Requested module parsing preserves order and rejects invalid requests", {
  expect_equal(
    parse_requested_modules(" qc,alpha shared\tkreport "),
    c("qc", "alpha", "shared", "kreport")
  )
  expect_error(parse_requested_modules("qc,qc"), "Duplicate module name")
  expect_error(parse_requested_modules("   "), "must contain at least one")
})

test_that("Configuration YAML root must be a named mapping", {
  root <- tempfile("config_root_")
  dir.create(root)
  scalar <- file.path(root, "scalar.yml")
  empty <- file.path(root, "empty.yml")
  writeLines("not-a-mapping", scalar)
  writeLines(character(0), empty)
  expect_error(load_config(scalar), "named mapping")
  expect_error(load_config(empty), "named mapping")
})

test_that("CLI --output-dir overrides base_dir and all derived paths", {
  config_path <- file.path("..", "..", "config.yml")
  override_dir <- file.path(tempdir(), "test_override_output")

  cfg <- load_config(config_path, cli_opts = list(output_dir = override_dir))

  expect_equal(normalizePath(cfg$output$base_dir, winslash = "/", mustWork = FALSE),
               normalizePath(override_dir, winslash = "/", mustWork = FALSE))
  expect_equal(cfg$output$dirs$qc, file.path(cfg$output$base_dir, "01_QC"))
  expect_equal(cfg$output$dirs$alpha, file.path(cfg$output$base_dir, "02_Alpha_Diversity"))

  # Relative CLI output_dir must resolve against getwd(), not config directory
  cfg_rel <- load_config(config_path, cli_opts = list(output_dir = "relative_test_out"))
  expect_equal(normalizePath(cfg_rel$output$base_dir, winslash = "/", mustWork = FALSE),
               normalizePath(file.path(getwd(), "relative_test_out"), winslash = "/", mustWork = FALSE))
})

test_that("load_config errors on missing config file", {
  expect_error(load_config("non_existent_config.yml"), "Configuration file not found")
})

test_that("unknown config keys fail closed", {
  bad <- get_default_config()
  expect_error(merge_config(bad, list(alhpa = list())), "Unknown configuration key.*alhpa")
  expect_error(
    merge_config(bad, list(krona = list(renderer = "unexpected"))),
    "Unknown configuration key.*krona.renderer"
  )
})

test_that("Krona configuration values fail closed", {
  bad_enabled <- get_default_config()
  bad_enabled$krona$enabled <- "yes"
  expect_error(validate_config(bad_enabled), "krona.enabled")

  bad_render <- get_default_config()
  bad_render$krona$render_html <- NA
  expect_error(validate_config(bad_render), "krona.render_html")

  bad_executable <- get_default_config()
  bad_executable$krona$executable <- "  "
  expect_error(validate_config(bad_executable), "krona.executable")
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
