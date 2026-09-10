# =============================================================================
# Process-level release behavior
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))

repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
runner <- file.path(repo_root, "analysis", "00_run_pipeline.R")
rscript <- Sys.which("Rscript")

write_process_config <- function(root, sample_names = "S1", mode = "auto",
                                 metadata = NULL, assignments = NULL,
                                 unresolved_policy = "warn", classifier = "minimap2") {
  abundance <- create_temp_abundance(root, n_species = 3, sample_names = sample_names)
  cache <- file.path(root, "taxonomy cache.json")
  writeLines("{}", cache)
  cfg <- get_default_config()
  cfg$mode <- mode
  cfg$input$abundance_table <- normalizePath(abundance, winslash = "/")
  cfg$input$metadata <- metadata
  cfg$input$params_json <- normalizePath(
    create_temp_params(root, classifier = classifier), winslash = "/"
  )
  cfg$input$assignments <- assignments
  cfg$taxonomy$cache <- normalizePath(cache, winslash = "/")
  cfg$taxonomy$unresolved_policy <- unresolved_policy
  cfg$output$base_dir <- file.path(root, "configured output")
  path <- file.path(root, "process config.yml")
  yaml::write_yaml(cfg, path)
  normalizePath(path, winslash = "/")
}

run_pipeline_process <- function(args, wd) {
  child_env <- Sys.getenv()
  child_env[["R_LIBS_USER"]] <- .libPaths()[1]
  child_env[["RENV_CONFIG_AUTOLOADER_ENABLED"]] <- "FALSE"
  processx::run(
    rscript,
    c(runner, "--allow-dirty", "--allow-unlocked", args),
    wd = wd,
    env = child_env,
    error_on_status = FALSE,
    echo = FALSE,
    timeout = 120000
  )
}

test_that("runner works outside the repository and validate-only is mutation-free", {
  root <- tempfile("process path with spaces ")
  dir.create(root)
  config <- write_process_config(root)
  output <- file.path(root, "validation output with spaces")
  expect_false(file.exists(output))
  result <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--validate-only"),
    wd = tempdir()
  )
  expect_equal(result$status, 0L, info = paste(result$stderr, result$stdout))
  expect_match(result$stdout, "Zero filesystem mutations performed")
  expect_false(file.exists(output))
})

test_that("invalid modules, assignments, and metadata return non-zero", {
  root <- tempfile("invalid_process_")
  dir.create(root)

  config <- write_process_config(root)
  invalid_module <- run_pipeline_process(
    c("--config", config, "--modules", "not_a_module"), tempdir()
  )
  expect_gt(invalid_module$status, 0L)
  expect_match(invalid_module$stderr, "Unknown module")

  bad_assignment <- file.path(root, "bad assignment.tsv")
  writeLines("C\tread1\t123\t1500", bad_assignment)
  config <- write_process_config(
    root,
    assignments = list(S1 = normalizePath(bad_assignment, winslash = "/"))
  )
  invalid_assignment <- run_pipeline_process(
    c("--config", config, "--validate-only"), tempdir()
  )
  expect_gt(invalid_assignment$status, 0L)
  expect_match(invalid_assignment$stderr, "expected exactly 5")

  bad_metadata <- file.path(root, "bad metadata.tsv")
  writeLines(c("SampleID\tGroup", "Wrong1\tA", "Wrong2\tB"), bad_metadata)
  config <- write_process_config(
    root,
    sample_names = c("S1", "S2"),
    mode = "cohort",
    metadata = normalizePath(bad_metadata, winslash = "/")
  )
  invalid_metadata <- run_pipeline_process(
    c("--config", config, "--validate-only"), tempdir()
  )
  expect_gt(invalid_metadata$status, 0L)
  expect_match(invalid_metadata$stderr, "Metadata SampleID mismatch")
})

test_that("invalid module requests fail before malformed inputs or output creation", {
  root <- tempfile("invalid_module_fast_fail_")
  dir.create(root)
  malformed <- file.path(root, "assignments.tsv")
  writeLines("malformed", malformed)
  config <- write_process_config(
    root, assignments = list(S1 = normalizePath(malformed, winslash = "/"))
  )
  output <- file.path(root, "should_not_exist")
  result <- run_pipeline_process(
    c("--config", config, "--modules", "not_a_module", "--output-dir", output), tempdir()
  )
  expect_gt(result$status, 0L)
  expect_match(result$stderr, "Unknown module")
  expect_false(grepl("expected exactly 5", result$stderr))
  expect_false(file.exists(output))
})

test_that("duplicate module requests fail during configuration", {
  root <- tempfile("duplicate_module_process_")
  dir.create(root)
  config <- write_process_config(root)
  result <- run_pipeline_process(c("--config", config, "--modules", "qc,qc"), tempdir())
  expect_gt(result$status, 0L)
  expect_match(result$stderr, "Duplicate module name")
})

test_that("Krona opt-in requires the kreport module before output mutation", {
  root <- tempfile("krona_dependency_")
  dir.create(root)
  config <- write_process_config(root)
  output <- file.path(root, "krona output")

  result <- run_pipeline_process(
    c("--config", config, "--output-dir", output,
      "--modules", "composition", "--krona", "--validate-only"),
    tempdir()
  )

  expect_gt(result$status, 0L)
  expect_match(result$stderr, "Krona export requires the 'kreport' module")
  expect_false(file.exists(output))
})

test_that("Kraken2 is rejected from params before assignment schema parsing", {
  root <- tempfile("kraken_contract_")
  dir.create(root)
  bad_assignment <- file.path(root, "six_field.tsv")
  writeLines("C\tread1\t123\t1500\tA:1\tBacteria", bad_assignment)
  config <- write_process_config(
    root, assignments = list(S1 = normalizePath(bad_assignment, winslash = "/")),
    classifier = "kraken2"
  )
  result <- run_pipeline_process(c("--config", config, "--validate-only"), tempdir())
  expect_gt(result$status, 0L)
  expect_match(result$stderr, "supports minimap2 only")
  expect_false(grepl("expected exactly 5", result$stderr))
})


test_that("keep-going records failure, executes later modules, and exits non-zero", {
  root <- tempfile("keep_going_")
  dir.create(root)
  config <- write_process_config(root, unresolved_policy = "error")
  output <- file.path(root, "keep going output")
  result <- run_pipeline_process(
    c("--config", config, "--output-dir", output,
      "--modules", "kreport,composition", "--keep-going"),
    tempdir()
  )
  expect_gt(result$status, 0L)
  expect_match(result$stderr, "E_KREPORT_PREFLIGHT")
  expect_false(file.exists(file.path(output, "run_manifest.json")))
})

test_that("fail-fast runs record later requested modules as not_run", {
  root <- tempfile("not_run_process_")
  dir.create(root)
  config <- write_process_config(root, unresolved_policy = "error")
  output <- file.path(root, "fail fast output")
  result <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "kreport,composition,alpha"),
    tempdir()
  )
  expect_gt(result$status, 0L)
  expect_match(result$stderr, "E_KREPORT_PREFLIGHT")
  expect_false(file.exists(file.path(output, "run_manifest.json")))
})

test_that("module failure rollback removes partial writes before publishing a failed manifest", {
  root <- tempfile("module_failure_rollback_")
  dir.create(root)
  config <- write_process_config(root)
  output <- file.path(root, "injected failure output")
  withr::local_envvar(WF16S_INJECT_MODULE_FAILURE = "qc", WF16S_TEST_MODE = "1")
  result <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "qc,alpha"), tempdir()
  )
  expect_gt(result$status, 0L)
  manifest_path <- file.path(output, "run_manifest.json")
  expect_true(file.exists(manifest_path), info = paste(result$stderr, result$stdout))
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  expect_identical(manifest$run_status, "failed")
  expect_identical(manifest$modules$qc$status, "failed")
  expect_length(manifest$modules$qc$outputs, 0L)
  expect_identical(manifest$modules$alpha$status, "not_run")
  physical <- list.files(output, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  physical <- physical[!dir.exists(file.path(output, physical))]
  expect_setequal(physical, c("resolved_config.yml", "session_info.txt", "run_manifest.json"))
})

test_that("failed first-run output retries only with --overwrite", {
  root <- tempfile("failed_first_run_retry_")
  dir.create(root)
  config <- write_process_config(root)
  output <- file.path(root, "retry output")
  failed <- withr::with_envvar(c(
    WF16S_INJECT_MODULE_FAILURE = "qc", WF16S_TEST_MODE = "1"
  ), run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "qc"), tempdir()
  ))
  expect_gt(failed$status, 0L)
  expect_identical(jsonlite::fromJSON(file.path(output, "run_manifest.json"), simplifyVector = FALSE)$run_status,
                   "failed")

  without_overwrite <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "qc"), tempdir()
  )
  expect_gt(without_overwrite$status, 0L)
  expect_match(without_overwrite$stderr, "non-empty; pass --overwrite")

  retried <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "qc", "--overwrite"), tempdir()
  )
  expect_equal(retried$status, 0L, info = paste(retried$stderr, retried$stdout))
  expect_identical(jsonlite::fromJSON(file.path(output, "run_manifest.json"), simplifyVector = FALSE)$run_status,
                   "completed")
})

test_that("keep-going continues after an injected module write and preserves only declared outputs", {
  root <- tempfile("module_failure_keep_going_")
  dir.create(root)
  config <- write_process_config(root)
  output <- file.path(root, "injected keep-going output")
  withr::local_envvar(WF16S_INJECT_MODULE_FAILURE = "qc", WF16S_TEST_MODE = "1")
  result <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "qc,alpha", "--keep-going"),
    tempdir()
  )
  expect_gt(result$status, 0L)
  manifest <- jsonlite::fromJSON(file.path(output, "run_manifest.json"), simplifyVector = FALSE)
  expect_identical(manifest$run_status, "failed")
  expect_identical(manifest$modules$qc$status, "failed")
  expect_identical(manifest$modules$alpha$status, "completed")
  owned <- unlist(manifest$owned_outputs, use.names = FALSE)
  physical <- list.files(output, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  physical <- physical[!dir.exists(file.path(output, physical))]
  expect_setequal(physical, owned)
})

test_that("failed overwrite preserves the previously completed output tree", {
  root <- tempfile("overwrite_rollback_")
  dir.create(root)
  config <- write_process_config(root)
  output <- file.path(root, "overwrite output")
  first <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "kreport"), tempdir()
  )
  expect_equal(first$status, 0L, info = paste(first$stderr, first$stdout))
  before <- readBin(file.path(output, "run_manifest.json"), "raw",
                    n = file.info(file.path(output, "run_manifest.json"))$size)
  withr::local_envvar(WF16S_INJECT_MODULE_FAILURE = "kreport", WF16S_TEST_MODE = "1")
  second <- run_pipeline_process(
    c("--config", config, "--output-dir", output, "--modules", "kreport", "--overwrite"),
    tempdir()
  )
  expect_gt(second$status, 0L)
  after <- readBin(file.path(output, "run_manifest.json"), "raw",
                   n = file.info(file.path(output, "run_manifest.json"))$size)
  expect_identical(after, before)
})
