# =============================================================================
# Process-level release behavior
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))

repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
runner <- file.path(repo_root, "analysis", "00_run_pipeline.R")
rscript <- Sys.which("Rscript")

write_process_config <- function(root, sample_names = "S1", mode = "auto",
                                 metadata = NULL, assignments = NULL,
                                 unresolved_policy = "warn") {
  abundance <- create_temp_abundance(root, n_species = 3, sample_names = sample_names)
  cache <- file.path(root, "taxonomy cache.json")
  jsonlite::write_json(list(), cache, auto_unbox = TRUE)
  cfg <- get_default_config()
  cfg$mode <- mode
  cfg$input$abundance_table <- normalizePath(abundance, winslash = "/")
  cfg$input$metadata <- metadata
  cfg$input$params_json <- NULL
  cfg$input$assignments <- assignments
  cfg$taxonomy$cache <- normalizePath(cache, winslash = "/")
  cfg$taxonomy$unresolved_policy <- unresolved_policy
  cfg$output$base_dir <- file.path(root, "configured output")
  path <- file.path(root, "process config.yml")
  yaml::write_yaml(cfg, path)
  normalizePath(path, winslash = "/")
}

run_pipeline_process <- function(args, wd) {
  processx::run(
    rscript,
    c(runner, args),
    wd = wd,
    error_on_status = FALSE,
    echo = FALSE
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
  manifest <- jsonlite::fromJSON(
    file.path(output, "run_manifest.json"), simplifyVector = FALSE
  )
  expect_equal(manifest$run_status, "failed")
  expect_equal(manifest$modules$kreport$status, "failed")
  expect_equal(manifest$modules$composition$status, "completed")
})
