# =============================================================================
# Cross-process transaction boundaries
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "atomic_io.R"))
source(file.path("helper-fixtures.R"))

repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
runner <- file.path(repo_root, "analysis", "00_run_pipeline.R")
rscript <- Sys.which("Rscript")

run_transaction_process <- function(args, wd, env = character(0)) {
  child_env <- Sys.getenv()
  child_env[["R_LIBS_USER"]] <- .libPaths()[1]
  child_env[["RENV_CONFIG_AUTOLOADER_ENABLED"]] <- "FALSE"
  if (length(env)) child_env[names(env)] <- env
  processx::run(
    rscript,
    c(runner, "--allow-dirty", "--allow-unlocked", args),
    wd = wd,
    env = child_env,
    timeout = 120000,
    error_on_status = FALSE,
    echo = FALSE
  )
}

write_transaction_config <- function(root, abundance, cache, assignments = NULL, output) {
  cfg <- get_default_config()
  cfg$input$abundance_table <- normalizePath(abundance, winslash = "/")
  cfg$input$params_json <- normalizePath(create_temp_params(root), winslash = "/")
  cfg$input$assignments <- assignments
  cfg$taxonomy$cache <- normalizePath(cache, winslash = "/")
  cfg$output$base_dir <- output
  path <- file.path(root, "transaction_config.yml")
  yaml::write_yaml(cfg, path)
  normalizePath(path, winslash = "/")
}

test_that("a held output lock rejects a subprocess without mutating the output tree", {
  root <- tempfile("output_lock_process_")
  dir.create(root)
  abundance <- create_temp_abundance(root, n_species = 3L, sample_names = "S1")
  cache <- file.path(root, "taxonomy.json")
  writeLines("{}", cache)
  output <- file.path(root, "owned output")
  dir.create(output)
  sentinel <- file.path(output, "owned_sentinel.txt")
  sentinel_bytes <- charToRaw("owned output must remain unchanged")
  writeBin(sentinel_bytes, sentinel)
  config <- write_transaction_config(root, abundance, cache, output = output)

  held_lock <- acquire_output_lock(output, timeout_ms = 1000)
  on.exit(release_output_lock(held_lock), add = TRUE)
  result <- run_transaction_process(c("--config", config, "--output-dir", output,
                                      "--modules", "qc"), tempdir())

  expect_gt(result$status, 0L)
  expect_match(paste(result$stderr, result$stdout), "E_OUTPUT_BUSY")
  expect_identical(readBin(sentinel, "raw", n = file.info(sentinel)$size), sentinel_bytes)
  expect_identical(list.files(output, recursive = TRUE, all.files = TRUE, no.. = TRUE),
                   "owned_sentinel.txt")
})

test_that("a later failure rolls back a deferred local taxonomy refresh and cleans transaction state", {
  root <- tempfile("taxonomy_rollback_process_")
  dir.create(root)
  lineage <- paste(c("Bacteria", "Bacillati", "Bacillota", "Bacilli", "Bacillales",
                     "Bacillaceae", "Bacillus", "Bacillus cereus"), collapse = ";")
  abundance <- file.path(root, "abundance.tsv")
  writeLines(c(
    "tax\tS1\ttotal",
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t1\t1",
    paste(lineage, "2", "2", sep = "\t")
  ), abundance)
  cache <- file.path(root, "taxonomy.json")
  lineage_parts <- strsplit(lineage, ";", fixed = TRUE)[[1]]
  parents <- stats::setNames(as.list(seq_len(7L)),
                             vapply(seq_len(7L), function(i) paste(lineage_parts[seq_len(i)], collapse = ";"), character(1)))
  jsonlite::write_json(parents, cache, auto_unbox = TRUE, pretty = TRUE)
  cache_before <- readBin(cache, "raw", n = file.info(cache)$size)

  assignments <- file.path(root, "S1_assignments.tsv")
  writeLines(c(
    "C\tread_1\t1386\t0|1500\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus cereus",
    "C\tread_2\t1386\t0|1500\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus cereus",
    "U\tread_3\t0\t1500\tUnclassified"
  ), assignments)
  output <- file.path(root, "refresh output")
  config <- write_transaction_config(root, abundance, cache,
                                     assignments = list(S1 = normalizePath(assignments, winslash = "/")),
                                     output = output)

  result <- run_transaction_process(
    c("--config", config, "--output-dir", output, "--modules", "kreport,qc",
      "--refresh-taxonomy"),
    tempdir(),
    env = c(WF16S_TEST_MODE = "1", WF16S_INJECT_MODULE_FAILURE = "qc")
  )

  expect_gt(result$status, 0L)
  expect_match(paste(result$stderr, result$stdout), "Injected failure after module 'qc'")
  expect_identical(readBin(cache, "raw", n = file.info(cache)$size), cache_before)
  expect_false(file.exists(get_taxonomy_journal_path(cache)))
  expect_false(file.exists(get_output_journal_path(output)))
  expect_length(list.files(dirname(cache), pattern = "^\\.wf16s_tax_backup_", all.files = TRUE), 0L)
})
