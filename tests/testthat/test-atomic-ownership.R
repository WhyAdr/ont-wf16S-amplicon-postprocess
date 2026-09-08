# =============================================================================
# Output ownership migration and staging tests
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "manifest.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "preflight.R"))
source(file.path("..", "..", "analysis", "utils", "atomic_io.R"))

write_legacy_output <- function(root, add_extra = FALSE) {
  dir.create(root, recursive = TRUE)
  qc <- file.path(root, "01_QC", "qc.tsv")
  alpha <- file.path(root, "02_Alpha_Diversity", "alpha.tsv")
  dir.create(dirname(qc), recursive = TRUE)
  dir.create(dirname(alpha), recursive = TRUE)
  writeLines("qc", qc)
  writeLines("alpha", alpha)
  writeLines("resolved", file.path(root, "resolved_config.yml"))
  writeLines("session", file.path(root, "session_info.txt"))
  if (add_extra) writeLines("unexpected", file.path(root, "stale.tsv"))
  manifest <- list(
    pipeline = "ont-wf16s-postprocess",
    output_root = normalizePath(root, winslash = "/"),
    run_status = "completed",
    modules = list(
      qc = list(outputs = qc),
      alpha = list(outputs = alpha)
    )
  )
  jsonlite::write_json(manifest, file.path(root, "run_manifest.json"), auto_unbox = TRUE)
}

test_that("legacy ownership is derived and stale owned products are not preserved", {
  root <- tempfile("legacy_ownership_")
  write_legacy_output(root)
  prior <- validate_prior_output(root, overwrite = TRUE)
  expect_true(isTRUE(prior$ownership_migrated))
  expect_true(all(c("01_QC/qc.tsv", "02_Alpha_Diversity/alpha.tsv") %in%
                    unlist(prior$owned_outputs, use.names = FALSE)))

  stage <- prepare_run_staging(root)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  preserve_unowned_outputs(root, stage, prior)
  expect_false(file.exists(file.path(stage, "02_Alpha_Diversity", "alpha.tsv")))
})

test_that("legacy migration rejects physical files without declared ownership", {
  root <- tempfile("legacy_ownership_reject_")
  write_legacy_output(root, add_extra = TRUE)
  expect_error(validate_prior_output(root, overwrite = TRUE), "E_OUTPUT_MIGRATION_REQUIRED")
})

test_that("output overlap checks are bidirectional and include discovered bamstats", {
  root <- tempfile("output_overlap_")
  cfg <- get_default_config()
  cfg$output$base_dir <- file.path(root, "results")
  cfg$input$wf16s_output_root <- file.path(root, "results", "upstream")
  expect_error(validate_output_root(cfg, normalizePath("..", winslash = "/")), "E_OUTPUT_UNSAFE")

  cfg$input$wf16s_output_root <- NULL
  discovered <- file.path(root, "results", "sample", "bamstats.readstats.tsv.gz")
  expect_error(validate_output_root(cfg, normalizePath("..", winslash = "/"), discovered),
               "E_OUTPUT_UNSAFE")
})

test_that("atomic_replace preserves original file if writer fails", {
  tf <- tempfile("atomic_test_")
  writeLines("initial_content", tf)
  expect_error(
    atomic_replace(tf, function(temp) stop("simulated writer failure")),
    "simulated writer failure"
  )
  expect_true(file.exists(tf))
  expect_equal(readLines(tf), "initial_content")
})

test_that("atomic_replace preserves original file on hash mismatch", {
  tf <- tempfile("atomic_hash_test_")
  writeLines("initial_content", tf)
  expect_error(
    atomic_replace(tf, function(temp) writeLines("new_content", temp),
                   expected_sha256 = paste(rep("0", 64), collapse = "")),
    "failed expected SHA-256 verification"
  )
  expect_true(file.exists(tf))
  expect_equal(readLines(tf), "initial_content")
})

test_that("acquire_output_lock enforces exclusivity with E_OUTPUT_BUSY", {
  root <- tempfile("lock_test_")
  l1 <- acquire_output_lock(root, timeout_ms = 1000)
  on.exit(release_output_lock(l1), add = TRUE)
  expect_error(acquire_output_lock(root, timeout_ms = 50), "E_OUTPUT_BUSY")
  release_output_lock(l1)
  l2 <- acquire_output_lock(root, timeout_ms = 1000)
  release_output_lock(l2)
})

test_that("recover_publication_journal restores prior completed run when final is missing", {
  root <- tempfile("journal_rec_")
  parent <- dirname(root)
  backup <- tempfile("journal_bak_", tmpdir = parent)
  write_legacy_output(backup)
  write_publication_journal(root, stage = NULL, backup = backup, phase = "prior_moved")

  expect_false(dir.exists(root))
  expect_true(dir.exists(backup))

  recover_publication_journal(root)

  expect_true(dir.exists(root))
  expect_false(dir.exists(backup))
  expect_false(file.exists(get_output_journal_path(root)))
})

test_that("recover_publication_journal cleans up backup when final is already valid", {
  root <- tempfile("journal_clean_")
  parent <- dirname(root)
  backup <- tempfile("journal_bak_", tmpdir = parent)
  write_legacy_output(root)
  write_legacy_output(backup)
  write_publication_journal(root, stage = NULL, backup = backup, phase = "stage_published")

  expect_true(dir.exists(root))
  expect_true(dir.exists(backup))

  recover_publication_journal(root)

  expect_true(dir.exists(root))
  expect_false(dir.exists(backup))
  expect_false(file.exists(get_output_journal_path(root)))
})

test_that("recover_publication_journal fails closed with E_OUTPUT_RECOVERY_REQUIRED", {
  root <- tempfile("journal_fail_")
  parent <- dirname(root)
  backup <- tempfile("journal_bak_", tmpdir = parent)
  dir.create(backup, recursive = TRUE)
  # Backup has no valid manifest
  write_publication_journal(root, stage = NULL, backup = backup, phase = "prior_moved")

  expect_error(recover_publication_journal(root), "E_OUTPUT_RECOVERY_REQUIRED")
  unlink(backup, recursive = TRUE, force = TRUE)
  remove_publication_journal(root)
})
