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
