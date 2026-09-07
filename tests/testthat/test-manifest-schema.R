# =============================================================================
# Unit Tests: Manifest Schema v2
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "manifest.R"))

make_manifest_fixture <- function(samples, modules = c("qc"), warnings = character(0)) {
  records <- stats::setNames(lapply(modules, function(module_name) {
    list(
      status = "completed",
      outputs = json_array(character(0)),
      warnings = json_array(character(0)),
      error = NULL,
      reason = NULL,
      start_time = "2026-09-07T00:00:00Z",
      end_time = "2026-09-07T00:00:01Z",
      duration_seconds = 1
    )
  }), modules)
  list(
    schema_version = 2L,
    config_schema_version = 1L,
    samples = json_array(samples),
    command = json_array(c("Rscript", "analysis/00_run_pipeline.R")),
    cli = list(modules = json_array(modules)),
    inputs = list(assignments = json_array(list()), bamstats = NULL),
    modules = records,
    warnings = json_array(warnings),
    environment = list(
      locked = TRUE,
      lockfile = "renv.lock",
      lockfile_sha256 = paste(rep("a", 64L), collapse = "")
    ),
    package_versions = json_array(list(list(package = "yaml", version = "2.3.10")))
  )
}

test_that("manifest v2 retains array shape across cardinalities", {
  root <- tempfile("manifest_schema_")
  dir.create(root)
  for (samples in list("S1", c("S1", "S2"))) {
    path <- file.path(root, paste0("manifest_", length(samples), ".json"))
    write_manifest_v2(make_manifest_fixture(samples), path)
    parsed <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    expect_true(is_json_array(parsed$samples))
    expect_true(is_json_array(parsed$command))
    expect_true(is_json_array(parsed$cli$modules))
    expect_true(is_json_array(parsed$warnings))
    expect_true(is_json_array(parsed$package_versions))
    expect_true(is_json_array(parsed$inputs$assignments))
    expect_true(is_json_array(parsed$modules$qc$outputs))
    expect_true(is_json_array(parsed$modules$qc$warnings))
  }
})

test_that("manifest v2 permits explicit not_run records and rejects scalar arrays", {
  manifest <- make_manifest_fixture("S1")
  manifest$modules$alpha <- manifest_module_record(new_not_run_module_record("Pipeline stopped."))
  expect_no_error(validate_manifest_v2(manifest))

  manifest$samples <- jsonlite::unbox("S1")
  expect_error(validate_manifest_v2(manifest), "samples.*JSON array")
})

test_that("Krona provenance sample records retain array shape", {
  root <- tempfile("krona_schema_")
  dir.create(root)
  for (records in list(
    list(list(sample_id = "S1", tsv_path = "S1.krona.tsv")),
    list(list(sample_id = "S1", tsv_path = "S1.krona.tsv"),
         list(sample_id = "S2", tsv_path = "S2.krona.tsv"))
  )) {
    path <- file.path(root, paste0("krona_", length(records), ".json"))
    jsonlite::write_json(list(samples = json_array(records)), path,
                         auto_unbox = TRUE, null = "null")
    parsed <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    expect_true(is_json_array(parsed$samples))
  }
})
