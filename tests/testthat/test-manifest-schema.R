# =============================================================================
# Unit Tests: Manifest Schema v2
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "manifest.R"))

make_manifest_fixture <- function(samples, modules = c("qc"), warnings = character(0)) {
  output_root <- tempfile("manifest_fixture_root_")
  dir.create(output_root)
  writeLines("resolved", file.path(output_root, "resolved_config.yml"))
  writeLines("session", file.path(output_root, "session_info.txt"))
  registry <- c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport", "faprotax")
  records <- stats::setNames(lapply(registry, function(module_name) {
    if (!(module_name %in% modules)) {
      return(manifest_module_record(new_not_run_module_record(sprintf(
        "Module '%s' was not requested.", module_name))))
    }
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
  }), registry)
  list(
    pipeline = "ont-wf16s-postprocess",
    pipeline_version = "0.4.2",
    schema_version = 2L,
    schema_revision = 1L,
    config_schema_version = 1L,
    run_status = "completed",
    start_time = "2026-09-07T00:00:00Z",
    end_time = "2026-09-07T00:00:01Z",
    duration_seconds = 1,
    project_name = "fixture",
    mode = "single",
    seed = 42L,
    config_file = "fixture.yml",
    git_commit = paste(rep("a", 40L), collapse = ""),
    git_dirty = FALSE,
    source_digest_sha256 = paste(rep("b", 64L), collapse = ""),
    output_root = normalizePath(output_root, winslash = "/"),
    samples = json_array(samples),
    command = json_array(c("Rscript", "analysis/00_run_pipeline.R")),
    cli = list(validate_only = FALSE, keep_going = FALSE, overwrite = FALSE,
               allow_unlocked = FALSE, allow_dirty = FALSE, online_preflight = FALSE,
               refresh_taxonomy = FALSE, krona = FALSE, modules = json_array(modules)),
    inputs = list(assignments = json_array(list()), bamstats = json_array(list())),
    owned_outputs = json_array(c("resolved_config.yml", "session_info.txt")),
    modules = records,
    warnings = json_array(warnings),
    environment = list(
      locked = TRUE,
      lock_status = "synchronized",
      lockfile = "renv.lock",
      lockfile_sha256 = paste(rep("a", 64L), collapse = ""),
      r_discrepancies = json_array(character(0)),
      package_discrepancies = json_array(character(0)),
      library_discrepancies = json_array(character(0)),
      library_paths = json_array(c("C:/R/library")),
      project_library = "C:/R/library",
      package_locations = json_array(list(list(package = "yaml", version = "2.3.10", path = "C:/R/library/yaml")))
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

test_that("manifest v2 enforces the complete module registry and not_run state", {
  manifest <- make_manifest_fixture("S1")
  manifest$modules$faprotax <- NULL
  expect_error(validate_manifest_v2(manifest), "maintained registry")

  manifest <- make_manifest_fixture("S1")
  manifest$modules$alpha$outputs <- json_array("stale.tsv")
  expect_error(validate_manifest_v2(manifest), "not_run")

  manifest <- make_manifest_fixture("S1")
  manifest$modules$alpha$reason <- NULL
  expect_error(validate_manifest_v2(manifest), "invalid result object")
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
