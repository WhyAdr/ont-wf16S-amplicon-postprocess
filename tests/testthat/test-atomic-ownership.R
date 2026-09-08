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

write_revision2_output <- function(root, declared_root = root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  cfg_file <- file.path(root, "resolved_config.yml")
  session_file <- file.path(root, "session_info.txt")
  qc_file <- file.path(root, "01_QC", "read_qc_summary.tsv")
  dir.create(dirname(qc_file), recursive = TRUE, showWarnings = FALSE)
  writeLines("resolved", cfg_file)
  writeLines("session", session_file)
  writeLines("read_qc", qc_file)

  hash_file <- function(path) digest::digest(file = path, algo = "sha256")
  registry <- c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport", "faprotax")
  records <- stats::setNames(lapply(registry, function(module_name) {
    if (identical(module_name, "qc")) {
      return(list(
        status = "completed", outputs = json_array(normalizePath(
          file.path(declared_root, "01_QC", "read_qc_summary.tsv"), winslash = "/",
          mustWork = FALSE)),
        warnings = json_array(character(0)), error = NULL, reason = NULL,
        start_time = "2026-09-07T00:00:00Z", end_time = "2026-09-07T00:00:01Z",
        duration_seconds = 1
      ))
    }
    list(
      status = "not_run", outputs = json_array(character(0)), warnings = json_array(character(0)),
      error = NULL, reason = sprintf("Module '%s' was not requested.", module_name),
      start_time = NULL, end_time = NULL, duration_seconds = NULL
    )
  }), registry)
  artifact <- function(relative_path, path, producer_module = NULL) {
    list(
      relative_path = relative_path, size_bytes = as.numeric(file.info(path)$size),
      sha256 = hash_file(path), producer_module = producer_module
    )
  }
  fingerprint <- function(path) list(
    path = path, size_bytes = 1L, mtime_utc = "2026-09-07T00:00:00Z",
    sha256 = paste(rep("c", 64L), collapse = "")
  )
  manifest <- list(
    pipeline = "ont-wf16s-postprocess", pipeline_version = "0.4.3",
    schema_version = 2L, schema_revision = 2L, config_schema_version = 1L,
    run_status = "completed", start_time = "2026-09-07T00:00:00Z",
    end_time = "2026-09-07T00:00:01Z", duration_seconds = 1,
    project_name = "fixture", mode = "single", seed = 42L, config_file = "fixture.yml",
    git_commit = paste(rep("a", 40L), collapse = ""), git_dirty = FALSE,
    source_digest_sha256 = paste(rep("b", 64L), collapse = ""),
    output_root = normalizePath(declared_root, winslash = "/", mustWork = FALSE),
    samples = json_array("S1"),
    command = json_array(c("Rscript", "analysis/00_run_pipeline.R")),
    cli = list(validate_only = FALSE, keep_going = FALSE, overwrite = FALSE,
      allow_unlocked = FALSE, allow_dirty = FALSE, online_preflight = FALSE,
      refresh_taxonomy = FALSE, krona = FALSE, modules = json_array("qc")),
    inputs = list(config_file = fingerprint("fixture.yml"), abundance_table = fingerprint("abundance.tsv"),
      params_json = fingerprint("params.json"), taxonomy_cache = fingerprint("taxonomy.json"),
      assignments = json_array(list()), bamstats = json_array(list())),
    modules = records,
    owned_outputs = json_array(c("01_QC/read_qc_summary.tsv", "resolved_config.yml",
      "run_manifest.json", "session_info.txt")),
    preserved_unowned_outputs = json_array(character(0)),
    artifacts = json_array(list(
      artifact("01_QC/read_qc_summary.tsv", qc_file, "qc"),
      artifact("resolved_config.yml", cfg_file), artifact("session_info.txt", session_file)
    )),
    warnings = json_array(character(0)),
    environment = list(locked = TRUE, lock_status = "synchronized", lockfile = "renv.lock",
      lockfile_sha256 = paste(rep("a", 64L), collapse = ""),
      r_discrepancies = json_array(character(0)), package_discrepancies = json_array(character(0)),
      library_discrepancies = json_array(character(0)), library_paths = json_array("C:/R/library"),
      project_library = "C:/R/library",
      package_locations = json_array(list(list(
        package = "yaml", version = "2.3.10", path = "C:/R/library",
        description_sha256 = NULL, lock_source = NULL, lock_repository = NULL,
        lock_remote_type = NULL, lock_remote_host = NULL, lock_remote_repo = NULL,
        lock_remote_ref = NULL, lock_remote_sha = NULL, lock_hash = NULL
      )))),
    package_versions = json_array(list(list(package = "yaml", version = "2.3.10")))
  )
  validate_manifest_v2(manifest, physical_root = root)
  jsonlite::write_json(manifest, file.path(root, "run_manifest.json"), auto_unbox = TRUE,
                       pretty = TRUE, null = "null")
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

test_that("atomic_replace restores original bytes after injected post-publish corruption", {
  tf <- tempfile("atomic_post_publish_")
  original <- as.raw(c(0x00, 0x01, 0xfe, 0xff, 0x41))
  writeBin(original, tf)

  withr::with_envvar(c(
    WF16S_TEST_MODE = "1",
    WF16S_INJECT_ATOMIC_POST_PUBLISH_CORRUPTION = basename(tf)
  ), {
    expect_error(
      atomic_replace(tf, function(temp) writeBin(as.raw(c(0x42, 0x43, 0x44)), temp)),
      "Published file .* failed (size|hash) verification"
    )
  })

  expect_identical(readBin(tf, "raw", n = file.info(tf)$size), original)
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

test_that("filesystem identity folds case only on Windows", {
  upper <- "C:/case-sensitive/Output"
  lower <- "C:/case-sensitive/output"
  expect_false(paths_are_same(upper, lower, os_type = "unix"))
  expect_true(paths_are_same(upper, lower, os_type = "windows"))
  expect_false(path_is_descendant("/tmp/Stage/file.tsv", "/tmp/stage", os_type = "unix"))
  expect_true(path_is_descendant("C:/tmp/Stage/file.tsv", "C:/tmp/stage", os_type = "windows"))
})

test_that("taxonomy journal recovery is content-addressed and fail-closed", {
  root <- tempfile("taxonomy_txn_")
  dir.create(root)
  cache <- file.path(root, "taxonomy.json")
  backup <- file.path(root, ".taxonomy.backup")
  output <- file.path(root, "output")
  original <- charToRaw("{\"original\":1}\n")
  candidate <- charToRaw("{\"candidate\":2}\n")
  writeBin(original, cache)
  writeBin(original, backup)
  original_hash <- compute_file_hash(cache)
  candidate_file <- file.path(root, "candidate.json")
  writeBin(candidate, candidate_file)
  candidate_hash <- compute_file_hash(candidate_file)

  write_taxonomy_journal(cache, backup, output, original_hash, candidate_hash, "prepared")
  expect_true(recover_taxonomy_journal(cache))
  expect_false(file.exists(get_taxonomy_journal_path(cache)))
  expect_false(file.exists(backup))
  expect_identical(readBin(cache, "raw", n = file.info(cache)$size), original)

  writeBin(original, backup)
  writeBin(candidate, cache)
  write_taxonomy_journal(cache, backup, output, original_hash, candidate_hash,
                         "candidate_committed")
  expect_true(recover_taxonomy_journal(cache))
  expect_identical(readBin(cache, "raw", n = file.info(cache)$size), original)
  expect_false(file.exists(get_taxonomy_journal_path(cache)))

  writeBin(original, backup)
  writeLines("external revision", cache)
  write_taxonomy_journal(cache, backup, output, original_hash, candidate_hash,
                         "candidate_committed")
  expect_error(recover_taxonomy_journal(cache), "E_TAXONOMY_CACHE_CHANGED")
  expect_true(file.exists(get_taxonomy_journal_path(cache)))
  expect_true(file.exists(backup))
  remove_file_checked(get_taxonomy_journal_path(cache), "test journal")
  remove_file_checked(backup, "test backup")
})

test_that("recover_publication_journal restores prior completed run when final is missing", {
  root <- tempfile("journal_rec_")
  parent <- dirname(root)
  backup <- tempfile("journal_bak_", tmpdir = parent)
  write_revision2_output(backup, declared_root = root)
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
  write_revision2_output(root)
  write_revision2_output(backup, declared_root = root)
  write_publication_journal(root, stage = NULL, backup = backup, phase = "stage_published")

  expect_true(dir.exists(root))
  expect_true(dir.exists(backup))

  recover_publication_journal(root)

  expect_true(dir.exists(root))
  expect_false(dir.exists(backup))
  expect_false(file.exists(get_output_journal_path(root)))
})

test_that("recovery rejects minimal or tampered completed manifests", {
  root <- tempfile("journal_strict_manifest_")
  parent <- dirname(root)
  backup <- tempfile("journal_bak_", tmpdir = parent)
  write_legacy_output(root)
  write_revision2_output(backup)
  write_publication_journal(root, stage = NULL, backup = backup, phase = "stage_published")
  expect_error(recover_publication_journal(root), "E_OUTPUT_RECOVERY_REQUIRED")
  expect_true(dir.exists(backup))
  remove_publication_journal(root)

  unlink(root, recursive = TRUE, force = TRUE)
  write_revision2_output(root)
  writeLines("tampered", file.path(root, "01_QC", "read_qc_summary.tsv"))
  write_publication_journal(root, stage = NULL, backup = backup, phase = "stage_published")
  expect_error(recover_publication_journal(root), "E_OUTPUT_RECOVERY_REQUIRED")
  expect_true(dir.exists(backup))
  remove_publication_journal(root)

  unlink(root, recursive = TRUE, force = TRUE)
  write_revision2_output(root)
  writeLines("undeclared", file.path(root, "undeclared.tsv"))
  write_publication_journal(root, stage = NULL, backup = backup, phase = "stage_published")
  expect_error(recover_publication_journal(root), "E_OUTPUT_RECOVERY_REQUIRED")
  expect_true(dir.exists(backup))
  remove_publication_journal(root)
  unlink(backup, recursive = TRUE, force = TRUE)
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

test_that("prepare, publish, and cleanup module staging work correctly", {
  root <- tempfile("mod_stage_root_")
  stage <- file.path(root, "run_stage")
  dir.create(stage, recursive = TRUE)

  mod_stage <- prepare_module_staging(stage, "qc")
  expect_true(dir.exists(mod_stage))

  out_dir <- file.path(mod_stage, "01_QC")
  dir.create(out_dir)
  test_file <- file.path(out_dir, "qc.tsv")
  writeLines("qc content", test_file)

  publish_module_staging(mod_stage, stage, test_file)
  expect_true(file.exists(file.path(stage, "01_QC", "qc.tsv")))
  expect_false(dir.exists(mod_stage))
})

test_that("publish_module_staging rejects collisions without changing either file", {
  root <- tempfile("module_stage_collision_")
  stage <- file.path(root, "run_stage")
  dir.create(file.path(stage, "01_QC"), recursive = TRUE)
  staged_file <- file.path(stage, "01_QC", "qc.tsv")
  writeBin(charToRaw("prior-stage-bytes"), staged_file)

  module_stage <- prepare_module_staging(stage, "qc")
  module_file <- file.path(module_stage, "01_QC", "qc.tsv")
  dir.create(dirname(module_file), recursive = TRUE)
  writeBin(charToRaw("module-stage-bytes"), module_file)

  expect_error(publish_module_staging(module_stage, stage, module_file), "E_OUTPUT_COLLISION")
  expect_identical(readBin(staged_file, "raw", n = file.info(staged_file)$size), charToRaw("prior-stage-bytes"))
  expect_identical(readBin(module_file, "raw", n = file.info(module_file)$size), charToRaw("module-stage-bytes"))
  cleanup_module_staging(module_stage)
})

test_that("verify_physical_file_census enforces exact file set and detects mismatches", {
  stage <- tempfile("census_stage_")
  dir.create(file.path(stage, "01_QC"), recursive = TRUE)
  writeLines("qc", file.path(stage, "01_QC", "qc.tsv"))
  writeLines("session", file.path(stage, "session_info.txt"))
  writeLines("unowned", file.path(stage, "user_note.txt"))

  owned <- c("01_QC/qc.tsv", "session_info.txt", "run_manifest.json")
  preserved <- c("user_note.txt")

  expect_no_error(verify_physical_file_census(stage, owned, preserved))

  writeLines("extra", file.path(stage, "extra.txt"))
  expect_error(verify_physical_file_census(stage, owned, preserved), "E_CENSUS_MISMATCH")
  unlink(file.path(stage, "extra.txt"))

  expect_error(verify_physical_file_census(stage, c(owned, "missing.tsv"), preserved), "E_CENSUS_MISMATCH")
})
