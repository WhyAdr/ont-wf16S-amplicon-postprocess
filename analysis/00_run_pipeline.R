#!/usr/bin/env Rscript
# =============================================================================
# ONT wf-16s Amplicon Post-Processing Pipeline Runner
# =============================================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/") else
    normalizePath(file.path(getwd(), "analysis"), winslash = "/")
}

read_pipeline_version <- function(path) {
  if (!file.exists(path)) stop("VERSION is missing.", call. = FALSE)
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  text <- rawToChar(bytes)
  if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\n$", text)) {
    stop("VERSION must contain exactly one newline-terminated SemVer value.", call. = FALSE)
  }
  sub("\\n$", "", text)
}

main <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(dirname(script_dir), winslash = "/")
  options(wf16s.pipeline_root = repo_root)
  pipeline_version <- read_pipeline_version(file.path(repo_root, "VERSION"))

  fatal <- function(label, error) {
    stop(sprintf("[FATAL] %s: %s", label, conditionMessage(error)), call. = FALSE)
  }

  source(file.path(script_dir, "utils", "env_loader.R"))
  load_pipeline_environment(repo_root = repo_root, fatal_fn = fatal)

  source(file.path(script_dir, "utils", "dependencies.R"))
  check_dependencies()
  for (file in c("cli.R", "config.R", "io.R", "metrics.R", "plotting.R", "manifest.R",
                 "kreport.R", "atomic_io.R", "module_result.R", "preflight.R")) {
    source(file.path(script_dir, "utils", file))
  }
  for (file in sprintf("%02d_%s.R", 1:8, c("qc_diagnostics", "alpha_diversity",
      "beta_diversity", "taxa_composition", "ordination", "shared_taxa",
      "kreport_pavian", "faprotax"))) source(file.path(script_dir, file))

cli_opts <- parse_cli_args()
cfg <- tryCatch(load_config(cli_opts$config, cli_opts = cli_opts),
                error = function(e) fatal("Configuration error", e))
cfg$pipeline_root <- repo_root

module_registry <- list(
  qc = run_qc, alpha = run_alpha, beta = run_beta,
  composition = run_taxa_composition, ordination = run_ordination,
  shared = run_shared_taxa, kreport = run_kreport, faprotax = run_faprotax
)
requested_modules <- cfg$cli$modules
invalid_modules <- setdiff(requested_modules, names(module_registry))
if (length(invalid_modules)) fatal("Configuration error", simpleError(sprintf(
  "Unknown module(s) requested: %s", paste(invalid_modules, collapse = ", "))))
if (isTRUE(cfg$krona$enabled) && !("kreport" %in% requested_modules)) {
  fatal("Configuration error", simpleError("Krona export requires the 'kreport' module."))
}

# All release gates run before output mutation or large input parsing.
tryCatch(validate_output_root(cfg, repo_root), error = function(e) fatal("Preflight error", e))
source_info <- source_provenance(repo_root)
if (isTRUE(source_info$git_dirty) && !isTRUE(cfg$cli$allow_dirty)) {
  fatal("Preflight error", simpleError(
    "E_SOURCE_DIRTY: maintained tracked source files differ from HEAD; use --allow-dirty for development only"))
}
packages <- unique(c(RUNTIME_PACKAGES, get_module_packages(requested_modules)))
  lock_info <- read_lock_status(repo_root, packages)
  if (!identical(lock_info$lock_status, "synchronized") && !isTRUE(cfg$cli$allow_unlocked)) {
    discrepancies <- c(unlist(lock_info$r_discrepancies), unlist(lock_info$package_discrepancies),
                       unlist(lock_info$library_discrepancies))
  fatal("Preflight error", simpleError(sprintf("E_LOCK_MISMATCH: %s",
    paste(discrepancies, collapse = "; "))))
}
if (!identical(lock_info$lock_status, "synchronized")) {
  cat(sprintf("[WARNING] UNLOCKED DEVELOPMENT RUN: lock status is %s.\n", lock_info$lock_status),
      file = stderr())
}
check_module_dependencies(requested_modules)

taxonomy_cache_path <- cfg$taxonomy$cache
taxonomy_lock <- NULL
taxonomy_journal_path <- get_taxonomy_journal_path(taxonomy_cache_path)
if (isTRUE(cfg$cli$validate_only) && file.exists(taxonomy_journal_path)) {
  fatal("Taxonomy recovery error", simpleError(sprintf(
    "E_TAXONOMY_RECOVERY_REQUIRED: pending taxonomy transaction at '%s'; validation cannot mutate it.",
    taxonomy_journal_path)))
}
if (!isTRUE(cfg$cli$validate_only) &&
    (identical(cfg$taxonomy$network_mode, "refresh") || file.exists(taxonomy_journal_path))) {
  taxonomy_lock <- tryCatch(acquire_taxonomy_lock(taxonomy_cache_path),
                            error = function(e) fatal("Taxonomy lock error", e))
  on.exit(release_taxonomy_lock(taxonomy_lock), add = TRUE)
  tryCatch(recover_taxonomy_journal(taxonomy_cache_path),
           error = function(e) fatal("Taxonomy recovery error", e))
  if (!identical(cfg$taxonomy$network_mode, "refresh")) {
    release_taxonomy_lock(taxonomy_lock)
    taxonomy_lock <- NULL
  }
}

  initial_inventory <- tryCatch(inventory_inputs(cfg), error = function(e) fatal("Preflight error", e))
  context <- tryCatch(build_context(cfg), error = function(e) fatal("Input validation error", e))
  tryCatch(validate_output_root(cfg, repo_root,
    extra_paths = unlist(context$bamstats, use.names = FALSE)),
    error = function(e) fatal("Preflight error", e))
  full_inventory <- tryCatch(inventory_inputs(cfg, context$bamstats),
                           error = function(e) fatal("Preflight error", e))
tryCatch(assert_inputs_unchanged(initial_inventory), error = function(e) fatal("Preflight error", e))
context$input_inventory <- full_inventory
preflight_warnings <- tryCatch(run_module_preflight(context, requested_modules),
                              error = function(e) fatal("Module preflight error", e))
tryCatch(assert_inputs_unchanged(full_inventory), error = function(e) fatal("Preflight error", e))

if (isTRUE(cfg$cli$validate_only)) {
  summary <- list(
    status = "validated", pipeline_version = pipeline_version, mode = context$mode,
    samples = unname(context$samples), requested_modules = unname(requested_modules),
    total_reads = sum(context$sample_stats$TotalReads), lock_status = lock_info$lock_status,
    git_dirty = source_info$git_dirty, input_fingerprints = full_inventory,
    warnings = json_array(preflight_warnings)
  )
  cat(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE, null = "null"), "\n")
  if (length(preflight_warnings) > 0L) {
    for (pw in preflight_warnings) {
      cat(sprintf("[WARNING] PREFLIGHT: %s\n", pw), file = stderr())
    }
  }
  cat("Validation check PASSED. Zero filesystem mutations performed.\n")
  return(invisible(0L))
}

final_root <- normalizePath(cfg$output$base_dir, winslash = "/", mustWork = FALSE)

# Acquire exclusive output lock before reading/staging prior output
output_lock <- tryCatch(acquire_output_lock(final_root),
                        error = function(e) fatal("Output lock error", e))
on.exit({
  release_output_lock(output_lock)
}, add = TRUE)

# Recover publication journal if a previous run crashed mid-publication
tryCatch(recover_publication_journal(final_root),
         error = function(e) fatal("Output recovery error", e))

prior_manifest <- tryCatch(validate_prior_output(final_root, cfg$cli$overwrite),
                           error = function(e) fatal("Output validation error", e))

taxonomy_cache_state <- "unchanged" # unchanged, candidate_committed, restored, published
taxonomy_cache_backup_file <- NULL
taxonomy_cache_original_sha256 <- full_inventory$taxonomy_cache$sha256
taxonomy_cache_original_size <- full_inventory$taxonomy_cache$size_bytes
taxonomy_cache_committed_sha256 <- NULL

if (identical(cfg$taxonomy$network_mode, "refresh") && !is.null(taxonomy_cache_path) && file.exists(taxonomy_cache_path)) {
  taxonomy_cache_backup_file <- tempfile(pattern = ".wf16s_tax_backup_",
                                         tmpdir = dirname(taxonomy_cache_path))
  if (!file.copy(taxonomy_cache_path, taxonomy_cache_backup_file, overwrite = TRUE)) {
    fatal("Taxonomy backup error", simpleError(sprintf("Could not create backup of '%s'.", taxonomy_cache_path)))
  }
}

# Register cache cleanup immediately after backup creation so an error while
# preparing the run stage cannot strand unjournaled recovery material.
on.exit({
  if (identical(taxonomy_cache_state, "candidate_committed")) {
    tryCatch(restore_taxonomy_cache(), error = function(e) {
      warning(sprintf("Taxonomy source-cache rollback failed: %s", conditionMessage(e)),
              call. = FALSE)
    })
  }
  if (!is.null(taxonomy_cache_backup_file) && file.exists(taxonomy_cache_backup_file)) {
    tryCatch(cleanup_taxonomy_journal(taxonomy_cache_path, taxonomy_cache_backup_file),
             error = function(e) warning(conditionMessage(e), call. = FALSE))
  }
}, add = TRUE)

update_taxonomy_provenance <- function(root, restored = FALSE) {
  if (is.null(root)) return(invisible(FALSE))
  provenance_path <- file.path(root, "07_Kreport", "taxonomy_provenance.json")
  if (!file.exists(provenance_path)) return(invisible(FALSE))
  provenance <- jsonlite::fromJSON(provenance_path, simplifyVector = FALSE)
  provenance$source_cache_updated <- !isTRUE(restored)
  provenance$source_cache_commit_deferred <- FALSE
  provenance$source_cache_rolled_back <- isTRUE(restored)
  provenance$source_cache_sha256_after <- if (isTRUE(restored)) {
    taxonomy_cache_original_sha256
  } else {
    taxonomy_cache_committed_sha256
  }
  provenance$source_cache_sha256_committed <- provenance$source_cache_sha256_after
  atomic_write_json(provenance, provenance_path)
  invisible(TRUE)
}

restore_taxonomy_cache <- function() {
  if (!identical(cfg$taxonomy$network_mode, "refresh") || is.null(taxonomy_cache_backup_file) || !file.exists(taxonomy_cache_backup_file)) {
    return(invisible(TRUE))
  }
  if (file.exists(taxonomy_cache_path)) {
    cur_info <- file.info(taxonomy_cache_path)
    if (identical(as.numeric(cur_info$size), as.numeric(taxonomy_cache_original_size))) {
      cur_hash <- compute_file_hash(taxonomy_cache_path)
      if (identical(cur_hash, taxonomy_cache_original_sha256)) {
        taxonomy_cache_state <<- "restored"
        update_taxonomy_provenance(current_module_stage, restored = TRUE)
        update_taxonomy_provenance(stage, restored = TRUE)
        cleanup_taxonomy_journal(taxonomy_cache_path, taxonomy_cache_backup_file)
        taxonomy_cache_backup_file <<- NULL
        return(invisible(TRUE))
      }
    }
  }
  current_hash <- compute_file_hash(taxonomy_cache_path)
  if (is.null(taxonomy_cache_committed_sha256) ||
      !identical(current_hash, taxonomy_cache_committed_sha256)) {
    stop(sprintf(
      "E_TAXONOMY_CACHE_CHANGED: refusing rollback because taxonomy cache '%s' no longer matches this run's committed candidate.",
      taxonomy_cache_path), call. = FALSE)
  }
  atomic_replace(taxonomy_cache_path, function(temp) {
    if (!file.copy(taxonomy_cache_backup_file, temp, overwrite = TRUE)) {
      stop(sprintf("Failed to copy taxonomy backup from '%s'.", taxonomy_cache_backup_file), call. = FALSE)
    }
  }, expected_sha256 = taxonomy_cache_original_sha256)
  restored_hash <- compute_file_hash(taxonomy_cache_path)
  if (!identical(restored_hash, taxonomy_cache_original_sha256)) {
    stop("Taxonomy cache restoration failed SHA-256 verification.", call. = FALSE)
  }
  taxonomy_cache_state <<- "restored"
  update_taxonomy_provenance(current_module_stage, restored = TRUE)
  update_taxonomy_provenance(stage, restored = TRUE)
  cleanup_taxonomy_journal(taxonomy_cache_path, taxonomy_cache_backup_file)
  taxonomy_cache_backup_file <<- NULL
  invisible(TRUE)
}

taxonomy_cache_expected_hash <- function() {
  if (!identical(cfg$taxonomy$network_mode, "refresh")) return(NULL)
  if (identical(taxonomy_cache_state, "candidate_committed")) {
    taxonomy_cache_committed_sha256
  } else {
    full_inventory$taxonomy_cache$sha256
  }
}

stage <- prepare_run_staging(final_root)
stage_active <- TRUE
current_module_stage <- NULL
on.exit({
  if (!is.null(current_module_stage) && dir.exists(current_module_stage)) {
    cleanup_module_staging(current_module_stage)
  }
  if (stage_active && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
}, add = TRUE)

preserved_unowned_outputs <- tryCatch(preserve_unowned_outputs(final_root, stage, prior_manifest),
                                      error = function(e) fatal("Output staging error", e))

rebase_output <- function(cfg, root) {
  cfg$output$base_dir <- root
  cfg$output$dirs <- list(
    qc = file.path(root, "01_QC"), alpha = file.path(root, "02_Alpha_Diversity"),
    beta = file.path(root, "03_Beta_Diversity"), composition = file.path(root, "04_Taxa_Composition"),
    ordination = file.path(root, "05_Ordination"), shared_taxa = file.path(root, "06_Shared_Taxa"),
    kreport = file.path(root, "07_Kreport"), faprotax = file.path(root, "08_FAPROTAX")
  )
  cfg$output$manifest_file <- file.path(root, "run_manifest.json")
  cfg$output$resolved_config_file <- file.path(root, "resolved_config.yml")
  cfg$output$session_info_file <- file.path(root, "session_info.txt")
  cfg
}
cfg <- rebase_output(cfg, stage)
context$config <- cfg

cat(sprintf("ONT wf-16s post-processing %s | %s | %d sample(s)\n",
            pipeline_version, context$mode, length(context$samples)))
start_time <- Sys.time()
module_results <- stats::setNames(lapply(names(module_registry), function(name) {
  new_not_run_module_record(sprintf("Module '%s' was not requested.", name))
}), names(module_registry))
any_failed <- FALSE

for (module_name in requested_modules) {
  cat(sprintf(">>> Executing module [%s]...\n", module_name))
  module_stage <- tryCatch(prepare_module_staging(stage, module_name),
                           error = function(e) fatal(sprintf("Module [%s] staging preparation failed", module_name), e))
  current_module_stage <- module_stage

  module_cfg <- rebase_output(cfg, module_stage)
  module_context <- context
  module_context$config <- module_cfg

  started <- Sys.time()
  captured_warnings <- character(0)
  result <- tryCatch(withCallingHandlers(module_registry[[module_name]](module_context), warning = function(w) {
    captured_warnings <<- c(captured_warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(status = "failed", error = conditionMessage(e), outputs = character(0)))

  result <- tryCatch(validate_module_result(result, module_name, module_stage), error = function(e) {
    list(status = "failed", error = paste("Module result contract violation:", conditionMessage(e)),
         outputs = character(0))
  })

  injected_failure <- Sys.getenv("WF16S_INJECT_MODULE_FAILURE", unset = "")
  test_mode <- Sys.getenv("WF16S_TEST_MODE", unset = "")
  if (nzchar(injected_failure) && !identical(test_mode, "1")) {
    fatal("Execution error", simpleError(
      "WF16S_INJECT_MODULE_FAILURE is only permitted when WF16S_TEST_MODE=1."
    ))
  }
  if (identical(test_mode, "1") && identical(injected_failure, module_name)) {
    result <- list(status = "failed", error = sprintf(
      "Injected failure after module '%s' wrote its outputs.", module_name),
      outputs = character(0))
  }
  ended <- Sys.time()
  result$start_time <- utc_timestamp(started)
  result$end_time <- utc_timestamp(ended)
  result$duration_seconds <- as.numeric(difftime(ended, started, units = "secs"))
  returned_warnings <- result$warnings %||% character(0)
  if (!is.character(returned_warnings) || anyNA(returned_warnings)) {
    result$status <- "failed"
    result$error <- "Module returned invalid warnings: must be a character vector with no NA values."
    result$outputs <- character(0)
    result$warnings <- captured_warnings
  } else {
    result$warnings <- unique(c(returned_warnings, captured_warnings))
  }

  if (identical(module_name, "kreport") && identical(result$status, "completed") &&
      identical(cfg$taxonomy$network_mode, "refresh")) {
    candidate_cache <- file.path(module_stage, "07_Kreport", "resolved_taxonomy_cache.json")
    if (!file.exists(candidate_cache)) {
      fatal("Taxonomy commit error", simpleError("Resolver candidate cache is missing."))
    }
    taxonomy_cache_committed_sha256 <- compute_file_hash(candidate_cache)
    current_cache_sha256 <- compute_file_hash(taxonomy_cache_path)
    if (!identical(current_cache_sha256, taxonomy_cache_original_sha256)) {
      fatal("Taxonomy commit error", simpleError(sprintf(
        "E_TAXONOMY_CACHE_CHANGED: taxonomy cache '%s' changed before deferred commit.",
        taxonomy_cache_path)))
    }
    tryCatch(write_taxonomy_journal(
      taxonomy_cache_path, taxonomy_cache_backup_file, final_root,
      taxonomy_cache_original_sha256, taxonomy_cache_committed_sha256, "prepared"
    ), error = function(e) fatal("Taxonomy journal error", e))
    tryCatch(atomic_replace(taxonomy_cache_path, function(temp) {
      if (!file.copy(candidate_cache, temp, overwrite = TRUE)) {
        stop("Could not copy deferred taxonomy candidate.", call. = FALSE)
      }
    }, expected_sha256 = taxonomy_cache_committed_sha256),
    error = function(e) fatal("Taxonomy commit error", e))
    taxonomy_cache_state <- "candidate_committed"
    tryCatch(write_taxonomy_journal(
      taxonomy_cache_path, taxonomy_cache_backup_file, final_root,
      taxonomy_cache_original_sha256, taxonomy_cache_committed_sha256,
      "candidate_committed"
    ), error = function(e) fatal("Taxonomy journal error", e))
    tryCatch(update_taxonomy_provenance(module_stage, restored = FALSE),
             error = function(e) fatal("Taxonomy provenance update error", e))
  }

  input_check <- tryCatch({
    assert_inputs_unchanged(full_inventory,
      allow_taxonomy_cache_change = identical(taxonomy_cache_state, "candidate_committed"),
      taxonomy_cache_expected_sha256 = taxonomy_cache_expected_hash())
    NULL
  }, error = function(e) e)
  if (!is.null(input_check)) {
    result$status <- "failed"
    result$error <- conditionMessage(input_check)
    result$outputs <- character(0)
  }

  if (identical(result$status, "failed")) {
    cleanup_module_staging(module_stage)
    current_module_stage <- NULL
    module_results[[module_name]] <- result
    any_failed <- TRUE
    cat(sprintf("ERROR in module [%s]: %s\n", module_name, result$error), file = stderr())
    if (!isTRUE(cfg$cli$keep_going)) {
      later <- requested_modules[match(module_name, requested_modules):length(requested_modules)]
      later <- setdiff(later, module_name)
      for (name in later) module_results[[name]] <- new_not_run_module_record(sprintf(
        "Pipeline stopped after failure in module '%s'.", module_name))
      break
    }
  } else {
    tryCatch({
      publish_module_staging(module_stage, stage, result$outputs)
    }, error = function(e) fatal(sprintf("Publishing module [%s] outputs to stage failed", module_name), e))
    current_module_stage <- NULL
    to_stage <- function(path) sub(paste0("^", gsub("([][{}()+*^$|\\?.])", "\\\\\\1", module_stage)), stage, path)
    result$outputs <- to_stage(result$outputs)
    module_results[[module_name]] <- result
  }
}

if (any_failed && identical(cfg$taxonomy$network_mode, "refresh")) {
  tryCatch({
    restore_taxonomy_cache()
  }, error = function(e) fatal("Taxonomy source-cache rollback failed", e))
}

end_time <- Sys.time()
lock_locations <- lock_info$package_locations %||% list()
if (length(lock_locations)) {
  deps <- stats::setNames(vapply(lock_locations, function(record) as.character(record$version), character(1)),
                          vapply(lock_locations, function(record) as.character(record$package), character(1)))
} else {
  deps <- get_dependency_versions(packages)
}
session_lines <- c(
  "=== System & Interpreter ===", sprintf("R version: %s", R.version.string),
  sprintf("Platform:  %s", R.version$platform),
  sprintf("Run time:  %s to %s", utc_timestamp(start_time), utc_timestamp(end_time)), "",
  "=== Package Versions ===", sprintf("  %-15s: %s", names(deps), deps), "",
  "=== Full sessionInfo() ===", capture.output(print(sessionInfo()))
)
atomic_write_lines(session_lines, cfg$output$session_info_file)
manifest_cfg <- cfg
manifest_cfg$output <- rebase_output(cfg, final_root)$output
atomic_write_yaml(manifest_cfg, cfg$output$resolved_config_file)

to_final <- function(path) sub(paste0("^", gsub("([][{}()+*^$|\\?.])", "\\\\\\1", stage)),
                               final_root, path)
for (name in names(module_results)) {
  module_results[[name]]$outputs <- to_final(module_results[[name]]$outputs %||% character(0))
}

taxonomy_provenance_path <- file.path(stage, "07_Kreport", "taxonomy_provenance.json")
taxonomy_provenance <- if (file.exists(taxonomy_provenance_path))
  jsonlite::fromJSON(taxonomy_provenance_path, simplifyVector = FALSE) else list()
unresolved_path <- file.path(stage, "07_Kreport", "unresolved_taxids.tsv")
unresolved_count <- if (file.exists(unresolved_path)) max(0L, length(readLines(unresolved_path)) - 1L) else NA_integer_
python <- tryCatch(find_python(), error = function(e) NA_character_)
python_version <- if (!is.na(python)) trimws(paste(unlist(processx::run(python, "--version",
  error_on_status = FALSE)[c("stdout", "stderr")]), collapse = " ")) else NA_character_

# Map relative paths and producer modules
module_output_relpaths <- list()
for (name in names(module_results)) {
  paths <- module_results[[name]]$outputs %||% character(0)
  for (p in paths) {
    rel <- substring(p, nchar(final_root) + 2L)
    module_output_relpaths[[rel]] <- name
  }
}

owned <- sort(unique(c("resolved_config.yml", "session_info.txt", names(module_output_relpaths))))
owned_with_manifest <- sort(c(owned, "run_manifest.json"))

# Verify physical stage census before building artifacts array and manifest
tryCatch(verify_physical_file_census(stage, owned_with_manifest, preserved_unowned_outputs),
         error = function(e) fatal("Physical stage census error", e))

# Build artifacts array for all owned files excluding run_manifest.json
artifacts_list <- lapply(sort(owned), function(rel_path) {
  phys_path <- file.path(stage, rel_path)
  if (!file.exists(phys_path)) {
    fatal("Manifest artifact error", simpleError(sprintf("Owned output '%s' missing from stage.", rel_path)))
  }
  producer <- module_output_relpaths[[rel_path]] %||% NULL
  list(
    relative_path = rel_path,
    size_bytes = as.numeric(file.info(phys_path)$size),
    sha256 = compute_file_hash(phys_path),
    producer_module = producer
  )
})

manifest <- list(
  pipeline = "ont-wf16s-postprocess", pipeline_version = pipeline_version,
  git_commit = source_info$git_commit, git_dirty = source_info$git_dirty,
  source_digest_sha256 = source_info$source_digest_sha256,
  schema_version = 2L, schema_revision = 2L, config_schema_version = cfg$schema_version,
  run_status = if (any_failed) "failed" else "completed",
  start_time = utc_timestamp(start_time), end_time = utc_timestamp(end_time),
  duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
  project_name = cfg$project_name, mode = context$mode, seed = cfg$seed,
  samples = json_array(context$samples), config_file = cfg$config_file,
  output_root = final_root, command = json_array(commandArgs(trailingOnly = FALSE)),
  cli = utils::modifyList(cfg$cli, list(modules = json_array(requested_modules))),
  inputs = full_inventory, upstream_contract = context$upstream_contract,
  modules = lapply(module_results, manifest_module_record),
  owned_outputs = json_array(owned_with_manifest),
  preserved_unowned_outputs = json_array(sort(preserved_unowned_outputs)),
  artifacts = json_array(artifacts_list),
  warnings = json_array(unique(c(context$warnings, if (!lock_info$locked) "Unlocked development run" else character(0),
    if (isTRUE(source_info$git_dirty)) "Dirty-source development run" else character(0),
    unlist(lapply(module_results, function(x) x$warnings %||% character(0)), use.names = FALSE)))),
  taxonomy = list(network_mode = cfg$taxonomy$network_mode,
    unresolved_policy = cfg$taxonomy$unresolved_policy, unresolved_count = unresolved_count,
    conflicts_count = taxonomy_provenance$conflicts_count %||% NA_integer_,
    resolution_source_counts = taxonomy_provenance$resolution_source_counts %||% NULL),
  interpreter = list(r = R.version.string, platform = R.version$platform, python = python_version),
  environment = lock_info,
  package_versions = json_array(lapply(names(deps), function(package) list(
    package = package, version = deps[[package]]
  )))
)

tryCatch(write_manifest_v2(manifest, cfg$output$manifest_file, physical_root = stage),
         error = function(e) fatal("Manifest publication failed", e))
tryCatch(assert_inputs_unchanged(full_inventory,
  allow_taxonomy_cache_change = identical(taxonomy_cache_state, "candidate_committed"),
  taxonomy_cache_expected_sha256 = taxonomy_cache_expected_hash(),
  allow_taxonomy_mtime_change = identical(taxonomy_cache_state, "restored")),
  error = function(e) fatal("Pre-publication input check failed", e))

if (any_failed && !is.null(prior_manifest)) {
  unlink(stage, recursive = TRUE, force = TRUE)
  stage_active <- FALSE
  stop("[FATAL] Transaction aborted; the previous completed output was preserved.", call. = FALSE)
}
tryCatch(publish_staged_run(stage, final_root),
         error = function(e) fatal("Output publication failed", e))
stage_active <- FALSE
if (identical(taxonomy_cache_state, "candidate_committed")) {
  taxonomy_cache_state <- "published"
  tryCatch({
    cleanup_taxonomy_journal(taxonomy_cache_path, taxonomy_cache_backup_file)
    taxonomy_cache_backup_file <- NULL
  }, error = function(e) fatal("Taxonomy transaction cleanup failed", e))
}
if (any_failed) {
  stop(sprintf("[FATAL] Failed run manifest published: %s", file.path(final_root, "run_manifest.json")),
       call. = FALSE)
}
cat(sprintf("Pipeline completed transactionally. Manifest: %s\n", file.path(final_root, "run_manifest.json")))
invisible(0L)
}

status <- tryCatch(main(), error = function(error) {
  cat(sprintf("%s\n", conditionMessage(error)), file = stderr())
  1L
})
quit(status = status, save = "no")
