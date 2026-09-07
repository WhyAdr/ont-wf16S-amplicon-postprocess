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

script_dir <- get_script_dir()
repo_root <- normalizePath(dirname(script_dir), winslash = "/")
options(wf16s.pipeline_root = repo_root)
pipeline_version <- read_pipeline_version(file.path(repo_root, "VERSION"))

source(file.path(script_dir, "utils", "dependencies.R"))
check_dependencies()
for (file in c("cli.R", "config.R", "io.R", "metrics.R", "plotting.R", "manifest.R",
               "kreport.R", "atomic_io.R", "module_result.R", "preflight.R")) {
  source(file.path(script_dir, "utils", file))
}
for (file in sprintf("%02d_%s.R", 1:8, c("qc_diagnostics", "alpha_diversity",
    "beta_diversity", "taxa_composition", "ordination", "shared_taxa",
    "kreport_pavian", "faprotax"))) source(file.path(script_dir, file))

fatal <- function(label, error) {
  if (exists("taxonomy_cache_committed", inherits = TRUE) &&
      exists("restore_taxonomy_cache", inherits = TRUE) &&
      !isTRUE(taxonomy_cache_committed)) {
    tryCatch({
      restore_taxonomy_cache()
      taxonomy_cache_committed <<- TRUE
    }, error = function(restore_error) {
      cat(sprintf("[FATAL] Could not restore taxonomy source cache: %s\n",
                  conditionMessage(restore_error)), file = stderr())
    })
  }
  cat(sprintf("[FATAL] %s: %s\n", label, conditionMessage(error)), file = stderr())
  quit(status = 1L)
}

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
  discrepancies <- c(unlist(lock_info$r_discrepancies), unlist(lock_info$package_discrepancies))
  fatal("Preflight error", simpleError(sprintf("E_LOCK_MISMATCH: %s",
    paste(discrepancies, collapse = "; "))))
}
if (!identical(lock_info$lock_status, "synchronized")) {
  cat(sprintf("[WARNING] UNLOCKED DEVELOPMENT RUN: lock status is %s.\n", lock_info$lock_status),
      file = stderr())
}
check_module_dependencies(requested_modules)

initial_inventory <- tryCatch(inventory_inputs(cfg), error = function(e) fatal("Preflight error", e))
context <- tryCatch(build_context(cfg), error = function(e) fatal("Input validation error", e))
full_inventory <- tryCatch(inventory_inputs(cfg, context$bamstats),
                           error = function(e) fatal("Preflight error", e))
tryCatch(assert_inputs_unchanged(initial_inventory), error = function(e) fatal("Preflight error", e))
context$input_inventory <- full_inventory
tryCatch(run_module_preflight(context, requested_modules),
         error = function(e) fatal("Module preflight error", e))
tryCatch(assert_inputs_unchanged(full_inventory), error = function(e) fatal("Preflight error", e))

if (isTRUE(cfg$cli$validate_only)) {
  summary <- list(
    status = "validated", pipeline_version = pipeline_version, mode = context$mode,
    samples = unname(context$samples), requested_modules = unname(requested_modules),
    total_reads = sum(context$sample_stats$TotalReads), lock_status = lock_info$lock_status,
    git_dirty = source_info$git_dirty, input_fingerprints = full_inventory
  )
  cat(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE, null = "null"), "\n")
  cat("Validation check PASSED. Zero filesystem mutations performed.\n")
  quit(status = 0L)
}

final_root <- normalizePath(cfg$output$base_dir, winslash = "/", mustWork = FALSE)
prior_manifest <- tryCatch(validate_prior_output(final_root, cfg$cli$overwrite),
                           error = function(e) fatal("Output validation error", e))

taxonomy_cache_path <- cfg$taxonomy$cache
taxonomy_cache_original <- if (identical(cfg$taxonomy$network_mode, "refresh")) {
  readBin(taxonomy_cache_path, "raw", n = file.info(taxonomy_cache_path)$size)
} else {
  NULL
}
taxonomy_cache_committed <- !identical(cfg$taxonomy$network_mode, "refresh")
restore_taxonomy_cache <- function() {
  if (!identical(cfg$taxonomy$network_mode, "refresh") || is.null(taxonomy_cache_original)) return(invisible(TRUE))
  atomic_replace(taxonomy_cache_path, function(temp) writeBin(taxonomy_cache_original, temp))
  invisible(TRUE)
}

stage <- prepare_run_staging(final_root)
stage_active <- TRUE
on.exit({
  if (!taxonomy_cache_committed) restore_taxonomy_cache()
  if (stage_active && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
}, add = TRUE)
tryCatch(preserve_unowned_outputs(final_root, stage, prior_manifest),
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
  started <- Sys.time()
  warnings <- character(0)
  result <- tryCatch(withCallingHandlers(module_registry[[module_name]](context), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(status = "failed", error = conditionMessage(e), outputs = character(0)))
  result <- tryCatch(validate_module_result(result, module_name, stage), error = function(e) {
    list(status = "failed", error = paste("Module result contract violation:", conditionMessage(e)),
         outputs = character(0))
  })
  ended <- Sys.time()
  result$start_time <- utc_timestamp(started)
  result$end_time <- utc_timestamp(ended)
  result$duration_seconds <- as.numeric(difftime(ended, started, units = "secs"))
  result$warnings <- unique(warnings)
  module_results[[module_name]] <- result
  tryCatch(assert_inputs_unchanged(full_inventory,
    allow_taxonomy_cache_change = identical(cfg$taxonomy$network_mode, "refresh")),
    error = function(e) {
      result$status <<- "failed"; result$error <<- conditionMessage(e); result$outputs <<- character(0)
      module_results[[module_name]] <<- result
    })
  if (identical(result$status, "failed")) {
    any_failed <- TRUE
    cat(sprintf("ERROR in module [%s]: %s\n", module_name, result$error), file = stderr())
    if (!isTRUE(cfg$cli$keep_going)) {
      later <- requested_modules[match(module_name, requested_modules):length(requested_modules)]
      later <- setdiff(later, module_name)
      for (name in later) module_results[[name]] <- new_not_run_module_record(sprintf(
        "Pipeline stopped after failure in module '%s'.", module_name))
      break
    }
  }
}

if (any_failed && !taxonomy_cache_committed) {
  tryCatch({
    restore_taxonomy_cache()
    taxonomy_cache_committed <- TRUE
  }, error = function(e) fatal("Taxonomy source-cache rollback failed", e))
}

end_time <- Sys.time()
deps <- get_dependency_versions(packages)
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

owned <- unique(c("resolved_config.yml", "session_info.txt", unlist(lapply(module_results, function(x) {
  paths <- x$outputs %||% character(0)
  ifelse(startsWith(tolower(paths), paste0(tolower(final_root), "/")),
         substring(paths, nchar(final_root) + 2L), paths)
}), use.names = FALSE)))
manifest <- list(
  pipeline = "ont-wf16s-postprocess", pipeline_version = pipeline_version,
  git_commit = source_info$git_commit, git_dirty = source_info$git_dirty,
  source_digest_sha256 = source_info$source_digest_sha256,
  schema_version = 2L, schema_revision = 1L, config_schema_version = cfg$schema_version,
  run_status = if (any_failed) "failed" else "completed",
  start_time = utc_timestamp(start_time), end_time = utc_timestamp(end_time),
  duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
  project_name = cfg$project_name, mode = context$mode, seed = cfg$seed,
  samples = json_array(context$samples), config_file = cfg$config_file,
  output_root = final_root, command = json_array(commandArgs(trailingOnly = FALSE)),
  cli = utils::modifyList(cfg$cli, list(modules = json_array(requested_modules))),
  inputs = full_inventory, upstream_contract = context$upstream_contract,
  modules = lapply(module_results, manifest_module_record),
  owned_outputs = json_array(sort(c(owned, "run_manifest.json"))),
  warnings = json_array(unique(c(context$warnings, if (!lock_info$locked) "Unlocked development run" else character(0),
    if (isTRUE(source_info$git_dirty)) "Dirty-source development run" else character(0),
    unlist(lapply(module_results, function(x) x$warnings %||% character(0)), use.names = FALSE)))),
  taxonomy = list(network_mode = cfg$taxonomy$network_mode,
    unresolved_policy = cfg$taxonomy$unresolved_policy, unresolved_count = unresolved_count,
    conflicts_count = taxonomy_provenance$conflicts_count %||% NA_integer_,
    resolution_source_counts = taxonomy_provenance$resolution_source_counts %||% NULL),
  interpreter = list(r = R.version.string, platform = R.version$platform, python = python_version),
  environment = lock_info,
  package_versions = json_array(lapply(names(deps), function(package) list(package = package, version = deps[[package]])))
)
tryCatch(write_manifest_v2(manifest, cfg$output$manifest_file, physical_root = stage),
         error = function(e) fatal("Manifest publication failed", e))
tryCatch(assert_inputs_unchanged(full_inventory,
  allow_taxonomy_cache_change = identical(cfg$taxonomy$network_mode, "refresh")),
  error = function(e) fatal("Pre-publication input check failed", e))
if (any_failed && !is.null(prior_manifest)) {
  unlink(stage, recursive = TRUE, force = TRUE)
  stage_active <- FALSE
  cat("[FATAL] Transaction aborted; the previous completed output was preserved.\n", file = stderr())
  quit(status = 1L)
}
tryCatch(publish_staged_run(stage, final_root),
         error = function(e) fatal("Output publication failed", e))
stage_active <- FALSE
if (!any_failed) taxonomy_cache_committed <- TRUE
if (any_failed) {
  cat(sprintf("[FATAL] Failed run manifest published: %s\n", file.path(final_root, "run_manifest.json")), file = stderr())
  quit(status = 1L)
}
cat(sprintf("Pipeline completed transactionally. Manifest: %s\n", file.path(final_root, "run_manifest.json")))
