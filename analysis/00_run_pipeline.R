#!/usr/bin/env Rscript
# =============================================================================
# ONT wf-16s Amplicon Post-Processing Pipeline Runner
# =============================================================================

# Determine script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(getwd(), "analysis"), winslash = "/", mustWork = FALSE)
  }
}

script_dir <- get_script_dir()
repo_root <- normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE)
version_file <- file.path(repo_root, "VERSION")
pipeline_version <- trimws(readLines(version_file, n = 1L, warn = FALSE))
if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", pipeline_version)) {
  stop("VERSION must contain a semantic version such as 0.1.0", call. = FALSE)
}

# Bootstrap dependency checking before sourcing files that attach packages.
source(file.path(script_dir, "utils", "dependencies.R"))
check_dependencies()

# Source all remaining utilities
source(file.path(script_dir, "utils", "cli.R"))
source(file.path(script_dir, "utils", "config.R"))
source(file.path(script_dir, "utils", "io.R"))
source(file.path(script_dir, "utils", "metrics.R"))
source(file.path(script_dir, "utils", "plotting.R"))
source(file.path(script_dir, "utils", "kreport.R"))

# Source all analysis modules
source(file.path(script_dir, "01_qc_diagnostics.R"))
source(file.path(script_dir, "02_alpha_diversity.R"))
source(file.path(script_dir, "03_beta_diversity.R"))
source(file.path(script_dir, "04_taxa_composition.R"))
source(file.path(script_dir, "05_ordination.R"))
source(file.path(script_dir, "06_shared_taxa.R"))
source(file.path(script_dir, "07_kreport_pavian.R"))

start_time <- Sys.time()

# 1. Parse CLI options
cli_opts <- parse_cli_args()

# 2. Load and resolve configuration
cfg <- tryCatch({
  load_config(cli_opts$config, cli_opts = cli_opts)
}, error = function(e) {
  cat(sprintf("[FATAL] Configuration error: %s\n", e$message), file = stderr())
  quit(status = 1)
})
cfg$pipeline_root <- repo_root

# 3. Build and validate shared context
context <- tryCatch({
  build_context(cfg)
}, error = function(e) {
  cat(sprintf("[FATAL] Input validation error: %s\n", e$message), file = stderr())
  quit(status = 1)
})

# Module registry and request validation must happen before any output mutation.
module_registry <- list(
  qc          = run_qc,
  alpha       = run_alpha,
  beta        = run_beta,
  composition = run_taxa_composition,
  ordination  = run_ordination,
  shared      = run_shared_taxa,
  kreport     = run_kreport
)

requested_modules <- cfg$cli$modules
invalid_modules <- setdiff(requested_modules, names(module_registry))
if (length(invalid_modules) > 0) {
  cat(sprintf("[FATAL] Unknown module(s) requested: %s\nAvailable: %s\n",
              paste(invalid_modules, collapse = ", "),
              paste(names(module_registry), collapse = ", ")), file = stderr())
  quit(status = 1)
}

# 4. Handle --validate-only
if (cfg$cli$validate_only) {
  cat("=== ONT wf-16s Pipeline Validation Check ===\n")
  cat(sprintf("Config file:      %s\n", cfg$config_file))
  cat(sprintf("Mode:             %s\n", context$mode))
  cat(sprintf("Samples (%d):      %s\n", length(context$samples), paste(context$samples, collapse = ", ")))
  cat(sprintf("Total reads:      %s\n", format(sum(context$sample_stats$TotalReads), big.mark = ",")))
  cat(sprintf("Classified reads: %s\n", format(sum(context$sample_stats$ClassifiedReads), big.mark = ",")))
  cat(sprintf("Abundance SHA256: %s\n", context$file_hashes$abundance_table))
  cat("Validation check PASSED. Zero filesystem mutations performed.\n")
  quit(status = 0)
}

# 5. Overwrite check: protect every known module output, including partial runs
# that failed before a manifest could be written.
known_existing <- c(
  c(cfg$output$manifest_file, cfg$output$resolved_config_file, cfg$output$session_info_file)[
    file.exists(c(cfg$output$manifest_file, cfg$output$resolved_config_file, cfg$output$session_info_file))
  ],
  unlist(lapply(cfg$output$dirs, function(path) {
    if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE) else character(0)
  }), use.names = FALSE)
)
if (!cfg$cli$overwrite && length(known_existing) > 0L) {
  cat(sprintf(
    "[FATAL] Refusing to overwrite %d existing pipeline output(s); first path: '%s'\nUse --overwrite to allow replacement.\n",
    length(known_existing), known_existing[1]
  ), file = stderr())
  quit(status = 1)
}

# Create base output directory
dir.create(cfg$output$base_dir, recursive = TRUE, showWarnings = FALSE)

cat("=============================================================================\n")
cat(sprintf("ONT wf-16s Amplicon Post-Processing Pipeline\n"))
cat(sprintf("Project: %s | Mode: %s | Samples: %d\n",
            cfg$project_name, context$mode, length(context$samples)))
cat(sprintf("Output root: %s\n", cfg$output$base_dir))
cat("=============================================================================\n")

module_results <- list()
any_failed <- FALSE

for (mod_name in requested_modules) {
  cat(sprintf("\n>>> Executing module [%s]...\n", mod_name))
  mod_fn <- module_registry[[mod_name]]
  mod_start <- Sys.time()

  mod_res <- tryCatch({
    mod_fn(context)
  }, error = function(e) {
    cat(sprintf("ERROR in module [%s]: %s\n", mod_name, e$message), file = stderr())
    list(
      status = "failed",
      error = e$message,
      outputs = character(0)
    )
  })

  mod_end <- Sys.time()
  mod_res$start_time <- format(mod_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  mod_res$end_time <- format(mod_end, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  mod_res$duration_seconds <- as.numeric(difftime(mod_end, mod_start, units = "secs"))
  module_results[[mod_name]] <- mod_res

  if (mod_res$status == "failed") {
    any_failed <- TRUE
    if (!cfg$cli$keep_going) {
      cat(sprintf("\n[FATAL] Pipeline stopped due to failure in module [%s]. (Use --keep-going to continue past errors)\n", mod_name), file = stderr())
      break
    }
  } else if (mod_res$status == "skipped") {
    cat(sprintf("Module [%s] SKIPPED: %s\n", mod_name, mod_res$reason %||% "prerequisites not met"))
  } else {
    cat(sprintf("Module [%s] COMPLETED (%d output artifacts produced)\n",
                mod_name, length(mod_res$outputs)))
  }
}

end_time <- Sys.time()
overall_status <- if (any_failed) "failed" else "completed"

deps <- get_dependency_versions()
session_lines <- c(
  "=== System & Interpreter ===",
  sprintf("R version: %s", R.version.string),
  sprintf("Platform:  %s", R.version$platform),
  sprintf("Run time:  %s to %s", start_time, end_time),
  "",
  "=== Package Versions ===",
  sprintf("  %-15s: %s", names(deps), deps),
  "",
  "=== Full sessionInfo() ===",
  capture.output(print(sessionInfo()))
)
writeLines(session_lines, cfg$output$session_info_file)

# 8. Write resolved config
yaml::write_yaml(cfg, cfg$output$resolved_config_file)

# 9. Build and write run manifest
input_meta <- list(
  abundance_table = list(
    path = cfg$input$abundance_table,
    size_bytes = if (file.exists(cfg$input$abundance_table)) file.info(cfg$input$abundance_table)$size else 0,
    mtime = if (file.exists(cfg$input$abundance_table)) as.character(file.info(cfg$input$abundance_table)$mtime) else NA,
    sha256 = context$file_hashes$abundance_table
  ),
  metadata = if (!is.null(cfg$input$metadata) && file.exists(cfg$input$metadata)) list(
    path = cfg$input$metadata,
    size_bytes = file.info(cfg$input$metadata)$size,
    mtime = as.character(file.info(cfg$input$metadata)$mtime),
    sha256 = context$file_hashes$metadata
  ) else NULL,
  params_json = if (!is.null(cfg$input$params_json) && file.exists(cfg$input$params_json)) list(
    path = cfg$input$params_json,
    size_bytes = file.info(cfg$input$params_json)$size,
    mtime = as.character(file.info(cfg$input$params_json)$mtime),
    sha256 = context$file_hashes$params_json
  ) else NULL,
  taxonomy_cache = if (!is.null(cfg$taxonomy$cache) && file.exists(cfg$taxonomy$cache)) list(
    path = cfg$taxonomy$cache,
    size_bytes = file.info(cfg$taxonomy$cache)$size,
    mtime = as.character(file.info(cfg$taxonomy$cache)$mtime),
    sha256 = context$file_hashes$taxonomy_cache
  ) else NULL,
  assignments = if (length(context$assignments) > 0L) {
    lapply(names(context$assignments), function(sample_id) {
      path <- context$assignments[[sample_id]]
      list(
        sample_id = sample_id,
        path = path,
        size_bytes = file.info(path)$size,
        mtime = as.character(file.info(path)$mtime),
        sha256 = context$file_hashes[[paste0("assignment_", sample_id)]]
      )
    })
  } else NULL
)

unresolved_file <- file.path(cfg$output$dirs$kreport, "unresolved_taxids.tsv")
unresolved_count <- if (file.exists(unresolved_file)) {
  max(0L, length(readLines(unresolved_file, warn = FALSE)) - 1L)
} else {
  NA_integer_
}

taxonomy_provenance_file <- file.path(cfg$output$dirs$kreport, "taxonomy_provenance.json")
taxonomy_provenance <- if (file.exists(taxonomy_provenance_file)) {
  jsonlite::fromJSON(taxonomy_provenance_file, simplifyVector = FALSE)
} else {
  NULL
}

python_cmd <- tryCatch(find_python(), error = function(e) NA_character_)
python_version <- if (!is.na(python_cmd)) {
  tryCatch({
    probe <- processx::run(python_cmd, "--version", error_on_status = FALSE)
    trimws(paste(c(probe$stdout, probe$stderr), collapse = " "))
  }, error = function(e) NA_character_)
} else {
  NA_character_
}

git_commit <- tryCatch({
  result <- processx::run(
    "git", c("-c", paste0("safe.directory=", repo_root),
             "-C", repo_root, "rev-parse", "HEAD"),
    error_on_status = FALSE
  )
  if (identical(result$status, 0L)) trimws(result$stdout) else NULL
}, error = function(e) NULL)

manifest <- list(
  pipeline = "ont-wf16s-postprocess",
  pipeline_version = pipeline_version,
  git_commit = git_commit,
  schema_version = cfg$schema_version,
  run_status = overall_status,
  start_time = format(start_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  end_time = format(end_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
  project_name = cfg$project_name,
  mode = context$mode,
  seed = cfg$seed,
  samples = context$samples,
  config_file = cfg$config_file,
  output_root = cfg$output$base_dir,
  command = commandArgs(trailingOnly = FALSE),
  cli = cfg$cli,
  inputs = input_meta,
  modules = module_results,
  warnings = context$warnings,
  taxonomy = list(
    network_mode = cfg$taxonomy$network_mode,
    unresolved_policy = cfg$taxonomy$unresolved_policy,
    unresolved_count = unresolved_count,
    conflicts_count = taxonomy_provenance$conflicts_count %||% NA_integer_,
    resolution_source_counts = taxonomy_provenance$resolution_source_counts %||% NULL
  ),
  interpreter = list(
    r = R.version.string,
    platform = R.version$platform,
    python = python_version
  ),
  package_versions = deps
)

jsonlite::write_json(manifest, cfg$output$manifest_file, pretty = TRUE, auto_unbox = TRUE)

cat("\n=============================================================================\n")
cat(sprintf("Pipeline finished with status: [%s]\n", toupper(overall_status)))
cat(sprintf("Manifest written to: %s\n", cfg$output$manifest_file))
cat(sprintf("Resolved config to: %s\n", cfg$output$resolved_config_file))
cat("=============================================================================\n")

if (any_failed) {
  quit(status = 1)
} else {
  quit(status = 0)
}
