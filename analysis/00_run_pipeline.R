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

# Source all utilities
source(file.path(script_dir, "utils", "dependencies.R"))
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

# 1. Dependency check
check_dependencies()

# 2. Parse CLI options
cli_opts <- parse_cli_args()

# 3. Load and resolve configuration
cfg <- tryCatch({
  load_config(cli_opts$config, cli_opts = cli_opts)
}, error = function(e) {
  cat(sprintf("[FATAL] Configuration error: %s\n", e$message), file = stderr())
  quit(status = 1)
})

# 4. Build and validate shared context
context <- tryCatch({
  build_context(cfg)
}, error = function(e) {
  cat(sprintf("[FATAL] Input validation error: %s\n", e$message), file = stderr())
  quit(status = 1)
})

# 5. Handle --validate-only
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

# 6. Overwrite check
if (!cfg$cli$overwrite && file.exists(cfg$output$manifest_file)) {
  cat(sprintf(
    "[FATAL] Output file already exists: '%s'\nUse --overwrite to allow replacing existing outputs.\n",
    cfg$output$manifest_file
  ), file = stderr())
  quit(status = 1)
}

# Create base output directory
dir.create(cfg$output$base_dir, recursive = TRUE, showWarnings = FALSE)

# Module registry
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

# 7. Write session info
sink(cfg$output$session_info_file)
cat("=== System & Interpreter ===\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform:  %s\n", R.version$platform))
cat(sprintf("Run time:  %s to %s\n\n", start_time, end_time))
cat("=== Package Versions ===\n")
deps <- get_dependency_versions()
for (pkg in names(deps)) {
  cat(sprintf("  %-15s: %s\n", pkg, deps[[pkg]]))
}
cat("\n=== Full sessionInfo() ===\n")
print(sessionInfo())
sink()

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
  ) else NULL
)

manifest <- list(
  pipeline = "ont-wf16s-postprocess",
  schema_version = cfg$schema_version,
  run_status = overall_status,
  start_time = format(start_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  end_time = format(end_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
  project_name = cfg$project_name,
  mode = context$mode,
  seed = cfg$seed,
  samples = context$samples,
  inputs = input_meta,
  modules = module_results,
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
