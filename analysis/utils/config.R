# =============================================================================
# Configuration Loader and Path Resolver
# =============================================================================

suppressMessages(library(yaml))

is_absolute_path <- function(p) {
  if (is.null(p) || is.na(p) || length(p) == 0 || !nzchar(p)) return(FALSE)
  grepl("^(/|\\\\|[A-Za-z]:[/\\])", p)
}

resolve_path <- function(p, base_dir) {
  if (is.null(p) || is.na(p) || length(p) == 0 || !nzchar(p)) return(NULL)
  if (is_absolute_path(p)) {
    normalizePath(p, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(base_dir, p), winslash = "/", mustWork = FALSE)
  }
}

get_default_config <- function() {
  list(
    schema_version = 1L,
    project_name = "ONT_wf16s_Postprocess",
    mode = "auto",
    seed = 42L,
    input = list(
      abundance_table = "output_AAy/abundance_table_species.tsv",
      metadata = NULL,
      params_json = "output_AAy/params.json",
      tax_column = "tax",
      aggregate_columns = c("total"),
      include_samples = NULL,
      assignments = NULL
    ),
    output = list(
      base_dir = "output"
    ),
    qc = list(
      display_min_length = 1200L,
      display_max_length = 1800L,
      target_min_length = NULL,
      target_max_length = NULL
    ),
    alpha = list(
      rarefaction_points = 25L,
      resample_depth = 50000L,
      resample_fraction_cap = 0.90,
      resample_iterations = 100L
    ),
    composition = list(
      top_n_taxa = 15L,
      heatmap_rank = "genus",
      heatmap_transform = "log10_relative"
    ),
    beta = list(
      distances = c("bray", "jaccard"),
      permutations = 999L,
      strata_column = NULL,
      resampling = list(
        enabled = FALSE,
        iterations = 100L,
        depth_fraction_of_minimum = 0.75
      )
    ),
    shared_taxa = list(
      rank = "species",
      minimum_count = 1L,
      group_prevalence = 0.5
    ),
    taxonomy = list(
      cache = "output_AAy/taxonomy_cache.json",
      network_mode = "cache_only",
      unresolved_policy = "warn",
      email_env = "NCBI_EMAIL",
      api_key_env = "NCBI_API_KEY"
    )
  )
}

merge_config <- function(default_cfg, user_cfg) {
  merged <- default_cfg
  for (key in names(user_cfg)) {
    if (is.list(user_cfg[[key]]) && is.list(merged[[key]])) {
      merged[[key]] <- merge_config(merged[[key]], user_cfg[[key]])
    } else {
      merged[[key]] <- user_cfg[[key]]
    }
  }
  merged
}

load_config <- function(config_path = "config.yml", cli_opts = list()) {
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found: '%s'", config_path), call. = FALSE)
  }
  
  config_file_abs <- normalizePath(config_path, winslash = "/", mustWork = TRUE)
  config_dir <- dirname(config_file_abs)
  
  raw_yaml <- yaml::read_yaml(config_file_abs)
  default_cfg <- get_default_config()
  cfg <- merge_config(default_cfg, raw_yaml)
  
  # CLI overrides
  if (!is.null(cli_opts$output_dir) && nzchar(cli_opts$output_dir)) {
    cfg$output$base_dir <- cli_opts$output_dir
  } else if (!is.null(cli_opts[["output-dir"]]) && nzchar(cli_opts[["output-dir"]])) {
    cfg$output$base_dir <- cli_opts[["output-dir"]]
  }
  
  if (isTRUE(cli_opts$refresh_taxonomy) || isTRUE(cli_opts[["refresh-taxonomy"]])) {
    cfg$taxonomy$network_mode <- "refresh"
  }
  
  cfg$cli <- list(
    validate_only = isTRUE(cli_opts$validate_only) || isTRUE(cli_opts[["validate-only"]]),
    keep_going = isTRUE(cli_opts$keep_going) || isTRUE(cli_opts[["keep-going"]]),
    overwrite = isTRUE(cli_opts$overwrite),
    modules = if (!is.null(cli_opts$modules) && nzchar(cli_opts$modules)) {
      strsplit(cli_opts$modules, "[, ]+")[[1]]
    } else {
      c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport")
    }
  )
  
  # Path resolution against config_dir
  cfg$input$abundance_table <- resolve_path(cfg$input$abundance_table, config_dir)
  cfg$input$metadata <- resolve_path(cfg$input$metadata, config_dir)
  cfg$input$params_json <- resolve_path(cfg$input$params_json, config_dir)
  
  if (!is.null(cfg$input$assignments) && is.list(cfg$input$assignments)) {
    for (s in names(cfg$input$assignments)) {
      cfg$input$assignments[[s]] <- resolve_path(cfg$input$assignments[[s]], config_dir)
    }
  }
  
  cfg$taxonomy$cache <- resolve_path(cfg$taxonomy$cache, config_dir)
  
  # Resolve base output directory
  cfg$output$base_dir <- resolve_path(cfg$output$base_dir, config_dir)
  
  # Derive all module output directories
  base_out <- cfg$output$base_dir
  cfg$output$dirs <- list(
    qc = file.path(base_out, "01_QC"),
    alpha = file.path(base_out, "02_Alpha_Diversity"),
    beta = file.path(base_out, "03_Beta_Diversity"),
    composition = file.path(base_out, "04_Taxa_Composition"),
    ordination = file.path(base_out, "05_Ordination"),
    shared_taxa = file.path(base_out, "06_Shared_Taxa"),
    kreport = file.path(base_out, "07_Kreport")
  )
  
  cfg$output$manifest_file <- file.path(base_out, "run_manifest.json")
  cfg$output$resolved_config_file <- file.path(base_out, "resolved_config.yml")
  cfg$output$session_info_file <- file.path(base_out, "session_info.txt")
  
  cfg$config_file <- config_file_abs
  cfg$config_dir <- config_dir
  
  cfg
}
