# =============================================================================
# Configuration Loader and Path Resolver
# =============================================================================

suppressMessages(library(yaml))

`%||%` <- function(x, fallback) {
  if (is.null(x) || length(x) == 0L) fallback else x
}

assert_scalar_number <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE,
                                 lower_open = FALSE, upper_open = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid && integer) valid <- abs(x - round(x)) <= sqrt(.Machine$double.eps)
  if (valid) valid <- if (lower_open) x > lower else x >= lower
  if (valid) valid <- if (upper_open) x < upper else x <= upper
  if (!valid) stop(sprintf("Invalid configuration value '%s'.", name), call. = FALSE)
  invisible(TRUE)
}

assert_nonempty_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(sprintf("'%s' must be one non-empty string.", name), call. = FALSE)
  }
  invisible(TRUE)
}

validate_config <- function(cfg) {
  if (!identical(as.integer(cfg$schema_version), 1L)) {
    stop("Unsupported schema_version; expected 1.", call. = FALSE)
  }
  assert_nonempty_string(cfg$project_name, "project_name")
  if (!is.character(cfg$mode) || length(cfg$mode) != 1L ||
      !cfg$mode %in% c("auto", "single", "cohort")) {
    stop("'mode' must be one of: auto, single, cohort.", call. = FALSE)
  }
  assert_scalar_number(cfg$seed, "seed", lower = 0, integer = TRUE)

  assert_nonempty_string(cfg$input$abundance_table, "input.abundance_table")
  assert_nonempty_string(cfg$input$tax_column, "input.tax_column")
  assert_nonempty_string(cfg$output$base_dir, "output.base_dir")
  if (!is.null(cfg$input$aggregate_columns) &&
      (!is.character(cfg$input$aggregate_columns) || anyNA(cfg$input$aggregate_columns) ||
       any(!nzchar(cfg$input$aggregate_columns)) || anyDuplicated(cfg$input$aggregate_columns))) {
    stop("'input.aggregate_columns' must contain unique, non-empty strings.", call. = FALSE)
  }
  if (cfg$input$tax_column %in% cfg$input$aggregate_columns) {
    stop("'input.tax_column' cannot also be listed in 'input.aggregate_columns'.", call. = FALSE)
  }
  if (!is.null(cfg$input$include_samples) &&
      (!is.character(cfg$input$include_samples) || anyNA(cfg$input$include_samples) ||
       any(!nzchar(cfg$input$include_samples)) || anyDuplicated(cfg$input$include_samples))) {
    stop("'input.include_samples' must contain unique, non-empty sample IDs.", call. = FALSE)
  }
  if (!is.null(cfg$input$assignments)) {
    if (!is.list(cfg$input$assignments) || is.null(names(cfg$input$assignments)) ||
        any(!nzchar(names(cfg$input$assignments))) || anyDuplicated(names(cfg$input$assignments))) {
      stop("'input.assignments' must be a named mapping of unique SampleID -> path.", call. = FALSE)
    }
    valid_paths <- vapply(cfg$input$assignments, function(path) {
      is.character(path) && length(path) == 1L && !is.na(path) && nzchar(trimws(path))
    }, logical(1))
    if (!all(valid_paths)) stop("Every 'input.assignments' value must be one non-empty path.", call. = FALSE)
  }

  assert_scalar_number(cfg$qc$display_min_length, "qc.display_min_length", lower = 1, integer = TRUE)
  assert_scalar_number(cfg$qc$display_max_length, "qc.display_max_length", lower = 1, integer = TRUE)
  if (!is.null(cfg$qc$target_min_length)) {
    assert_scalar_number(cfg$qc$target_min_length, "qc.target_min_length", lower = 1, integer = TRUE)
  }
  if (!is.null(cfg$qc$target_max_length)) {
    assert_scalar_number(cfg$qc$target_max_length, "qc.target_max_length", lower = 1, integer = TRUE)
  }
  if (cfg$qc$display_min_length >= cfg$qc$display_max_length) {
    stop("'qc.display_min_length' must be smaller than 'qc.display_max_length'.", call. = FALSE)
  }
  assert_scalar_number(cfg$alpha$rarefaction_points, "alpha.rarefaction_points", lower = 2, integer = TRUE)
  assert_scalar_number(cfg$alpha$resample_depth, "alpha.resample_depth", lower = 1, integer = TRUE)
  assert_scalar_number(cfg$alpha$resample_fraction_cap, "alpha.resample_fraction_cap",
                       lower = 0, upper = 1, lower_open = TRUE)
  assert_scalar_number(cfg$alpha$resample_iterations, "alpha.resample_iterations", lower = 1, integer = TRUE)
  assert_scalar_number(cfg$composition$top_n_taxa, "composition.top_n_taxa", lower = 1, integer = TRUE)
  assert_nonempty_string(cfg$composition$heatmap_rank, "composition.heatmap_rank")
  assert_nonempty_string(cfg$composition$heatmap_transform, "composition.heatmap_transform")
  if (!cfg$composition$heatmap_rank %in% c("phylum", "class", "order", "family", "genus", "species")) {
    stop("'composition.heatmap_rank' is not a supported rank.", call. = FALSE)
  }
  if (!cfg$composition$heatmap_transform %in% c("log10_relative", "none")) {
    stop("'composition.heatmap_transform' must be 'log10_relative' or 'none'.", call. = FALSE)
  }
  if (!is.character(cfg$beta$distances) || length(cfg$beta$distances) == 0L ||
      any(!cfg$beta$distances %in% c("bray", "jaccard")) || anyDuplicated(cfg$beta$distances)) {
    stop("'beta.distances' must contain unique values drawn from: bray, jaccard.", call. = FALSE)
  }
  assert_scalar_number(cfg$beta$permutations, "beta.permutations", lower = 1, integer = TRUE)
  assert_scalar_number(cfg$beta$minimum_count, "beta.minimum_count", lower = 1, integer = TRUE)
  if (!is.logical(cfg$beta$resampling$enabled) || length(cfg$beta$resampling$enabled) != 1L ||
      is.na(cfg$beta$resampling$enabled)) {
    stop("'beta.resampling.enabled' must be true or false.", call. = FALSE)
  }
  assert_scalar_number(cfg$beta$resampling$iterations, "beta.resampling.iterations", lower = 1, integer = TRUE)
  assert_scalar_number(cfg$beta$resampling$depth_fraction_of_minimum,
                       "beta.resampling.depth_fraction_of_minimum",
                       lower = 0, upper = 1, lower_open = TRUE)
  if (!is.null(cfg$beta$strata_column) &&
      (!is.character(cfg$beta$strata_column) || length(cfg$beta$strata_column) != 1L ||
       !nzchar(cfg$beta$strata_column))) {
    stop("'beta.strata_column' must be null or one non-empty metadata column name.", call. = FALSE)
  }
  assert_nonempty_string(cfg$shared_taxa$rank, "shared_taxa.rank")
  if (!cfg$shared_taxa$rank %in% c("superkingdom", "kingdom", "phylum", "class",
                                   "order", "family", "genus", "species")) {
    stop("'shared_taxa.rank' is not a supported rank.", call. = FALSE)
  }
  assert_scalar_number(cfg$shared_taxa$minimum_count, "shared_taxa.minimum_count", lower = 1, integer = TRUE)
  assert_scalar_number(cfg$shared_taxa$group_prevalence, "shared_taxa.group_prevalence",
                       lower = 0, upper = 1, lower_open = TRUE)
  assert_nonempty_string(cfg$taxonomy$cache, "taxonomy.cache")
  assert_nonempty_string(cfg$taxonomy$network_mode, "taxonomy.network_mode")
  assert_nonempty_string(cfg$taxonomy$unresolved_policy, "taxonomy.unresolved_policy")
  assert_nonempty_string(cfg$taxonomy$email_env, "taxonomy.email_env")
  assert_nonempty_string(cfg$taxonomy$api_key_env, "taxonomy.api_key_env")
  if (!cfg$taxonomy$network_mode %in% c("cache_only", "refresh")) {
    stop("'taxonomy.network_mode' must be 'cache_only' or 'refresh'.", call. = FALSE)
  }
  if (!cfg$taxonomy$unresolved_policy %in% c("warn", "error")) {
    stop("'taxonomy.unresolved_policy' must be 'warn' or 'error'.", call. = FALSE)
  }
  invisible(cfg)
}

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
      minimum_count = 1L,
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

merge_config <- function(default_cfg, user_cfg, path = "") {
  unknown <- setdiff(names(user_cfg), names(default_cfg))
  if (length(unknown) > 0L) {
    qualified <- paste0(path, unknown)
    stop(sprintf("Unknown configuration key(s): %s", paste(qualified, collapse = ", ")), call. = FALSE)
  }
  merged <- default_cfg
  for (key in names(user_cfg)) {
    if (is.list(user_cfg[[key]]) && is.list(merged[[key]])) {
      merged[[key]] <- merge_config(merged[[key]], user_cfg[[key]], paste0(path, key, "."))
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
  cli_refresh <- isTRUE(cli_opts$refresh_taxonomy) ||
    isTRUE(cli_opts[["refresh-taxonomy"]])

  if (identical(cfg$taxonomy$network_mode, "refresh") && !cli_refresh) {
    stop(paste(
      "YAML cannot enable taxonomy refresh.",
      "Keep taxonomy.network_mode: cache_only and pass --refresh-taxonomy explicitly."
    ), call. = FALSE)
  }

  # CLI overrides
  if (!is.null(cli_opts$output_dir) && nzchar(cli_opts$output_dir)) {
    cfg$output$base_dir <- cli_opts$output_dir
  } else if (!is.null(cli_opts[["output-dir"]]) && nzchar(cli_opts[["output-dir"]])) {
    cfg$output$base_dir <- cli_opts[["output-dir"]]
  }

  if (cli_refresh) {
    cfg$taxonomy$network_mode <- "refresh"
  }

  validate_config(cfg)

  cfg$cli <- list(
    validate_only = isTRUE(cli_opts$validate_only) || isTRUE(cli_opts[["validate-only"]]),
    keep_going = isTRUE(cli_opts$keep_going) || isTRUE(cli_opts[["keep-going"]]),
    overwrite = isTRUE(cli_opts$overwrite),
    refresh_taxonomy = cli_refresh,
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
