# =============================================================================
# Configuration Loader and Path Resolver
# =============================================================================

suppressMessages(library(yaml))

`%||%` <- function(x, fallback) {
  if (is.null(x) || length(x) == 0L) fallback else x
}

SUPPORTED_COMPOSITION_RANKS <- c(
  "phylum", "class", "order", "family", "genus", "species"
)

assert_scalar_number <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE,
                                 lower_open = FALSE, upper_open = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid && integer) valid <- x == floor(x)
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

parse_requested_modules <- function(value) {
  defaults <- c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport")
  if (is.null(value)) return(defaults)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop("'--modules' must contain at least one module name.", call. = FALSE)
  }

  if (!identical(value, trimws(value)) || grepl("(^,|,$|,,|,\\s|\\s,)", value, perl = TRUE)) {
    stop("'--modules' must be a strict comma-separated list without whitespace or empty entries.",
         call. = FALSE)
  }
  modules <- strsplit(value, ",", fixed = TRUE)[[1]]
  if (any(!grepl("^[a-z][a-z0-9_]*$", modules))) {
    stop("'--modules' contains an invalid module name.", call. = FALSE)
  }
  duplicates <- unique(modules[duplicated(modules)])
  if (length(duplicates) > 0L) {
    stop(sprintf("Duplicate module name(s): %s", paste(duplicates, collapse = ", ")),
         call. = FALSE)
  }
  modules
}

validate_config <- function(cfg) {
  if (!is.numeric(cfg$schema_version) || length(cfg$schema_version) != 1L ||
      is.na(cfg$schema_version) || !is.finite(cfg$schema_version) || cfg$schema_version != 1) {
    stop("Unsupported schema_version; expected 1.", call. = FALSE)
  }
  assert_nonempty_string(cfg$project_name, "project_name")
  if (!is.character(cfg$mode) || length(cfg$mode) != 1L ||
      !cfg$mode %in% c("auto", "single", "cohort")) {
    stop("'mode' must be one of: auto, single, cohort.", call. = FALSE)
  }
  assert_scalar_number(cfg$seed, "seed", lower = 0,
                       upper = .Machine$integer.max, integer = TRUE)

  assert_nonempty_string(cfg$input$abundance_table, "input.abundance_table")
  assert_nonempty_string(cfg$input$params_json, "input.params_json")
  assert_nonempty_string(cfg$input$tax_column, "input.tax_column")
  for (key in c("phylogenetic_tree", "phylogenetic_tip_map")) {
    if (!is.null(cfg$input[[key]])) {
      assert_nonempty_string(cfg$input[[key]], paste0("input.", key))
    }
  }
  if (is.null(cfg$input$phylogenetic_tree) && !is.null(cfg$input$phylogenetic_tip_map)) {
    stop("'input.phylogenetic_tip_map' requires 'input.phylogenetic_tree'.", call. = FALSE)
  }
  assert_nonempty_string(cfg$output$base_dir, "output.base_dir")
  if (!is.list(cfg$krona)) {
    stop("'krona' must be a configuration mapping.", call. = FALSE)
  }
  if (!is.logical(cfg$krona$enabled) || length(cfg$krona$enabled) != 1L ||
      is.na(cfg$krona$enabled)) {
    stop("'krona.enabled' must be true or false.", call. = FALSE)
  }
  if (!is.logical(cfg$krona$render_html) || length(cfg$krona$render_html) != 1L ||
      is.na(cfg$krona$render_html)) {
    stop("'krona.render_html' must be true or false.", call. = FALSE)
  }
  assert_nonempty_string(cfg$krona$executable, "krona.executable")
  assert_nonempty_string(cfg$krona$html_renderer, "krona.html_renderer")
  if (!cfg$krona$html_renderer %in% c("builtin", "kronatools", "auto")) {
    stop("'krona.html_renderer' must be 'builtin', 'kronatools', or 'auto'.", call. = FALSE)
  }
  if (!is.list(cfg$pavian)) {
    stop("'pavian' must be a configuration mapping.", call. = FALSE)
  }
  if (!is.logical(cfg$pavian$enabled) || length(cfg$pavian$enabled) != 1L ||
      is.na(cfg$pavian$enabled)) {
    stop("'pavian.enabled' must be true or false.", call. = FALSE)
  }
  if (!is.logical(cfg$pavian$render_html) || length(cfg$pavian$render_html) != 1L ||
      is.na(cfg$pavian$render_html)) {
    stop("'pavian.render_html' must be true or false.", call. = FALSE)
  }
  if (!is.null(cfg$input$wf16s_output_root)) {
    assert_nonempty_string(cfg$input$wf16s_output_root, "input.wf16s_output_root")
  }
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
  if (!is.null(cfg$qc$target_min_length) && !is.null(cfg$qc$target_max_length) &&
      cfg$qc$target_min_length >= cfg$qc$target_max_length) {
    stop("'qc.target_min_length' must be smaller than 'qc.target_max_length'.", call. = FALSE)
  }
  if (cfg$qc$display_min_length >= cfg$qc$display_max_length) {
    stop("'qc.display_min_length' must be smaller than 'qc.display_max_length'.", call. = FALSE)
  }
  if (!is.null(cfg$resource_budgets)) {
    if (!is.list(cfg$resource_budgets)) stop("'resource_budgets' must be a mapping.", call. = FALSE)
    for (key in c("max_distance_cells", "max_permutation_cells", "max_rarefaction_rows")) {
      if (!is.null(cfg$resource_budgets[[key]])) {
        assert_scalar_number(cfg$resource_budgets[[key]], paste0("resource_budgets.", key),
                             lower = 1, integer = TRUE)
      }
    }
  }
  assert_scalar_number(cfg$alpha$rarefaction_points, "alpha.rarefaction_points", lower = 2, upper = 1000, integer = TRUE)
  assert_scalar_number(cfg$alpha$resample_depth, "alpha.resample_depth", lower = 1,
                       upper = .Machine$integer.max, integer = TRUE)
  assert_scalar_number(cfg$alpha$resample_fraction_cap, "alpha.resample_fraction_cap",
                       lower = 0, upper = 1, lower_open = TRUE)
  assert_scalar_number(cfg$alpha$resample_iterations, "alpha.resample_iterations", lower = 1, upper = 10000, integer = TRUE)
  validate_alpha_order_config <- function(value, name) {
    if (!is.numeric(value) || !length(value) || anyNA(value) ||
        any(!is.finite(value)) || any(value < 0) || anyDuplicated(value)) {
      stop(sprintf("'%s' must contain unique, finite, non-negative q orders.", name),
           call. = FALSE)
    }
  }
  validate_alpha_order_config(cfg$alpha$hill_orders, "alpha.hill_orders")
  validate_alpha_order_config(cfg$alpha$renyi_orders, "alpha.renyi_orders")
  if (!any(cfg$alpha$renyi_orders == 1)) {
    stop("'alpha.renyi_orders' must include q=1 for the fixed 02d figure contract.",
         call. = FALSE)
  }
  assert_scalar_number(cfg$composition$top_n_taxa, "composition.top_n_taxa", lower = 1, upper = 10000, integer = TRUE)
  validate_rank_vector <- function(value, name) {
    if (!is.character(value) || length(value) == 0L || anyNA(value) ||
        any(!nzchar(value)) || anyDuplicated(value) ||
        any(!value %in% SUPPORTED_COMPOSITION_RANKS)) {
      stop(sprintf("'%s' must contain unique supported rank names.", name), call. = FALSE)
    }
    invisible(TRUE)
  }
  legacy_heatmap_rank <- cfg$composition[["heatmap_rank"]]
  if (!is.null(legacy_heatmap_rank)) {
    validate_rank_vector(legacy_heatmap_rank, "composition.heatmap_rank")
    if (length(legacy_heatmap_rank) != 1L) {
      stop("'composition.heatmap_rank' must be null or one supported rank.", call. = FALSE)
    }
  }
  validate_rank_vector(cfg$composition$stacked_bar_ranks, "composition.stacked_bar_ranks")
  validate_rank_vector(cfg$composition$heatmap_ranks, "composition.heatmap_ranks")
  assert_scalar_number(cfg$composition$stacked_bar_max_taxa,
                       "composition.stacked_bar_max_taxa", lower = 1, upper = 10000,
                       integer = TRUE)
  assert_scalar_number(cfg$composition$stacked_bar_min_taxa,
                       "composition.stacked_bar_min_taxa", lower = 1, upper = 10000,
                       integer = TRUE)
  if (cfg$composition$stacked_bar_min_taxa > cfg$composition$stacked_bar_max_taxa) {
    stop("'composition.stacked_bar_min_taxa' must be less than or equal to 'composition.stacked_bar_max_taxa'.",
         call. = FALSE)
  }
  assert_scalar_number(cfg$composition$stacked_bar_min_mean_relative,
                       "composition.stacked_bar_min_mean_relative", lower = 0, upper = 1)
  assert_scalar_number(cfg$composition$heatmap_top_n_taxa,
                       "composition.heatmap_top_n_taxa", lower = 1, upper = 10000,
                       integer = TRUE)
  if (!is.logical(cfg$composition$heatmap_include_other) ||
      length(cfg$composition$heatmap_include_other) != 1L ||
      is.na(cfg$composition$heatmap_include_other)) {
    stop("'composition.heatmap_include_other' must be true or false.", call. = FALSE)
  }
  assert_nonempty_string(cfg$composition$heatmap_transform, "composition.heatmap_transform")
  if (!cfg$composition$heatmap_transform %in% c("log10_relative", "none")) {
    stop("'composition.heatmap_transform' must be 'log10_relative' or 'none'.", call. = FALSE)
  }
  assert_scalar_number(cfg$faprotax$top_n_functions, "faprotax.top_n_functions",
                       lower = 1, upper = 10000, integer = TRUE)
  if (!is.character(cfg$beta$distances) || length(cfg$beta$distances) == 0L ||
      any(!cfg$beta$distances %in% c("bray", "jaccard")) || anyDuplicated(cfg$beta$distances)) {
    stop("'beta.distances' must contain unique values drawn from: bray, jaccard.", call. = FALSE)
  }
  assert_nonempty_string(cfg$beta$primary_distance, "beta.primary_distance")
  if (sum(cfg$beta$distances == cfg$beta$primary_distance) != 1L) {
    stop("'beta.primary_distance' must occur exactly once in 'beta.distances'.", call. = FALSE)
  }
  assert_scalar_number(cfg$beta$permutations, "beta.permutations", lower = 1, upper = 1000000, integer = TRUE)
  assert_scalar_number(cfg$beta$minimum_count, "beta.minimum_count", lower = 1, integer = TRUE)
  if (!is.logical(cfg$beta$resampling$enabled) || length(cfg$beta$resampling$enabled) != 1L ||
      is.na(cfg$beta$resampling$enabled)) {
    stop("'beta.resampling.enabled' must be true or false.", call. = FALSE)
  }
  assert_scalar_number(cfg$beta$resampling$iterations, "beta.resampling.iterations", lower = 1, upper = 10000, integer = TRUE)
  assert_scalar_number(cfg$beta$resampling$depth_fraction_of_minimum,
                       "beta.resampling.depth_fraction_of_minimum",
                       lower = 0, upper = 1, lower_open = TRUE)
  assert_scalar_number(cfg$beta$resampling$minimum_success_fraction,
                       "beta.resampling.minimum_success_fraction",
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
      abundance_table = "wf16s-inputs/output_AAy/abundance_table_species.tsv",
      metadata = NULL,
      params_json = "wf16s-inputs/output_AAy/params.json",
      wf16s_output_root = NULL,
      tax_column = "tax",
      aggregate_columns = c("total"),
      include_samples = NULL,
      assignments = NULL,
      phylogenetic_tree = NULL,
      phylogenetic_tip_map = NULL
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
      resample_iterations = 100L,
      hill_orders = c(0, 1, 2),
      renyi_orders = c(1)
    ),
    composition = list(
      top_n_taxa = 15L,
      stacked_bar_ranks = c("phylum", "family", "genus"),
      stacked_bar_max_taxa = 15L,
      stacked_bar_min_taxa = 5L,
      stacked_bar_min_mean_relative = 0.005,
      heatmap_rank = NULL,
      heatmap_ranks = c("phylum", "family", "genus"),
      heatmap_top_n_taxa = 10L,
      heatmap_include_other = TRUE,
      heatmap_transform = "log10_relative"
    ),
    faprotax = list(
      top_n_functions = 20L
    ),
    krona = list(
      enabled = FALSE,
      render_html = TRUE,
      html_renderer = "builtin",
      executable = "ktImportText"
    ),
    pavian = list(
      enabled = FALSE,
      render_html = TRUE
    ),
    beta = list(
      distances = c("bray", "jaccard"),
      primary_distance = "bray",
      permutations = 999L,
      minimum_count = 1L,
      strata_column = NULL,
      resampling = list(
        enabled = FALSE,
        iterations = 100L,
        depth_fraction_of_minimum = 0.75,
        minimum_success_fraction = 0.50
      )
    ),
    shared_taxa = list(
      rank = "species",
      minimum_count = 1L,
      group_prevalence = 0.5
    ),
    resource_budgets = list(
      max_distance_cells = 10000000L,
      max_permutation_cells = 10000000L,
      max_rarefaction_rows = 10000000L
    ),
    taxonomy = list(
      cache = "wf16s-inputs/output_AAy/taxonomy_cache.json",
      network_mode = "cache_only",
      unresolved_policy = "warn",
      email_env = "NCBI_EMAIL",
      api_key_env = "NCBI_API_KEY"
    )
  )
}

migrate_composition_config <- function(raw_yaml) {
  if (is.null(raw_yaml$composition) || !is.list(raw_yaml$composition)) {
    return(raw_yaml)
  }

  raw_composition <- raw_yaml$composition
  has_legacy <- "heatmap_rank" %in% names(raw_composition)
  has_new <- "heatmap_ranks" %in% names(raw_composition)
  legacy_value <- raw_composition[["heatmap_rank"]]
  legacy_non_null <- has_legacy && !is.null(legacy_value)

  if (legacy_non_null && has_new) {
    stop(paste(
      "Configuration supplies both 'composition.heatmap_rank' and",
      "'composition.heatmap_ranks'; remove the legacy scalar key or the new vector."
    ), call. = FALSE)
  }
  if (legacy_non_null) {
    warning(paste(
      "'composition.heatmap_rank' is deprecated; use 'composition.heatmap_ranks'.",
      "The scalar value was migrated to a one-element rank vector."
    ), call. = FALSE)
    raw_composition$heatmap_ranks <- legacy_value
    raw_composition$heatmap_rank <- NULL
  }
  raw_yaml$composition <- raw_composition
  raw_yaml
}

merge_config <- function(default_cfg, user_cfg, path = "") {
  unknown <- setdiff(names(user_cfg), names(default_cfg))
  if (length(unknown) > 0L) {
    qualified <- paste0(path, unknown)
    stop(sprintf("Unknown configuration key(s): %s", paste(qualified, collapse = ", ")), call. = FALSE)
  }
  merged <- default_cfg
  for (key in names(user_cfg)) {
    if (is.null(user_cfg[[key]])) {
      # Preserve explicit null keys. Besides retaining the distinction between
      # an omitted and a deprecated setting, this prevents R's partial `$`
      # matching from confusing heatmap_rank with heatmap_ranks.
      merged[key] <- list(NULL)
    } else if (is.list(user_cfg[[key]]) && is.list(merged[[key]])) {
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
  if (!is.list(raw_yaml) || is.null(names(raw_yaml))) {
    stop("Configuration YAML root must be a named mapping.", call. = FALSE)
  }
  raw_yaml <- migrate_composition_config(raw_yaml)
  cfg <- merge_config(default_cfg, raw_yaml)
  cli_refresh <- isTRUE(cli_opts$refresh_taxonomy) ||
    isTRUE(cli_opts[["refresh-taxonomy"]])
  cli_krona <- isTRUE(cli_opts$krona) || isTRUE(cli_opts[["krona"]])
  cli_pavian <- isTRUE(cli_opts$pavian) || isTRUE(cli_opts[["pavian"]])

  if (identical(cfg$taxonomy$network_mode, "refresh") && !cli_refresh) {
    stop(paste(
      "YAML cannot enable taxonomy refresh.",
      "Keep taxonomy.network_mode: cache_only and pass --refresh-taxonomy explicitly."
    ), call. = FALSE)
  }

  # CLI overrides
  cli_output_dir_provided <- FALSE
  if (!is.null(cli_opts$output_dir) && nzchar(cli_opts$output_dir)) {
    cfg$output$base_dir <- cli_opts$output_dir
    cli_output_dir_provided <- TRUE
  } else if (!is.null(cli_opts[["output-dir"]]) && nzchar(cli_opts[["output-dir"]])) {
    cfg$output$base_dir <- cli_opts[["output-dir"]]
    cli_output_dir_provided <- TRUE
  }

  if (cli_refresh) {
    cfg$taxonomy$network_mode <- "refresh"
  }
  if (cli_krona) {
    cfg$krona$enabled <- TRUE
  }
  if (cli_pavian) {
    cfg$pavian$enabled <- TRUE
  }

  validate_config(cfg)

  requested_modules <- parse_requested_modules(cli_opts$modules)
  cfg$cli <- list(
    validate_only = isTRUE(cli_opts$validate_only) || isTRUE(cli_opts[["validate-only"]]),
    keep_going = isTRUE(cli_opts$keep_going) || isTRUE(cli_opts[["keep-going"]]),
    overwrite = isTRUE(cli_opts$overwrite),
    allow_unlocked = isTRUE(cli_opts$allow_unlocked) || isTRUE(cli_opts[["allow-unlocked"]]),
    allow_dirty = isTRUE(cli_opts$allow_dirty) || isTRUE(cli_opts[["allow-dirty"]]),
    allow_large_workload = isTRUE(cli_opts$allow_large_workload) || isTRUE(cli_opts[["allow-large-workload"]]),
    online_preflight = isTRUE(cli_opts$online_preflight) || isTRUE(cli_opts[["online-preflight"]]),
    refresh_taxonomy = cli_refresh,
    krona = isTRUE(cfg$krona$enabled),
    pavian = isTRUE(cfg$pavian$enabled),
    modules = requested_modules
  )

  # Path resolution against config_dir
  cfg$input$abundance_table <- resolve_path(cfg$input$abundance_table, config_dir)
  cfg$input$metadata <- resolve_path(cfg$input$metadata, config_dir)
  cfg$input$params_json <- resolve_path(cfg$input$params_json, config_dir)
  cfg$input$wf16s_output_root <- resolve_path(cfg$input$wf16s_output_root, config_dir)
  cfg$input$phylogenetic_tree <- resolve_path(cfg$input$phylogenetic_tree, config_dir)
  cfg$input$phylogenetic_tip_map <- resolve_path(cfg$input$phylogenetic_tip_map, config_dir)

  if (!is.null(cfg$input$assignments) && is.list(cfg$input$assignments)) {
    for (s in names(cfg$input$assignments)) {
      cfg$input$assignments[[s]] <- resolve_path(cfg$input$assignments[[s]], config_dir)
    }
  }

  cfg$taxonomy$cache <- resolve_path(cfg$taxonomy$cache, config_dir)

  # Resolve base output directory
  base_dir_context <- if (cli_output_dir_provided) getwd() else config_dir
  cfg$output$base_dir <- resolve_path(cfg$output$base_dir, base_dir_context)

  # Derive all module output directories
  base_out <- cfg$output$base_dir
  cfg$output$dirs <- list(
    qc = file.path(base_out, "01_QC"),
    alpha = file.path(base_out, "02_Alpha_Diversity"),
    beta = file.path(base_out, "03_Beta_Diversity"),
    composition = file.path(base_out, "04_Taxa_Composition"),
    ordination = file.path(base_out, "05_Ordination"),
    shared_taxa = file.path(base_out, "06_Shared_Taxa"),
    kreport = file.path(base_out, "07_Kreport"),
    faprotax = file.path(base_out, "08_FAPROTAX")
  )

  cfg$output$manifest_file <- file.path(base_out, "run_manifest.json")
  cfg$output$resolved_config_file <- file.path(base_out, "resolved_config.yml")
  cfg$output$session_info_file <- file.path(base_out, "session_info.txt")

  cfg$config_file <- config_file_abs
  cfg$config_dir <- config_dir

  cfg
}
