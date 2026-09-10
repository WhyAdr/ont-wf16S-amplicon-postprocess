# =============================================================================
# Regression Tests: Section 2.3 Numerical Targets
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "metrics.R"))
source(file.path("..", "..", "analysis", "utils", "alpha_phylogeny.R"))
source(file.path("..", "..", "analysis", "utils", "plotting.R"))
source(file.path("..", "..", "analysis", "utils", "preflight.R"))
source(file.path("..", "..", "analysis", "02_alpha_diversity.R"))

test_that("Real Ambar Ayunda fixture reproduces exact Section 2.3 alpha regression metrics", {
  ab_path <- file.path("..", "..", "wf16s-inputs", "output_AAy", "abundance_table_species.tsv")
  skip_if_not(file.exists(ab_path), "Real abundance table not found")

  res <- read_abundance_table(ab_path)
  sample_col <- res$samples[1]
  counts <- res$count_matrix[, sample_col]
  class_counts <- counts[-res$unclass_index]

  alpha_df <- calc_alpha_indices(class_counts)
  vals <- setNames(alpha_df$Value, alpha_df$Metric)

  # Target values from Section 2.3:
  # Observed species richness: 1,836
  expect_equal(vals[["Observed species richness (S)"]], 1836)

  # Chao1: 2,983.851
  expect_equal(round(vals[["Chao1 (estimated richness)"]], 3), 2983.851, tolerance = 0.001)

  # Shannon: 4.603
  expect_equal(round(vals[["Shannon (H)"]], 3), 4.603, tolerance = 0.001)

  # Effective species number: 99.814
  expect_equal(round(vals[["Effective number of species (e^H)"]], 3), 99.814, tolerance = 0.001)

  # Simpson: 0.957
  expect_equal(round(vals[["Simpson's D (1-sum p^2)"]], 3), 0.957, tolerance = 0.001)

  # Inverse Simpson: 23.089
  expect_equal(round(vals[["Inverse Simpson"]], 3), 23.089, tolerance = 0.001)

  # Pielou evenness: 0.613
  expect_equal(round(vals[["Pielou's evenness (J)"]], 3), 0.613, tolerance = 0.001)

  # Fisher alpha: 334.543
  expect_equal(round(vals[["Fisher's alpha"]], 3), 334.543, tolerance = 0.001)

  # Berger-Parker dominance: 0.133
  expect_equal(round(vals[["Berger-Parker dominance"]], 3), 0.133, tolerance = 0.001)
})

test_that("Detected classified taxa count by rank matches Section 2.3 targets", {
  ab_path <- file.path("..", "..", "wf16s-inputs", "output_AAy", "abundance_table_species.tsv")
  skip_if_not(file.exists(ab_path), "Real abundance table not found")

  res <- read_abundance_table(ab_path)
  tax_df <- res$taxonomy[-res$unclass_index, ]

  # Expected: 31 phyla, 69 classes, 128 orders, 281 families, 867 genera, 1836 species
  expect_equal(length(unique(tax_df$phylum)), 31)
  expect_equal(length(unique(tax_df$class)), 69)
  expect_equal(length(unique(tax_df$order)), 128)
  expect_equal(length(unique(tax_df$family)), 281)
  expect_equal(length(unique(tax_df$genus)), 867)
  expect_equal(length(unique(tax_df$species)), 1836)
})

test_that("Ambar richness overview exposes the exact low-count tail", {
  ab_path <- file.path("..", "..", "wf16s-inputs", "output_AAy", "abundance_table_species.tsv")
  skip_if_not(file.exists(ab_path), "Real abundance table not found")
  parsed <- read_abundance_table(ab_path)
  class_matrix <- parsed$count_matrix[-parsed$unclass_index, , drop = FALSE]
  observed <- build_richness_overview(class_matrix, parsed$samples)
  expect_equal(observed$ClassifiedReads, 80556)
  expect_equal(observed$PositiveTaxa, 1836)
  expect_equal(observed$SingletonTaxa, 735)
  expect_equal(observed$TaxaLeq10, 1399)
  expect_equal(observed$ReadsInTaxaLeq10, 3456)
  expect_equal(observed$ReadsInTaxaLeq10Pct, 4.2901832265753, tolerance = 1e-10)
})

test_that("Seeded rarefaction resample output is byte-stable", {
  root <- tempfile("alpha_determinism_")
  dir.create(root)
  abundance <- create_temp_abundance(root, n_species = 8, sample_names = "S1")

  run_once <- function(output_name) {
    cfg <- get_default_config()
    cfg$input$abundance_table <- abundance
    cfg$input$params_json <- create_temp_params(root)
    cfg$alpha$resample_depth <- 100L
    cfg$alpha$resample_iterations <- 10L
    cfg$output$base_dir <- file.path(root, output_name)
    cfg$output$dirs <- list(alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
    run_alpha(build_context(cfg))
    file.path(cfg$output$dirs$alpha, c(
      "rarefaction_resamples.tsv", "rarefaction_resamples_long.tsv",
      "02b_resample_boxplots.png", "02c_resample_boxplots_estimators_evenness.png",
      "02d_resample_boxplots_evenness_entropy.png"
    ))
  }

  first <- run_once("first")
  second <- run_once("second")
  first_hashes <- vapply(first, digest::digest, character(1), file = TRUE, algo = "sha256")
  second_hashes <- vapply(second, digest::digest, character(1), file = TRUE, algo = "sha256")
  expect_identical(unname(first_hashes), unname(second_hashes))
})

test_that("Derived sample seeds are deterministic throughout the set.seed range", {
  samples <- c("S1", "sample with spaces", "another-sample")
  derived <- vapply(samples, function(sample_id) {
    derive_sample_seed(.Machine$integer.max, sample_id)
  }, integer(1))
  expect_true(all(is.finite(derived)))
  expect_true(all(derived >= 0L & derived <= .Machine$integer.max))
  expect_identical(
    derive_sample_seed(.Machine$integer.max, "S1"),
    derive_sample_seed(.Machine$integer.max, "S1")
  )
  expect_error(derive_sample_seed(.Machine$integer.max + 1, "S1"), "set.seed\\(\\) range")
})

test_that("registry metrics satisfy closed-form equal and skewed-count contracts", {
  equal <- calc_alpha_metric_long(c(a = 1, b = 1, c = 1))
  value <- setNames(equal$Value, equal$MetricKey)
  expect_equal(value[["richness"]], 3)
  expect_equal(value[["shannon"]], log(3))
  expect_equal(value[["ens"]], 3)
  expect_equal(value[["simpson"]], 2 / 3)
  expect_equal(value[["invsimpson"]], 3)
  expect_equal(value[["pielou"]], 1)
  expect_equal(value[["heip"]], 1)
  expect_equal(value[["evar"]], 1)
  expect_equal(value[["mcintosh"]], 1)
  expect_equal(value[["hill_q0"]], value[["richness"]])
  expect_equal(value[["hill_q1"]], value[["ens"]])
  expect_equal(value[["hill_q2"]], value[["invsimpson"]])
  expect_equal(value[["renyi_q1"]], value[["shannon"]])

  skewed <- calc_alpha_metric_long(c(a = 8, b = 1, c = 1))
  skew <- setNames(skewed$Value, skewed$MetricKey)
  p <- c(0.8, 0.1, 0.1)
  expected_h <- -sum(p * log(p))
  expect_equal(skew[["shannon"]], expected_h)
  expect_equal(skew[["berger_parker"]], 0.8)
  expect_equal(skew[["evar"]],
               1 - (2 / pi) * atan(mean((log(c(8, 1, 1)) - mean(log(c(8, 1, 1))))^2)))
})

test_that("undefined alpha domains are explicit and never NaN or infinite", {
  empty <- calc_alpha_metric_long(numeric(0))
  expect_false(any(empty$Valid))
  expect_true(all(is.na(empty$Value)))
  expect_true(all(nzchar(empty$Reason)))

  singleton <- calc_alpha_metric_long(c(only = 1))
  invalid <- singleton[!singleton$Valid, ]
  expect_true(all(invalid$MetricKey %in% c(
    "fisher_alpha", "chao1", "ace", "pielou", "heip", "evar", "mcintosh"
  )))
  expect_false(any(is.nan(singleton$Value)))
  expect_false(any(is.infinite(singleton$Value), na.rm = TRUE))
  expect_error(calc_alpha_metric_long(c(1, 1.5)), "non-negative integers")
  expect_error(calc_alpha_metric_long(c(1, -1)), "non-negative integers")

  large <- calc_alpha_metric_long(c(large = .Machine$integer.max + 1, small = 1))
  expect_equal(large$Value[large$MetricKey == "richness"], 2)
  expect_true(all(is.finite(large$Value[large$Valid])))
})

test_that("a depth-one sample retains a valid rarefaction sensitivity contract", {
  resampled <- calc_rarefaction_resamples(c(only = 1), subsample_depth = 1,
                                          n_iterations = 3, seed = 42)
  expect_equal(resampled$subsample_depth, rep(1L, 3))
  long <- attr(resampled, "alpha_long")
  expect_equal(unique(long$subsample_depth), 1L)
  expect_true(all(long$Value[long$MetricKey == "richness"] == 1))
})

test_that("phylogenetic alpha uses rooted branch lengths and explicit path mapping", {
  root <- tempfile("alpha_phylo_")
  dir.create(root)
  tree_path <- file.path(root, "tree.nwk")
  map_path <- file.path(root, "map.tsv")
  writeLines("((t1:1,t2:1):1,t3:2);", tree_path)
  map <- data.frame(TaxonPath = c("P1", "P2", "P3"),
                    TipLabel = c("t1", "t2", "t3"))
  write.table(map, map_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cfg <- get_default_config()
  cfg$input$phylogenetic_tree <- tree_path
  cfg$input$phylogenetic_tip_map <- map_path
  phylogeny <- prepare_alpha_phylogeny(cfg, map$TaxonPath)
  observed <- calc_phylogenetic_alpha(c(P1 = 1, P2 = 1, P3 = 1), phylogeny)
  value <- setNames(observed$Value, observed$MetricKey)
  expect_equal(value[["faith_pd"]], 5)
  expect_equal(value[["psr"]], 2.5)
  expect_equal(value[["pse"]], 5 / 6)
  expect_true(all(observed$Valid))

  writeLines("((t1:1,t2:1):1,unused:3);", tree_path)
  two_map <- data.frame(TaxonPath = c("P1", "P2"), TipLabel = c("t1", "t2"))
  write.table(two_map, map_path, sep = "\t", row.names = FALSE, quote = FALSE)
  rooted <- prepare_alpha_phylogeny(cfg, two_map$TaxonPath)
  rooted_value <- setNames(
    calc_phylogenetic_alpha(c(P1 = 1, P2 = 1), rooted)$Value,
    c("faith_pd", "psr", "pse")
  )
  expect_equal(rooted_value[["faith_pd"]], 3)
  expect_equal(rooted_value[["psr"]], 1)
  expect_equal(rooted_value[["pse"]], 0.5)
  expect_equal(rooted$pruned_tip_count, 1)
  expect_match(rooted$pruning_method, "original root retained")

  cfg$input$phylogenetic_tip_map <- NULL
  expect_error(prepare_alpha_phylogeny(cfg, map$TaxonPath), "Mapped phylogenetic tip")
  cfg$input$phylogenetic_tree <- NULL
  skipped <- prepare_alpha_phylogeny(cfg, map$TaxonPath)
  expect_false(skipped$enabled)
  expect_equal(skipped$skip$Status, "Skipped")
})

test_that("alpha module emits canonical long tables and fixed non-phylogenetic figures", {
  root <- tempfile("alpha_contract_")
  dir.create(root)
  cfg <- get_default_config()
  cfg$input$abundance_table <- create_temp_abundance(root, n_species = 8, sample_names = "S1")
  cfg$input$params_json <- create_temp_params(root)
  cfg$alpha$resample_depth <- 100L
  cfg$alpha$resample_iterations <- 4L
  cfg$output$base_dir <- file.path(root, "out")
  cfg$output$dirs <- list(alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
  result <- run_alpha(build_context(cfg))
  output_names <- basename(result$outputs)
  expect_true(all(c(
    "alpha_metric_definitions.tsv", "alpha_metric_status.tsv",
    "rarefaction_resamples.tsv", "rarefaction_resamples_long.tsv",
    "02b_resample_boxplots.png", "02c_resample_boxplots_estimators_evenness.png",
    "02d_resample_boxplots_evenness_entropy.png", "alpha_phylogeny_skipped.tsv",
    "alpha_plot_registry.tsv"
  ) %in% output_names))
  expect_false("02e_resample_boxplots_phylogenetic.png" %in% output_names)
  long <- read.delim(file.path(cfg$output$dirs$alpha, "rarefaction_resamples_long.tsv"),
                     check.names = FALSE)
  expect_identical(names(long), c(
    "SampleID", "iteration", "subsample_depth", "MetricKey", "MetricLabel",
    "Parameter", "Value", "Valid", "Reason"
  ))
  expect_equal(nrow(long), cfg$alpha$resample_iterations * 17L)
  expect_equal(anyDuplicated(long[c("SampleID", "iteration", "MetricKey")]), 0L)
  expect_true(all(is.finite(long$Value[long$Valid])))
  plot_contract <- read.delim(file.path(cfg$output$dirs$alpha, "alpha_plot_registry.tsv"),
                              check.names = FALSE)
  expect_identical(plot_contract$MetricKey[plot_contract$FigureKey == "02b"],
                   c("richness", "shannon", "ens", "simpson", "invsimpson", "fisher_alpha"))
})

test_that("a valid tree and exact tip map produce the complete 02e contract", {
  skip_if_not_installed("ape")
  root <- tempfile("alpha_phylo_module_")
  dir.create(root)
  abundance <- create_temp_abundance(root, n_species = 5, sample_names = "S1")
  parsed <- read_abundance_table(abundance)
  paths <- parsed$taxonomy$TaxonPath[-parsed$unclass_index]
  set.seed(406)
  tree <- ape::rtree(length(paths), tip.label = paste0("tip", seq_along(paths)))
  tree_path <- file.path(root, "tree.nwk")
  map_path <- file.path(root, "tip-map.tsv")
  ape::write.tree(tree, tree_path)
  write.table(data.frame(TaxonPath = paths, TipLabel = tree$tip.label), map_path,
              sep = "\t", row.names = FALSE, quote = FALSE)

  cfg <- get_default_config()
  cfg$input$abundance_table <- abundance
  cfg$input$params_json <- create_temp_params(root)
  cfg$input$phylogenetic_tree <- tree_path
  cfg$input$phylogenetic_tip_map <- map_path
  cfg$alpha$resample_depth <- 50L
  cfg$alpha$resample_iterations <- 3L
  cfg$output$base_dir <- file.path(root, "out")
  cfg$output$dirs <- list(alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
  context <- build_context(cfg)
  budget_cfg <- cfg
  budget_cfg$resource_budgets$max_rarefaction_rows <- 59L
  budget_context <- context
  budget_context$config <- budget_cfg
  expect_error(run_module_preflight(budget_context, "alpha"), "60 rows")
  result <- run_alpha(context)
  expect_true("02e_resample_boxplots_phylogenetic.png" %in% basename(result$outputs))
  expect_false("alpha_phylogeny_skipped.tsv" %in% basename(result$outputs))
  long <- read.delim(file.path(cfg$output$dirs$alpha, "rarefaction_resamples_long.tsv"),
                     check.names = FALSE)
  expect_true(all(c("faith_pd", "psr", "pse") %in% long$MetricKey))
  expect_equal(nrow(long), 3L * 20L)
  provenance <- read.delim(file.path(cfg$output$dirs$alpha, "alpha_phylogeny_status.tsv"),
                           check.names = FALSE)
  expect_equal(provenance$Status, "Completed")
  expect_match(provenance$TreeSHA256, "^[0-9a-f]{64}$")
})
