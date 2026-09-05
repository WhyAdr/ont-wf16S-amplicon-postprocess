# =============================================================================
# Regression Tests: Section 2.3 Numerical Targets
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "metrics.R"))
source(file.path("..", "..", "analysis", "utils", "plotting.R"))
source(file.path("..", "..", "analysis", "02_alpha_diversity.R"))

test_that("Real Ambar Ayunda fixture reproduces exact Section 2.3 alpha regression metrics", {
  ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
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
  ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
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

test_that("Seeded rarefaction resample output is byte-stable", {
  root <- tempfile("alpha_determinism_")
  dir.create(root)
  abundance <- create_temp_abundance(root, n_species = 8, sample_names = "S1")

  run_once <- function(output_name) {
    cfg <- get_default_config()
    cfg$input$abundance_table <- abundance
    cfg$alpha$resample_depth <- 100L
    cfg$alpha$resample_iterations <- 10L
    cfg$output$base_dir <- file.path(root, output_name)
    cfg$output$dirs <- list(alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
    run_alpha(build_context(cfg))
    file.path(cfg$output$dirs$alpha, "rarefaction_resamples.tsv")
  }

  first <- run_once("first")
  second <- run_once("second")
  expect_equal(
    digest::digest(file = first, algo = "sha256"),
    digest::digest(file = second, algo = "sha256")
  )
})
