# =============================================================================
# Unit Tests: Cohort Modules & Synthetic Multi-Sample Gating
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "metrics.R"))
source(file.path("..", "..", "analysis", "utils", "plotting.R"))
source(file.path("..", "..", "analysis", "02_alpha_diversity.R"))
source(file.path("..", "..", "analysis", "03_beta_diversity.R"))
source(file.path("..", "..", "analysis", "04_taxa_composition.R"))
source(file.path("..", "..", "analysis", "05_ordination.R"))
source(file.path("..", "..", "analysis", "06_shared_taxa.R"))

test_that("Single-sample context skips beta, ordination, and shared_taxa gracefully", {
  tmp <- tempdir()
  ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = c("S1"))

  cfg <- get_default_config()
  cfg$input$abundance_table <- ab_file
  cfg$input$params_json <- create_temp_params(tmp)
  cfg$output$base_dir <- file.path(tmp, "out_single")
  cfg$output$dirs <- list(
    beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"),
    alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"),
    ordination = file.path(cfg$output$base_dir, "05_Ordination"),
    shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa")
  )

  context <- build_context(cfg)
  expect_equal(context$mode, "single")

  # Run beta
  res_beta <- run_beta(context)
  expect_equal(res_beta$status, "skipped")
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "beta_diversity_skipped.tsv")))

  # Run ordination
  res_ord <- run_ordination(context)
  expect_equal(res_ord$status, "skipped")
  expect_true(file.exists(file.path(cfg$output$dirs$ordination, "ordination_skipped.tsv")))

  # Run shared taxa
  res_shared <- run_shared_taxa(context)
  expect_equal(res_shared$status, "skipped")
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "shared_taxa_skipped.tsv")))
})

test_that("Synthetic cohort (2 groups x 3 replicates) passes cohort gates and pairs PERMANOVA with betadisper", {
  tmp <- tempdir()
  sample_names <- c("Ctrl1", "Ctrl2", "Ctrl3", "Trt1", "Trt2", "Trt3")
  ab_file <- create_temp_abundance(tmp, n_species = 20, sample_names = sample_names)
  meta_file <- create_temp_metadata(tmp, sample_names = sample_names,
                                    groups = c("Control", "Control", "Control", "Treated", "Treated", "Treated"))

  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- ab_file
  cfg$input$metadata <- meta_file
  cfg$input$params_json <- create_temp_params(tmp)
  cfg$output$base_dir <- file.path(tmp, "out_cohort")
  cfg$output$dirs <- list(
    alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"),
    beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"),
    ordination = file.path(cfg$output$base_dir, "05_Ordination"),
    shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa")
  )

  context <- build_context(cfg)
  expect_equal(context$mode, "cohort")
  expect_equal(length(context$samples), 6)

  # Run cohort alpha, including the metadata-aware group test path.
  res_alpha <- run_alpha(context)
  expect_equal(res_alpha$status, "completed")
  alpha_out <- read.delim(file.path(cfg$output$dirs$alpha, "alpha_diversity.tsv"), check.names = FALSE)
  expect_equal(alpha_out$SampleID, sample_names)
  expect_true("Group" %in% names(alpha_out))
  expect_false(any(c("Group.x", "Group.y") %in% names(alpha_out)))

  # Run beta
  res_beta <- run_beta(context)
  expect_equal(res_beta$status, "completed")
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "permanova.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "betadisper.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "pcoa_scores_bray.tsv")))

  # Run ordination
  res_ord <- run_ordination(context)
  expect_equal(res_ord$status, "completed")
  expect_true(file.exists(file.path(cfg$output$dirs$ordination, "pca_scores.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$ordination, "pca_variance.tsv")))
  nmds_diag <- file.path(cfg$output$dirs$ordination, "nmds_diagnostics.tsv")
  expect_true(file.exists(nmds_diag))
  nmds <- read.delim(nmds_diag, check.names = FALSE)
  expect_equal(nmds$Converged[1], nmds$BestSolutionRepetitions[1] > 0)

  # Run shared taxa
  res_shared <- run_shared_taxa(context)
  expect_equal(res_shared$status, "completed")
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "core_taxa.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "unique_taxa.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "group_prevalence.tsv")))
})

test_that("Forced cohort mode rejects a one-sample table", {
  tmp <- tempfile("forced_cohort_")
  dir.create(tmp)
  ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = "Only1")
  meta_file <- create_temp_metadata(tmp, sample_names = "Only1", groups = "Control")
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- ab_file
  cfg$input$metadata <- meta_file
  cfg$input$params_json <- create_temp_params(tmp)
  expect_error(build_context(cfg), "at least 2")
})

test_that("Under-replicated cohort skips PERMANOVA with explicit reason", {
  tmp <- tempdir()
  # 2 groups with only 1 sample each
  sample_names <- c("Ctrl1", "Trt1")
  ab_file <- create_temp_abundance(tmp, n_species = 10, sample_names = sample_names)
  meta_file <- create_temp_metadata(tmp, sample_names = sample_names, groups = c("Control", "Treated"))

  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- ab_file
  cfg$input$metadata <- meta_file
  cfg$input$params_json <- create_temp_params(tmp)
  cfg$output$base_dir <- file.path(tmp, "out_underrep")
  cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))

  context <- build_context(cfg)
  res_beta <- run_beta(context)

  perm_tsv <- file.path(cfg$output$dirs$beta, "permanova.tsv")
  expect_true(file.exists(perm_tsv))
  perm_df <- read.delim(perm_tsv)
  expect_equal(perm_df$Status[1], "Skipped")
  expect_match(perm_df$Reason[1], "at least 2 groups with >= 2 samples")
})

test_that("Cohort distance results are invariant to sample and metadata order", {
  root <- tempfile("cohort_order_")
  dir.create(root)
  samples <- c("Ctrl1", "Ctrl2", "Ctrl3", "Trt1", "Trt2", "Trt3")
  abundance_a <- create_temp_abundance(root, n_species = 12, sample_names = samples)
  abundance_df <- read.delim(abundance_a, check.names = FALSE)
  permuted <- rev(samples)
  abundance_b <- file.path(root, "synthetic_abundance_permuted.tsv")
  write.table(
    abundance_df[, c("tax", permuted, "total")], abundance_b,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  groups <- c(Ctrl1 = "Control", Ctrl2 = "Control", Ctrl3 = "Control",
              Trt1 = "Treated", Trt2 = "Treated", Trt3 = "Treated")
  metadata_a <- create_temp_metadata(root, samples, unname(groups[samples]))
  metadata_b <- file.path(root, "metadata_permuted.tsv")
  write.table(
    data.frame(SampleID = permuted, Group = unname(groups[permuted])), metadata_b,
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  run_beta_fixture <- function(abundance, metadata, output_name) {
    cfg <- get_default_config()
    cfg$mode <- "cohort"
    cfg$input$abundance_table <- abundance
    cfg$input$metadata <- metadata
    cfg$input$params_json <- create_temp_params(root)
    cfg$beta$permutations <- 19L
    cfg$output$base_dir <- file.path(root, output_name)
    cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
    context <- build_context(cfg)
    expect_equal(run_beta(context)$status, "completed")
    result <- read.delim(
      file.path(cfg$output$dirs$beta, "distance_bray.tsv"),
      check.names = FALSE
    )
    result <- result[order(result$SampleID), c("SampleID", sort(samples))]
    rownames(result) <- NULL
    result
  }

  expect_equal(
    run_beta_fixture(abundance_a, metadata_a, "first"),
    run_beta_fixture(abundance_b, metadata_b, "second"),
    tolerance = 1e-12
  )
})

test_that("Shared taxa preserves numeric-like sample IDs and arbitrary group labels", {
  root <- tempfile("shared_identity_")
  dir.create(root)
  samples <- c("01", "02")
  abundance <- create_temp_abundance(root, n_species = 4, sample_names = samples)
  metadata <- create_temp_metadata(root, samples, c("0", "Group-one"))
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- abundance
  cfg$input$metadata <- metadata
  cfg$input$params_json <- create_temp_params(root)
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa"))
  context <- build_context(cfg)
  expect_identical(context$metadata$SampleID, samples)
  expect_identical(context$metadata$Group, c("0", "Group-one"))
  expect_equal(run_shared_taxa(context)$status, "completed")
  prevalence <- read.delim(file.path(cfg$output$dirs$shared_taxa, "group_prevalence.tsv"),
                           check.names = FALSE)
  expect_true(all(c("0", "Group-one") %in% names(prevalence)))
})

test_that("One-row composition heatmaps do not attempt row clustering", {
  root <- tempfile("one_row_heatmap_")
  dir.create(root)
  samples <- c("S1", "S2")
  abundance <- create_temp_abundance(root, n_species = 3, sample_names = samples)
  metadata <- create_temp_metadata(root, samples, c("Control", "Treated"))
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- abundance
  cfg$input$metadata <- metadata
  cfg$input$params_json <- create_temp_params(root)
  cfg$composition$top_n_taxa <- 1L
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(composition = file.path(cfg$output$base_dir, "04_Taxa_Composition"))
  context <- build_context(cfg)
  expect_equal(run_taxa_composition(context)$status, "completed")
  heatmap <- file.path(cfg$output$dirs$composition, "04_heatmap_genus.png")
  expect_true(file.exists(heatmap))
  expect_gt(file.info(heatmap)$size, 0)
})

test_that("NMDS repetition helper follows vegan metaMDS semantics", {
  expect_false(nmds_solution_repeated(0L))
  expect_true(nmds_solution_repeated(1L))
  expect_true(nmds_solution_repeated(2L))
  expect_false(nmds_solution_repeated(NULL))
  expect_false(nmds_solution_repeated(NA))
  expect_false(nmds_solution_repeated("not-a-number"))
})

test_that("Beta rarefaction stability records isolated Procrustes failures", {
  root <- tempfile("beta_resampling_failure_")
  dir.create(root)
  samples <- c("Ctrl1", "Ctrl2", "Ctrl3", "Trt1", "Trt2", "Trt3")
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- create_temp_abundance(root, n_species = 12, sample_names = samples)
  cfg$input$metadata <- create_temp_metadata(
    root, samples, c("Control", "Control", "Control", "Treated", "Treated", "Treated")
  )
  cfg$input$params_json <- create_temp_params(root)
  cfg$beta$permutations <- 19L
  cfg$beta$resampling$enabled <- TRUE
  cfg$beta$resampling$iterations <- 3L
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
  context <- build_context(cfg)

  calls <- 0L
  fail_once <- function(X, Y) {
    calls <<- calls + 1L
    if (calls == 1L) stop("injected Procrustes failure")
    vegan::procrustes(X, Y)
  }
  expect_equal(run_beta(context, procrustes_fn = fail_once)$status, "completed")
  diagnostics <- read.delim(
    file.path(cfg$output$dirs$beta, "pcoa_rarefaction_diagnostics.tsv"), check.names = FALSE
  )
  stability <- read.delim(
    file.path(cfg$output$dirs$beta, "pcoa_rarefaction_stability.tsv"), check.names = FALSE
  )
  expect_equal(diagnostics$SuccessfulIterations, 2L)
  expect_identical(as.character(diagnostics$FailedIterations), "1")
  expect_setequal(unique(stability$Iteration), c(2L, 3L))

  expect_error(
    run_beta(context, procrustes_fn = function(X, Y) stop("always fail")),
    "All beta-diversity rarefaction stability iterations failed"
  )
})
