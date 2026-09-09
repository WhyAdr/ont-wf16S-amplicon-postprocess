# =============================================================================
# Unit Tests: Deterministic Composition Contracts
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "plotting.R"))
source(file.path("..", "..", "analysis", "04_taxa_composition.R"))

composition_rank_fixture <- function() {
  data.frame(
    TaxonPath = c("Bacteria;P1;Shared", "Bacteria;P2;Shared", "Bacteria;P1;Rare"),
    Taxon = c("Shared", "Shared", "Rare"),
    S2 = c(0.50, 0.25, 0.25),
    S10 = c(0.25, 0.50, 0.25),
    S1 = c(0, 0, 0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("TaxonPath identity survives contextual display-label disambiguation", {
  mapping <- make_taxon_display_map(composition_rank_fixture())
  expect_equal(nrow(mapping), 3L)
  expect_equal(length(unique(mapping$DisplayTaxon)), 3L)
  expect_true(all(grepl("Shared", mapping$DisplayTaxon[1:2], fixed = TRUE)))
  expect_false(any(mapping$DisplayTaxon == "Other"))
})

test_that("Taxon selection uses valid-sample means and TaxonPath tie breaks", {
  table <- composition_rank_fixture()
  selected <- select_display_taxa(
    table, c("S2", "S10", "S1"), valid_samples = c("S2", "S10"),
    min_mean = 0.25, min_n = 2L, max_n = 2L
  )
  expect_equal(selected$TaxonPath, c("Bacteria;P1;Shared", "Bacteria;P2;Shared"))
  expect_equal(selected$SelectionOrder, c(1L, 2L))
})

test_that("Composition collapse conserves valid classified-read abundance", {
  table <- composition_rank_fixture()
  metadata <- data.frame(
    SampleID = c("S2", "S10", "S1"),
    Group = c("B", "A", "B"),
    ValidDenominator = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  collapsed <- collapse_rank_for_display(
    table, c("Bacteria;P1;Shared"), c("S2", "S10", "S1"), metadata
  )
  sums <- aggregate(RelativeAbundance ~ SampleID, collapsed, sum)
  expect_equal(sums$RelativeAbundance[sums$SampleID %in% c("S2", "S10")], c(1, 1), tolerance = 1e-12)
  expect_equal(sums$RelativeAbundance[sums$SampleID == "S1"], 0)
  expect_false(any(collapsed$ValidDenominator[collapsed$SampleID == "S1"]))
  expect_equal(max(collapsed$StackOrder), nrow(collapsed[collapsed$SampleID == "S2", , drop = FALSE]))
  expect_true(all(collapsed$TaxonPath[collapsed$IsOther] == "__OTHER__"))
})

test_that("Cohort order, arithmetic group means, and heatmap contracts are explicit", {
  table <- composition_rank_fixture()
  metadata <- data.frame(
    SampleID = c("S2", "S10", "S1"),
    Group = c("B", "A", "B"),
    ValidDenominator = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  collapsed <- collapse_rank_for_display(
    table, c("Bacteria;P1;Shared", "Bacteria;P2;Shared"),
    c("S2", "S10", "S1"), metadata
  )
  means <- summarize_group_means(collapsed)
  expect_equal(unique(means$Group), c("B", "A"))
  expect_equal(unique(means$SamplesUsed[means$Group == "B"]), 1L)
  expect_equal(unique(means$SamplesExcludedZeroClassified[means$Group == "B"]), 1L)

  order <- ordered_cohort_samples(c("S2", "S10", "S1"), metadata)
  expect_equal(order$sample_order, c("S2", "S1", "S10"))
  expect_equal(order$group_order, c("B", "A"))

  prepared <- prepare_heatmap_data(
    table, c("S2", "S10", "S1"), metadata, top_n = 2L,
    include_other = TRUE, transform = "log10_relative", rank = "genus"
  )
  expect_false(prepared$skipped)
  expect_equal(prepared$sample_order, c("S2", "S1", "S10"))
  expect_true(isTRUE(prepared$cohort))
  expect_false(any(prepared$sidecar$ColumnOrder > 3L))
  expect_true(all(is.finite(prepared$sidecar$TransformedValue)))
  expect_true(all(prepared$sidecar$PseudoCount > 0))

  none <- prepare_heatmap_data(
    table, c("S2", "S10", "S1"), metadata, top_n = 1L,
    include_other = FALSE, transform = "none", rank = "genus"
  )
  expect_true(all(none$sidecar$PseudoCount == 0))
})

test_that("Cohort composition emits all configured rank figures and sidecars", {
  root <- tempfile("composition_integration_")
  dir.create(root)
  samples <- c("S2", "S10", "S1")
  paths <- c(
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown",
    "Bacteria;Bacillota;P1;Class1;Order1;Family1;Genus1;Species1",
    "Bacteria;Pseudomonadota;P2;Class2;Order2;Family2;Genus2;Species2"
  )
  taxonomy <- as.data.frame(do.call(rbind, strsplit(paths, ";", fixed = TRUE)),
                            stringsAsFactors = FALSE)
  names(taxonomy) <- c("superkingdom", "kingdom", "phylum", "class", "order",
                       "family", "genus", "species")
  taxonomy$TaxonPath <- paths
  counts <- matrix(c(100, 20, 10, 20, 100, 0, 0, 0, 0), nrow = 3, byrow = TRUE,
                   dimnames = list(NULL, samples))
  context <- list(
    config = get_default_config(), samples = samples, count_matrix = counts,
    taxonomy = taxonomy, unclass_index = 1L,
    sample_stats = data.frame(
      SampleID = samples, TotalReads = colSums(counts),
      UnclassifiedReads = counts[1, ], ClassifiedReads = colSums(counts[-1, , drop = FALSE]),
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(SampleID = samples, Group = c("B", "A", "B"),
                          stringsAsFactors = FALSE)
  )
  context$config$mode <- "cohort"
  context$config$input$metadata <- "metadata.tsv"
  context$config$output$base_dir <- file.path(root, "out")
  context$config$output$dirs <- list(composition = file.path(root, "out", "04_Taxa_Composition"))

  result <- run_taxa_composition(context)
  expect_equal(result$status, "completed")
  expected <- c(
    sprintf("04_%s_stacked.png", c("phylum", "family", "genus")),
    sprintf("04_%s_stacked.tsv", c("phylum", "family", "genus")),
    sprintf("04_%s_group_mean_stacked.png", c("phylum", "family", "genus")),
    sprintf("04_%s_group_mean_stacked.tsv", c("phylum", "family", "genus")),
    sprintf("04_heatmap_%s.png", c("phylum", "family", "genus")),
    sprintf("04_heatmap_%s.tsv", c("phylum", "family", "genus"))
  )
  expect_true(all(file.exists(file.path(context$config$output$dirs$composition, expected))))
  stacked <- read.delim(file.path(context$config$output$dirs$composition, "04_phylum_stacked.tsv"),
                        check.names = FALSE)
  expect_identical(names(stacked), c(
    "Rank", "TaxonPath", "DisplayTaxon", "SampleID", "Group", "RelativeAbundance",
    "MeanRelativeAbundance", "IsOther", "ValidDenominator", "StackOrder", "SampleOrder"
  ))
  expect_true(any(stacked$SampleID == "S1" & !stacked$ValidDenominator))
})
