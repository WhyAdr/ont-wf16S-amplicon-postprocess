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

test_that("display labels are invariant to input-row permutation and literal Other", {
  table <- data.frame(
    TaxonPath = c("Bacteria;P1;Shared", "Bacteria;P2;Shared",
                  "Bacteria;P3;Other", "Bacteria;P4;Other [P3]"),
    stringsAsFactors = FALSE
  )
  expected <- make_taxon_display_map(table)
  set.seed(1046)
  for (i in seq_len(100L)) {
    permuted <- table[sample(seq_len(nrow(table))), , drop = FALSE]
    actual <- make_taxon_display_map(permuted)
    expect_equal(actual[order(actual$TaxonPath), , drop = FALSE],
                 expected[order(expected$TaxonPath), , drop = FALSE])
  }
  expect_length(unique(expected$DisplayTaxon), nrow(expected))
  expect_false(any(expected$DisplayTaxon == "Other"))
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
  expect_equal(unique(prepared$sidecar$PseudoCount), 0.125, tolerance = 1e-12)
  expect_equal(sum(prepared$sidecar$RelativeAbundance[prepared$sidecar$SampleID == "S1"]), 0)
  expect_true(all(is.na(prepared$matrix_plot[, "S1"])))
  expect_true(all(is.finite(prepared$matrix_transformed[, "S1"])))
  expect_false(any(prepared$sidecar$TaxonPath == "__OTHER__" &
                   prepared$sidecar$SampleID == "S1" &
                   prepared$sidecar$RelativeAbundance > 0))

  none <- prepare_heatmap_data(
    table, c("S2", "S10", "S1"), metadata, top_n = 1L,
    include_other = FALSE, transform = "none", rank = "genus"
  )
  expect_true(all(is.na(none$sidecar$PseudoCount)))
  expect_identical(none$sidecar$TransformedValue, none$sidecar$RelativeAbundance)
})

test_that("invalid heatmap denominators skip the complete heatmap", {
  table <- composition_rank_fixture()
  metadata <- data.frame(
    SampleID = "S1", Group = "empty", ValidDenominator = FALSE,
    stringsAsFactors = FALSE
  )
  prepared <- prepare_heatmap_data(
    table[, c("TaxonPath", "Taxon", "S1"), drop = FALSE], "S1", metadata,
    top_n = 2L, include_other = TRUE, transform = "log10_relative", rank = "genus"
  )
  expect_true(prepared$skipped)
  expect_equal(prepared$skip$Rank, "genus")
})

test_that("a single zero-classified module run emits only structured composition skips", {
  root <- tempfile("single_zero_composition_")
  paths <- c(
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis"
  )
  taxonomy <- as.data.frame(do.call(rbind, strsplit(paths, ";", fixed = TRUE)),
                            stringsAsFactors = FALSE)
  names(taxonomy) <- c("superkingdom", "kingdom", "phylum", "class", "order",
                       "family", "genus", "species")
  taxonomy$TaxonPath <- paths
  counts <- matrix(c(4, 0), ncol = 1, dimnames = list(NULL, "S_zero"))
  cfg <- get_default_config()
  cfg$output$base_dir <- root
  cfg$output$dirs <- list(composition = file.path(root, "04_Taxa_Composition"))
  context <- list(
    config = cfg, samples = "S_zero", count_matrix = counts,
    taxonomy = taxonomy, unclass_index = 1L,
    sample_stats = data.frame(
      SampleID = "S_zero", TotalReads = 4, UnclassifiedReads = 4,
      ClassifiedReads = 0, stringsAsFactors = FALSE
    ),
    metadata = NULL
  )
  result <- run_taxa_composition(context)
  expect_equal(result$status, "completed")
  expect_true(file.exists(file.path(
    cfg$output$dirs$composition, "04_legacy_single_sample_skipped.tsv"
  )))
  expect_false(any(file.exists(file.path(
    cfg$output$dirs$composition,
    c("04a_phylum_composition.png", "04b_family_composition.png",
      "04c_genus_composition.png", "04d_species_composition.png")
  ))))
  expect_true(all(file.exists(file.path(
    cfg$output$dirs$composition,
    c("04_phylum_stacked_skipped.tsv", "04_family_stacked_skipped.tsv",
      "04_genus_stacked_skipped.tsv", "04_heatmap_phylum_skipped.tsv",
      "04_heatmap_family_skipped.tsv", "04_heatmap_genus_skipped.tsv")
  ))))
})

test_that("composition validators reject inconsistent identity and order metadata", {
  table <- composition_rank_fixture()
  samples <- c("S2", "S10", "S1")
  metadata <- data.frame(
    SampleID = samples, Group = c("B", "A", "B"),
    ValidDenominator = c(TRUE, TRUE, FALSE), stringsAsFactors = FALSE
  )
  collapsed <- collapse_rank_for_display(
    table, c("Bacteria;P1;Shared", "Bacteria;P2;Shared"), samples, metadata
  )
  collapsed$Rank <- "genus"
  stacked <- collapsed[, c("Rank", composition_stacked_columns[composition_stacked_columns != "Rank"])]
  expected_sample_order <- c("S2", "S1", "S10")
  expect_silent(validate_stacked_sidecar(
    stacked, samples, "genus", expected_sample_order
  ))

  bad <- stacked
  bad$SampleOrder[bad$SampleID == "S2"] <- 3L
  bad$SampleOrder[bad$SampleID == "S10"] <- 1L
  expect_error(validate_stacked_sidecar(
    bad, samples, "genus", expected_sample_order
  ), "expected deterministic sample order")

  bad <- stacked
  s2 <- which(bad$SampleID == "S2")
  bad$ValidDenominator[s2[length(s2)]] <- FALSE
  expect_error(validate_stacked_sidecar(
    bad, samples, "genus", expected_sample_order
  ), "inconsistent within a sample")

  bad <- stacked
  s10 <- which(bad$SampleID == "S10")
  bad$StackOrder[s10] <- rev(bad$StackOrder[s10])
  expect_error(validate_stacked_sidecar(
    bad, samples, "genus", expected_sample_order
  ), "mapping differs across samples")

  groups <- summarize_group_means(collapsed)
  groups$Rank <- "genus"
  groups <- groups[, composition_group_columns]
  expect_silent(validate_group_sidecar(groups, collapsed, "genus"))
  bad_samples <- collapsed
  bad_samples$GroupOrder[bad_samples$SampleID == "S1"] <- 2L
  expect_error(validate_group_sidecar(groups, bad_samples, "genus"),
               "conflicting GroupOrder")
  expect_error(
    validate_group_sidecar(groups[groups$Group == "B", ], collapsed, "genus"),
    "does not cover"
  )
  bad <- groups
  bad$GroupOrder[bad$Group == "B"] <- 99L
  expect_error(validate_group_sidecar(bad, collapsed, "genus"),
               "order or TaxonPath coverage")

  cfg <- get_default_config()
  prepared <- prepare_heatmap_data(
    table, samples, metadata, top_n = 2L, include_other = TRUE,
    transform = "log10_relative", rank = "genus"
  )
  heatmap <- prepared$sidecar
  expect_silent(validate_heatmap_sidecar(
    heatmap, samples, cfg, "genus", prepared$sample_order
  ))

  bad <- heatmap
  bad$ColumnOrder[bad$SampleID == "S2"] <- 3L
  bad$ColumnOrder[bad$SampleID == "S10"] <- 1L
  expect_error(validate_heatmap_sidecar(
    bad, samples, cfg, "genus", prepared$sample_order
  ), "expected deterministic sample order")

  bad <- heatmap
  s1 <- which(bad$SampleID == "S1")
  bad$ValidDenominator[s1[length(s1)]] <- TRUE
  expect_error(validate_heatmap_sidecar(
    bad, samples, cfg, "genus", prepared$sample_order
  ), "inconsistent within a sample")

  bad <- heatmap
  s10 <- which(bad$SampleID == "S10")
  bad$RowOrder[s10] <- rev(bad$RowOrder[s10])
  expect_error(validate_heatmap_sidecar(
    bad, samples, cfg, "genus", prepared$sample_order
  ), "mapping differs across samples")
})

test_that("zero-valid groups retain NA means and an explicit invalid plot state", {
  table <- composition_rank_fixture()
  metadata <- data.frame(
    SampleID = c("S2", "S10", "S1"),
    Group = c("valid", "valid", "empty"),
    ValidDenominator = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  collapsed <- collapse_rank_for_display(
    table, c("Bacteria;P1;Shared", "Bacteria;P2;Shared"),
    metadata$SampleID, metadata
  )
  means <- summarize_group_means(collapsed)
  expect_true(all(is.na(means$MeanRelativeAbundance[means$Group == "empty"])))
  expect_equal(unique(means$SamplesUsed[means$Group == "empty"]), 0L)

  plot_data <- means
  plot_data$SampleID <- plot_data$Group
  plot_data$RelativeAbundance <- plot_data$MeanRelativeAbundance
  plot_data$ValidDenominator <- plot_data$SamplesUsed > 0L
  plot_data$SampleOrder <- plot_data$GroupOrder
  plot <- build_stacked_taxa_plot(
    plot_data[, c("TaxonPath", "DisplayTaxon", "SampleID", "Group",
                  "RelativeAbundance", "IsOther", "ValidDenominator",
                  "StackOrder", "GroupOrder")],
    composition_colors(unique(plot_data$TaxonPath), unique(plot_data$DisplayTaxon)),
    c("valid", "empty"), c("valid", "empty"), single = TRUE,
    rank = "genus group mean",
    subtitle = "Arithmetic mean of per-sample classified-read relative abundance; samples are not pooled by read depth.",
    invalid_groups = "empty"
  )
  expect_match(plot$labels$subtitle, "arithmetic mean", ignore.case = TRUE)
  expect_match(plot$labels$subtitle, "not pooled", ignore.case = TRUE)
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))))
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
  expect_equal(sum(stacked$RelativeAbundance[stacked$SampleID == "S1"]), 0)
  expect_false(any(stacked$TaxonPath == "__OTHER__" & stacked$SampleID == "S1" &
                   stacked$RelativeAbundance > 0))
})
