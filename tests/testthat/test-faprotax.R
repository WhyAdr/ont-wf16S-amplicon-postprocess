# =============================================================================
# Unit Tests: FAPROTAX 1.2.12 Integration
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "dependencies.R"))
source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "plotting.R"))
source(file.path("..", "..", "analysis", "08_faprotax.R"))

make_faprotax_context <- function(output_dir, include_prokaryotes = TRUE) {
  unclassified <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
  bacillus <- paste(c(
    "Bacteria", "Bacillati", "Bacillota", "Bacilli", "Bacillales",
    "Bacillaceae", "Bacillus", "Bacillus subtilis"
  ), collapse = ";")
  unmapped <- paste(c(
    "Bacteria", "Testkingdom", "Testphylum", "Testclass", "Testorder",
    "Testfamily", "Testgenus", "Test species"
  ), collapse = ";")
  eukaryote <- paste(c(
    "Eukaryota", "Opisthokonta", "Ascomycota", "Saccharomycetes",
    "Saccharomycetales", "Saccharomycetaceae", "Saccharomyces", "Saccharomyces cerevisiae"
  ), collapse = ";")

  lineages <- if (include_prokaryotes) {
    c(unclassified, bacillus, unmapped, eukaryote)
  } else {
    c(unclassified, eukaryote)
  }
  counts <- if (include_prokaryotes) c(5, 20, 15, 10) else c(5, 20)
  parsed <- do.call(rbind, strsplit(lineages, ";", fixed = TRUE))
  taxonomy <- as.data.frame(parsed, stringsAsFactors = FALSE)
  colnames(taxonomy) <- RANKS_8
  taxonomy$TaxonPath <- lineages
  count_matrix <- matrix(
    counts,
    ncol = 1,
    dimnames = list(lineages, "S1")
  )

  list(
    config = list(
      output = list(dirs = list(faprotax = output_dir)),
      faprotax = list(top_n_functions = 20L)
    ),
    samples = "S1",
    count_matrix = count_matrix,
    taxonomy = taxonomy,
    unclass_index = 1L,
    sample_stats = data.frame(
      SampleID = "S1",
      TotalReads = sum(counts),
      UnclassifiedReads = counts[1],
      ClassifiedReads = sum(counts) - counts[1],
      stringsAsFactors = FALSE
    ),
    metadata = NULL
  )
}

test_that("embedded FAPROTAX database is exactly 1.2.12", {
  runtime <- validate_faprotax_runtime()
  expect_identical(runtime$faprotax_version, "1.2.12")
  expect_true(utils::compareVersion(as.character(utils::packageVersion("microeco")),
                                    "2.3.0") >= 0)
})

test_that("FAPROTAX input maps superkingdom to Kingdom and excludes non-prokaryotes", {
  context <- make_faprotax_context(tempdir())
  prepared <- prepare_faprotax_input(context)

  expect_equal(nrow(prepared$counts), 2L)
  expect_equal(prepared$feature_ids, c("Taxon_000002", "Taxon_000003"))
  expect_equal(prepared$tax_table$Kingdom, c("Bacteria", "Bacteria"))
  expect_false(any(prepared$source_taxonomy$superkingdom == "Eukaryota"))
  expect_false(any(prepared$tax_table$Kingdom == "Bacillati"))
})

test_that("FAPROTAX outputs conserve explicit read denominators", {
  output_dir <- file.path(tempdir(), "faprotax_contract")
  context <- make_faprotax_context(output_dir)
  result <- run_faprotax(context)

  expect_equal(result$status, "completed")
  expect_length(result$outputs, 5L)
  expect_true(all(file.exists(result$outputs)))
  expect_identical(result$faprotax_version, "1.2.12")

  abundance <- read.delim(result$outputs[1], check.names = FALSE, stringsAsFactors = FALSE)
  coverage <- read.delim(result$outputs[2], check.names = FALSE, stringsAsFactors = FALSE)
  assignments <- read.delim(result$outputs[3], check.names = FALSE, stringsAsFactors = FALSE)

  expect_identical(names(abundance), c(
    "SampleID", "Function", "FunctionReadCount",
    "PctEligibleProkaryoticReads", "PctClassifiedReads", "PctTotalReads"
  ))
  expect_identical(names(coverage), c(
    "SampleID", "TotalReads", "UnclassifiedReads", "ClassifiedReads",
    "EligibleProkaryoticReads", "ExcludedNonProkaryoticClassifiedReads",
    "FunctionMappedReads", "FunctionUnmappedEligibleReads",
    "EligibleProkaryoticTaxa", "FunctionMappedTaxa"
  ))
  expect_identical(names(assignments), c("FeatureID", "Function", "TaxonPath"))

  expect_equal(coverage[1, -1], data.frame(
    TotalReads = 50,
    UnclassifiedReads = 5,
    ClassifiedReads = 45,
    EligibleProkaryoticReads = 35,
    ExcludedNonProkaryoticClassifiedReads = 10,
    FunctionMappedReads = 20,
    FunctionUnmappedEligibleReads = 15,
    EligibleProkaryoticTaxa = 2,
    FunctionMappedTaxa = 1,
    check.names = FALSE
  ))
  expect_true(all(abundance$FunctionReadCount <= 35))
  expect_true(all(abundance$PctEligibleProkaryoticReads >= 0 &
                    abundance$PctEligibleProkaryoticReads <= 100))
  expect_true(all(abundance$PctClassifiedReads >= 0 & abundance$PctClassifiedReads <= 100))
  expect_true(all(abundance$PctTotalReads >= 0 & abundance$PctTotalReads <= 100))
  expect_true(all(assignments$FeatureID == "Taxon_000002"))
  expect_true(all(assignments$TaxonPath ==
                    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis"))

  # Multiple functions share the same mapped Bacillus reads. Their percentages
  # are independent and are deliberately not a 100% composition.
  expect_gt(nrow(abundance), 1L)
  expect_gt(sum(abundance$FunctionReadCount), coverage$EligibleProkaryoticReads[1])
})

test_that("FAPROTAX skips cleanly when no positive-count prokaryotes exist", {
  output_dir <- file.path(tempdir(), "faprotax_skip")
  result <- run_faprotax(make_faprotax_context(output_dir, include_prokaryotes = FALSE))

  expect_equal(result$status, "skipped")
  expect_equal(length(result$outputs), 1L)
  expect_true(file.exists(result$outputs))
  skipped <- read.delim(result$outputs, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(skipped$Status, "skipped")
})
