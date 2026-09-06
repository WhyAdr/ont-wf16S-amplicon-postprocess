# =============================================================================
# Unit Tests: Shared Data Layer & I/O Validation
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))

test_that("producer contract accepts supported minimap2 and rejects unsafe alternatives", {
  root <- tempfile("params_contract_")
  dir.create(root)
  supported <- read_upstream_params(create_temp_params(root))
  expect_equal(supported$classifier, "minimap2")
  expect_error(read_upstream_params(create_temp_params(root, classifier = "kraken2")),
               "supports minimap2 only")
  expect_error(read_upstream_params(create_temp_params(root, database_set = "SILVA_138_1")),
               "bundled NCBI database sets only")
  expect_error(read_upstream_params(create_temp_params(root, taxonomic_rank = "G")),
               "expected species rank")
})

test_that("bamstats partition is one-to-one and conserves all C0 reads", {
  root <- tempfile("bamstats_contract_")
  dir.create(root)
  reads <- data.frame(
    status = c("C", "C", "C", "C"),
    read_id = c("positive", "identity", "coverage", "both"),
    taxid = c(123, 0, 0, 0), stringsAsFactors = FALSE
  )
  bamstats <- create_temp_bamstats(root, "S1", data.frame(
    name = reads$read_id, iden = c(95, 89, 95, 89),
    ref_coverage = c(95, 95, 89, 89)
  ))
  params <- read_upstream_params(create_temp_params(root))
  observed <- partition_minimap2_failures(reads, bamstats, params, "S1")
  expect_equal(unlist(observed), c(matched = 3, identity_only = 1,
                                   coverage_only = 1, both = 1))
})

test_that("Synthetic abundance table parses and validates correctly", {
  tmp <- tempdir()
  ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = c("S1", "S2"))

  res <- read_abundance_table(ab_file)
  expect_equal(length(res$samples), 2)
  expect_equal(res$samples, c("S1", "S2"))
  expect_equal(nrow(res$count_matrix), 6) # 1 unclass + 5 species
  expect_equal(ncol(res$taxonomy), 9) # 8 ranks + TaxonPath
})

test_that("Abundance table errors if total column does not equal row sum", {
  tmp <- tempdir()
  ab_file <- create_temp_abundance(tmp, n_species = 3, sample_names = c("S1"))

  # Corrupt total column
  df <- read.delim(ab_file, check.names = FALSE)
  df$total[1] <- df$total[1] + 999
  write.table(df, ab_file, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_error(read_abundance_table(ab_file), "does not equal sample row sums")
})

test_that("Abundance table errors if rank count is not 8", {
  tmp <- tempdir()
  file_path <- file.path(tmp, "bad_rank_ab.tsv")
  bad_df <- data.frame(
    tax = c("Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown",
            "Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus"), # 6 ranks
    S1 = c(10, 50),
    total = c(10, 50)
  )
  write.table(bad_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_error(read_abundance_table(file_path), "expected 8 ranks, found 6")
})

test_that("Assignments parser handles pipe and plain lengths and reconciles counts", {
  tmp <- tempdir()
  asgn_file <- create_temp_assignments(tmp, sample_id = "S1", n_classified = 40, n_unclassified = 10)

  reads <- read_assignments_file(
    asgn_file,
    sample_id = "S1",
    expected_total = 50,
    expected_classified = 40,
    expected_unclassified = 10
  )

  expect_equal(nrow(reads), 50)
  expect_equal(sum(reads$effective_classified), 40)
  expect_true(all(reads$read_length >= 1200 & reads$read_length <= 1600))
  expect_true(all(is.integer(reads$read_length)))
})

test_that("Assignments parser rejects non-five-field rows with a physical line number", {
  tmp <- tempfile("bad_assignment_")
  dir.create(tmp)
  bad_file <- file.path(tmp, "bad.tsv")
  writeLines(c("C\tread_1\t123\t0|1500\tBacteria|Example", "U\tread_2\t0\t1400"), bad_file)
  expect_error(read_assignments_file(bad_file, "S1"), "4 fields at line 2; expected exactly 5")
})

test_that("Real Ambar Ayunda fixture satisfies all Section 2.2 invariants", {
  ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
  asgn_path <- file.path("..", "..", "output_AAy", "reads_assignments",
                         "AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv")

  skip_if_not(file.exists(ab_path), "Real abundance table not found")
  skip_if_not(file.exists(asgn_path), "Real assignments file not found")

  ab_res <- read_abundance_table(ab_path)
  expect_equal(nrow(ab_res$count_matrix), 1837)
  expect_equal(ab_res$samples, "AmbarAyunda_minimap2_16S")

  sample_col <- ab_res$samples[1]
  total_reads <- sum(ab_res$count_matrix[, sample_col])
  unclass_reads <- ab_res$count_matrix[ab_res$unclass_index, sample_col]
  classified_reads <- total_reads - unclass_reads

  expect_equal(total_reads, 114056)
  expect_equal(classified_reads, 80556)
  expect_equal(unclass_reads, 33500)

  # Check assignments
  reads <- read_assignments_file(
    asgn_path,
    sample_id = sample_col,
    expected_total = 114056,
    expected_classified = 80556,
    expected_unclassified = 33500
  )

  expect_equal(nrow(reads), 114056)
  expect_equal(sum(reads$status == "C"), 89809)
  expect_equal(sum(reads$status == "U"), 24247)
  expect_equal(sum(reads$status == "C" & reads$taxid == 0), 9253)
  expect_equal(sum(reads$effective_classified), 80556)
})

test_that("Metadata validation aligns samples and detects discrepancies", {
  tmp <- tempdir()
  meta_file <- create_temp_metadata(tmp, sample_names = c("S2", "S1"), groups = c("GroupB", "GroupA"))

  meta_aligned <- read_metadata_table(meta_file, selected_samples = c("S1", "S2"))
  expect_equal(meta_aligned$SampleID, c("S1", "S2"))
  expect_equal(meta_aligned$Group, c("GroupA", "GroupB"))

  # Test missing sample
  expect_error(
    read_metadata_table(meta_file, selected_samples = c("S1", "S2", "S3")),
    "Missing from metadata: S3"
  )
})

test_that("Metadata rejects empty groups", {
  tmp <- tempfile("bad_metadata_")
  dir.create(tmp)
  meta_file <- file.path(tmp, "metadata.tsv")
  writeLines(c("SampleID\tGroup", "S1\t"), meta_file)
  expect_error(read_metadata_table(meta_file, "S1"), "empty Group")
})

test_that("Sample IDs reject unsafe names and post-sanitization collisions", {
  invalid_sets <- list(
    c("S1", "bad/name"),
    c("S1", "bad\\name"),
    c("S1", "."),
    c("S1", ".."),
    c("S1", "bad\nname"),
    c("A B", "A?B")
  )
  for (sample_ids in invalid_sets) {
    expect_error(validate_sample_ids(sample_ids), "Sample ID validation error")
  }
  expect_error(validate_sample_ids(c("S1", "")), "Empty or NA")
})
