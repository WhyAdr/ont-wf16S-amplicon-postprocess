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

test_that("producer contract requires whole-number length and abundance thresholds", {
  root <- tempfile("params_integer_contract_")
  dir.create(root)
  for (field in c("min_len", "max_len", "abundance_threshold")) {
    path <- create_temp_params(root)
    params <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    params[[field]] <- params[[field]] + 0.5
    jsonlite::write_json(params, path, auto_unbox = TRUE, null = "null")
    expect_error(read_upstream_params(path), sprintf("field '%s' must be a whole number", field))
  }
  params <- read_upstream_params(create_temp_params(root))
  expect_equal(params$min_len, 1300)
  path <- create_temp_params(root)
  one_point_zero <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  one_point_zero$abundance_threshold <- 1.0
  jsonlite::write_json(one_point_zero, path, auto_unbox = TRUE, null = "null")
  expect_no_error(read_upstream_params(path))
})

test_that("bamstats partition is one-to-one and conserves all C0 reads", {
  root <- tempfile("bamstats_contract_")
  dir.create(root)
  reads <- data.frame(
    status = c("C", "C", "C", "C", "U"),
    read_id = c("positive", "identity", "coverage", "both", "raw_u"),
    taxid = c(123, 0, 0, 0, 0), stringsAsFactors = FALSE
  )
  bamstats <- create_temp_bamstats(root, "S1", data.frame(
    name = reads$read_id, iden = c(95, 89, 95, 89, NaN),
    ref_coverage = c(95, 95, 89, 89, NaN)
  ))
  params <- read_upstream_params(create_temp_params(root))
  observed <- partition_minimap2_failures(reads, bamstats, params, "S1")
  expect_equal(unlist(observed), c(matched = 3, identity_only = 1,
                                   coverage_only = 1, both = 1))
})

test_that("bamstats partition rejects missing or non-finite values for status-C reads", {
  root <- tempfile("bamstats_contract_c_nan_")
  dir.create(root)
  reads <- data.frame(
    status = c("C"),
    read_id = c("corrupt"),
    taxid = c(0), stringsAsFactors = FALSE
  )
  bamstats <- create_temp_bamstats(root, "S1", data.frame(
    name = reads$read_id, iden = c(NaN),
    ref_coverage = c(95)
  ))
  params <- read_upstream_params(create_temp_params(root))
  expect_error(partition_minimap2_failures(reads, bamstats, params, "S1"),
               "must be finite numbers")
})

test_that("bamstats partition rejects NA, empty, or mismatched sample_name", {
  root <- tempfile("bamstats_contract_sample_na_")
  dir.create(root)
  reads <- data.frame(status = c("C"), read_id = c("r1"), taxid = c(123), stringsAsFactors = FALSE)
  con <- gzfile(file.path(root, "bamstats.readstats.tsv.gz"), open = "wt")
  writeLines(c("name\tsample_name\tiden\tref_coverage", "r1\t\t95\t95"), con)
  close(con)
  params <- read_upstream_params(create_temp_params(root))
  expect_error(
    partition_minimap2_failures(reads, file.path(root, "bamstats.readstats.tsv.gz"), params, "S1"),
    "consistently equal"
  )
})

test_that("discover_bamstats_files rejects missing or inconsistent sample_names", {
  root <- tempfile("bamstats_discover_bad_")
  sub_dir <- file.path(root, "sample.bamstats_results")
  dir.create(sub_dir, recursive = TRUE)
  con <- gzfile(file.path(sub_dir, "bamstats.readstats.tsv.gz"), open = "wt")
  writeLines(c("name\tsample_name\tiden\tref_coverage", "r1\tS1\t95\t95", "r2\tS2\t95\t95"), con)
  close(con)
  expect_error(discover_bamstats(root, c("S1", "S2")), "inconsistent sample names")
})

test_that("discover_bamstats rejects padded sample names", {
  root <- tempfile("bamstats_discover_padded_")
  sub_dir <- file.path(root, "sample.bamstats_results")
  dir.create(sub_dir, recursive = TRUE)
  con <- gzfile(file.path(sub_dir, "bamstats.readstats.tsv.gz"), open = "wt")
  writeLines(c("name\tsample_name\tiden\tref_coverage", "r1\tS1 \t95\t95"), con)
  close(con)
  expect_error(discover_bamstats(root, "S1"), "leading/trailing whitespace")
})

test_that("discover_bamstats rejects multiple bamstats files mapping to the same sample", {
  root <- tempfile("bamstats_duplicate_")
  dir.create(root)
  dir1 <- file.path(root, "run1.bamstats_results")
  dir2 <- file.path(root, "run2.bamstats_results")
  dir.create(dir1)
  dir.create(dir2)
  create_temp_bamstats(dir1, "S1", data.frame(name = "r1", iden = 95, ref_coverage = 95))
  create_temp_bamstats(dir2, "S1", data.frame(name = "r2", iden = 95, ref_coverage = 95))
  expect_error(discover_bamstats(root, "S1"), "Multiple bamstats files discovered for sample 'S1'")
})

test_that("bamstats preserves numeric-looking read and sample identifiers", {
  root <- tempfile("bamstats_numeric_identifiers_")
  dir.create(root)
  bams_dir <- file.path(root, "numeric.bamstats_results")
  dir.create(bams_dir)
  bamstats <- create_temp_bamstats(
    bams_dir, "01", data.frame(name = "0001", iden = 89, ref_coverage = 95)
  )
  reads <- data.frame(status = "C", read_id = "0001", taxid = 0,
                      stringsAsFactors = FALSE)
  params <- read_upstream_params(create_temp_params(root))
  discovered <- discover_bamstats(root, "01")
  expect_identical(basename(discovered[["01"]]), basename(bamstats))
  expect_equal(
    unlist(partition_minimap2_failures(reads, discovered[["01"]], params, "01")),
    c(matched = 1, identity_only = 1, coverage_only = 0, both = 0)
  )
})

test_that("bamstats helper preserves ordinary alphanumeric identifiers", {
  root <- tempfile("bamstats_alphanumeric_identifiers_")
  dir.create(root)
  path <- create_temp_bamstats(root, "Sample_A", data.frame(
    name = "read_0001", iden = 95, ref_coverage = 95
  ))
  parsed <- read_bamstats_table(path)
  expect_identical(parsed$name, "read_0001")
  expect_identical(parsed$sample_name, "Sample_A")
})

test_that("partition_minimap2_failures rejects duplicate read names within a single bamstats file", {
  root <- tempfile("bamstats_duplicate_name_")
  dir.create(root)
  reads <- data.frame(status = c("C"), read_id = c("r1"), taxid = c(123), stringsAsFactors = FALSE)
  con <- gzfile(file.path(root, "bamstats.readstats.tsv.gz"), open = "wt")
  writeLines(c("name\tsample_name\tiden\tref_coverage", "r1\tS1\t95\t95", "r1\tS1\t95\t95"), con)
  close(con)
  params <- read_upstream_params(create_temp_params(root))
  expect_error(
    partition_minimap2_failures(reads, file.path(root, "bamstats.readstats.tsv.gz"), params, "S1"),
    "Bamstats read names for 'S1' must be non-empty and unique"
  )
})

test_that("partition_minimap2_failures rejects bamstats missing status-C reads", {
  root <- tempfile("bamstats_missing_c_")
  dir.create(root)
  reads <- data.frame(
    status = c("C", "C"),
    read_id = c("read_present", "read_absent"),
    taxid = c(123, 0), stringsAsFactors = FALSE
  )
  bamstats <- create_temp_bamstats(root, "S1", data.frame(
    name = c("read_present"), iden = c(95), ref_coverage = c(95)
  ))
  params <- read_upstream_params(create_temp_params(root))
  expect_error(partition_minimap2_failures(reads, bamstats, params, "S1"),
               "Bamstats is missing 1 status-C read[(]s[)] for 'S1'")
})

test_that("partition_minimap2_failures rejects non-numeric or malformed metrics for status-C reads", {
  root <- tempfile("bamstats_malformed_num_")
  dir.create(root)
  reads <- data.frame(
    status = c("C"),
    read_id = c("r1"),
    taxid = c(0), stringsAsFactors = FALSE
  )
  con <- gzfile(file.path(root, "bamstats.readstats.tsv.gz"), open = "wt")
  writeLines(c("name\tsample_name\tiden\tref_coverage", "r1\tS1\tinvalid_string\t95"), con)
  close(con)
  params <- read_upstream_params(create_temp_params(root))
  expect_error(partition_minimap2_failures(reads, file.path(root, "bamstats.readstats.tsv.gz"), params, "S1"),
               "must be finite numbers")
})

test_that("partition_minimap2_failures rejects C+TaxID0 reads passing both thresholds", {
  root <- tempfile("bamstats_c0_passes_")
  dir.create(root)
  reads <- data.frame(
    status = c("C"),
    read_id = c("r1"),
    taxid = c(0), stringsAsFactors = FALSE
  )
  bamstats <- create_temp_bamstats(root, "S1", data.frame(
    name = c("r1"), iden = c(95), ref_coverage = c(95)
  ))
  params <- read_upstream_params(create_temp_params(root))
  expect_error(partition_minimap2_failures(reads, bamstats, params, "S1"),
               "At least one C[+]TaxID0 read for 'S1' passes both recorded thresholds")
})

test_that("partition_minimap2_failures rejects TaxID>0 reads failing recorded thresholds", {
  root <- tempfile("bamstats_positive_fails_")
  dir.create(root)
  reads <- data.frame(
    status = c("C"),
    read_id = c("r1"),
    taxid = c(123), stringsAsFactors = FALSE
  )
  bamstats <- create_temp_bamstats(root, "S1", data.frame(
    name = c("r1"), iden = c(85), ref_coverage = c(95)
  ))
  params <- read_upstream_params(create_temp_params(root))
  expect_error(partition_minimap2_failures(reads, bamstats, params, "S1"),
               "At least one TaxID>0 read for 'S1' fails the recorded thresholds")
})

test_that("QC module reports BamstatsAvailable TRUE when bamstats is mapped without assignments", {
  source(file.path("..", "..", "analysis", "01_qc_diagnostics.R"))
  root <- tempfile("qc_bamstats_no_asgn_")
  dir.create(root)
  bams_dir <- file.path(root, "sample.bamstats_results")
  dir.create(bams_dir)
  create_temp_bamstats(bams_dir, "S1", data.frame(name = "r1", iden = 95, ref_coverage = 95))

  cfg <- get_default_config()
  cfg$input$abundance_table <- create_temp_abundance(root, n_species = 3, sample_names = "S1")
  cfg$input$params_json <- create_temp_params(root)
  cfg$input$wf16s_output_root <- root
  cfg$input$assignments <- NULL
  cfg$output$base_dir <- file.path(root, "out")
  cfg$output$dirs <- list(qc = file.path(cfg$output$base_dir, "01_QC"))

  context <- build_context(cfg)
  expect_true(file.exists(context$bamstats[["S1"]]))

  res_qc <- run_qc(context)
  expect_equal(res_qc$status, "completed")

  inv <- read.delim(file.path(cfg$output$dirs$qc, "00_read_investigation.tsv"), check.names = FALSE)
  expect_equal(nrow(inv), 1L)
  expect_false(inv$AssignmentAvailable)
  expect_true(inv$BamstatsAvailable)
  expect_true(is.na(inv$BamstatsC0Matched))
  expect_true(is.na(inv$IdentityOnlyFailed))
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

test_that("Abundance table rejects empty or padded rank cells", {
  root <- tempfile("bad_rank_cells_")
  dir.create(root)
  bad_lineages <- c(
    "Bacteria;;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Species",
    "Bacteria; ;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Species",
    "Bacteria;Kingdom ;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Species"
  )
  for (lineage in bad_lineages) {
    path <- create_temp_abundance(root, n_species = 1, sample_names = "S1")
    table <- read.delim(path, check.names = FALSE)
    table$tax[2] <- lineage
    write.table(table, path, sep = "\t", row.names = FALSE, quote = FALSE)
    expect_error(read_abundance_table(path), "Lineage schema violation at row 3: rank")
  }
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

test_that("Assignments parser streams plain and gzipped files across chunk boundaries", {
  root <- tempfile("streamed_assignments_")
  dir.create(root)
  lines <- c(
    "C\tread_001\t123\t0|1501\tBacteria|Example",
    "",
    "U\tread_002\t0\t1402\tUnclassified",
    "C\tread_003\t0\t1|1303\tUnclassified",
    "C\tread_004\t456\t1504\tBacteria|Example"
  )
  plain_path <- file.path(root, "assignments.tsv")
  gz_path <- paste0(plain_path, ".gz")
  writeLines(lines, plain_path)
  gz_connection <- gzfile(gz_path, open = "wt")
  writeLines(lines, gz_connection)
  close(gz_connection)

  plain <- read_assignments_file(plain_path, "S1", expected_total = 4L, chunk_size = 2L)
  gzipped <- read_assignments_file(gz_path, "S1", expected_total = 4L, chunk_size = 3L)

  expect_identical(plain, gzipped)
  expect_identical(names(plain), c(
    "status", "read_id", "taxid", "len_field", "lineage", "read_length", "effective_classified"
  ))
  expect_type(plain$read_id, "character")
  expect_type(plain$taxid, "character")
  expect_type(plain$read_length, "integer")
  expect_identical(plain$read_length, c(1501L, 1402L, 1303L, 1504L))
  expect_identical(plain$effective_classified, c(TRUE, FALSE, FALSE, TRUE))
})

test_that("Assignments parser retains validation contracts while streaming", {
  root <- tempfile("streamed_assignment_validation_")
  dir.create(root)
  duplicate_path <- file.path(root, "duplicate.tsv")
  empty_id_path <- file.path(root, "empty_id.tsv")
  trailing_field_path <- file.path(root, "trailing_field.tsv")
  writeLines(c(
    "C\tread_001\t123\t1500\tBacteria",
    "U\tread_001\t0\t1400\tUnclassified"
  ), duplicate_path)
  writeLines("C\t\t123\t1500\tBacteria", empty_id_path)
  writeLines("C\tread_001\t123\t1500\t", trailing_field_path)

  expect_error(read_assignments_file(duplicate_path, "S1", chunk_size = 1L),
               "contains duplicate read ID: 'read_001'")
  expect_error(read_assignments_file(empty_id_path, "S1"), "unsafe read ID at line 1")
  expect_error(read_assignments_file(trailing_field_path, "S1"), "blank or unsafe lineage")
  expect_error(read_assignments_file(trailing_field_path, "S1", chunk_size = 0L),
               "chunk_size.*positive integer")
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
  expect_equal(sum(reads$status == "C" & reads$taxid == "0"), 9253)
  expect_equal(sum(reads$effective_classified), 80556)
  expect_identical(sum(reads$read_length), 171330353L)

  # This digest covers the complete typed return value, including column order
  # and values, so the streamed reader remains equivalent for the fixture.
  digest_path <- tempfile("assignments_typed_", fileext = ".rds")
  saveRDS(reads, digest_path, version = 2)
  expected_digest <- if (.Platform$OS.type == "windows") {
    "5f3b6f3227754fd0a739e2c5b4626caa"
  } else {
    "58a09628ccfbcb106b352ab247c08613"
  }
  expect_identical(unname(tools::md5sum(digest_path)), expected_digest)
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

test_that("Metadata preserves numeric-looking identity columns as exact strings", {
  root <- tempfile("metadata_identity_")
  dir.create(root)
  path <- file.path(root, "metadata.tsv")
  writeLines(c("SampleID\tGroup", "01\t0", "02\t1"), path)
  metadata <- read_metadata_table(path, c("01", "02"))
  expect_identical(metadata$SampleID, c("01", "02"))
  expect_identical(metadata$Group, c("0", "1"))
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
    c("A B", "A?B"),
    c("S1", " S2"),
    c("S1", "S2 "),
    c("Sample", "sample"),
    c("S1", "ends."),
    c("S1", "CON"),
    c("S1", "NUL.txt"),
    c("S1", "COM1"),
    c("S1", "LPT9.tsv")
  )
  for (sample_ids in invalid_sets) {
    expect_error(validate_sample_ids(sample_ids), "Sample ID validation error")
  }
  expect_error(validate_sample_ids(c("S1", "")), "Empty or NA")
})
