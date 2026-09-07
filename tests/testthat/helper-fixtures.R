# =============================================================================
# Testthat Helper: Synthetic Fixtures & Generators
# =============================================================================

create_temp_abundance <- function(dir, n_species = 10, sample_names = c("Sample1"),
                                  include_total = TRUE, unclass_reads = 100) {
  ranks_template <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp%02d"
  lineages <- vapply(seq_len(n_species), function(i) sprintf(ranks_template, i), character(1))

  unclass_lineage <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
  all_lineages <- c(unclass_lineage, lineages)

  set.seed(123)
  count_df <- data.frame(tax = all_lineages, stringsAsFactors = FALSE)

  for (idx in seq_along(sample_names)) {
    s <- sample_names[idx]
    base_lambda <- if (idx %% 2 == 1) 40 else 80
    sp_counts <- rpois(n_species, lambda = base_lambda)
    # Zero out a few species per sample to introduce presence/absence differences
    zero_idx <- unique(c((idx * 2) %% n_species + 1, (idx * 3) %% n_species + 1))
    sp_counts[zero_idx] <- 0L
    counts <- c(unclass_reads, sp_counts)
    count_df[[s]] <- counts
  }

  if (include_total) {
    count_df$total <- rowSums(as.matrix(count_df[, sample_names, drop = FALSE]))
  }

  file_path <- file.path(dir, "synthetic_abundance.tsv")
  write.table(count_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  file_path
}

create_temp_assignments <- function(dir, sample_id = "Sample1", n_classified = 50, n_unclassified = 10) {
  total <- n_classified + n_unclassified
  read_ids <- sprintf("read_%05d", seq_len(total))

  status <- c(rep("C", n_classified), rep("U", n_unclassified))
  taxids <- c(rep(1386L, n_classified), rep(0L, n_unclassified))
  len_fields <- c(
    sprintf("0|%d", sample(1400:1600, n_classified, replace = TRUE)),
    sprintf("%d", sample(1200:1500, n_unclassified, replace = TRUE))
  )
  lineages <- c(
    rep("Bacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus cereus", n_classified),
    rep("Unclassified", n_unclassified)
  )

  df <- data.frame(
    status = status,
    read_id = read_ids,
    taxid = taxids,
    len_field = len_fields,
    lineage = lineages,
    stringsAsFactors = FALSE
  )

  file_path <- file.path(dir, sprintf("%s_assignments.tsv", sample_id))
  write.table(df, file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  file_path
}

create_temp_metadata <- function(dir, sample_names = c("Sample1", "Sample2"), groups = c("Control", "Treated")) {
  df <- data.frame(
    SampleID = sample_names,
    Group = groups,
    stringsAsFactors = FALSE
  )
  file_path <- file.path(dir, "metadata.tsv")
  write.table(df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  file_path
}

create_temp_params <- function(dir, classifier = "minimap2",
                               database_set = "ncbi_16s_18s", taxonomic_rank = "S") {
  path <- file.path(dir, "params.json")
  database_sets <- list()
  database_sets[[database_set]] <- list(
    reference = sprintf("s3://ont-open-data/workflow-databases/wf-metagenomics-dbs/%s/%s.fna",
                        database_set, if (database_set == "ncbi_16s_18s") "ncbi_targeted_loci_16s_18s" else "ncbi_16s_18s_28s_ITS"),
    database = sprintf("s3://ont-open-data/workflow-databases/wf-metagenomics-dbs/%s/%s_kraken2.tar.gz",
                       database_set, if (database_set == "ncbi_16s_18s") "ncbi_targeted_loci" else "ncbi_16s_18s_28s_ITS"),
    ref2taxid = if (database_set == "ncbi_16s_18s")
      "s3://ont-open-data/workflow-databases/wf-metagenomics-dbs/ncbi_16s_18s/ref2taxid.targloci.tsv"
      else "s3://ont-open-data/workflow-databases/wf-metagenomics-dbs/ncbi_16s_18s_28s_ITS/ref2taxid.ncbi_16s_18s_28s_ITS.tsv",
    taxonomy = sprintf("s3://ont-open-data/workflow-databases/wf-metagenomics-dbs/%s/new_taxdump_2025-01-01.zip", database_set)
  )
  jsonlite::write_json(list(
    classifier = classifier, database_set = database_set,
    taxonomic_rank = taxonomic_rank, min_len = 1300, max_len = 1700,
    min_read_qual = 10, min_percent_identity = 90, min_ref_coverage = 90,
    abundance_threshold = 1, taxonomy = NULL, reference = NULL,
    ref2taxid = NULL, database = NULL, database_sets = database_sets,
    output_unclassified = TRUE, include_read_assignments = TRUE,
    wf = list(agent = "epi2melabs/test")
  ), path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  path
}

create_temp_bamstats <- function(dir, sample_id, data) {
  path <- file.path(dir, "bamstats.readstats.tsv.gz")
  data$sample_name <- sample_id
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(data[, c("name", "sample_name", "iden", "ref_coverage")], con,
              sep = "\t", row.names = FALSE, quote = FALSE, na = "nan")
  path
}
