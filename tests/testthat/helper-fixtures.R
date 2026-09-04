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
