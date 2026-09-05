#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: verify_release_run.R OUTPUT_DIR")
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
manifest <- jsonlite::fromJSON(
  file.path(root, "run_manifest.json"),
  simplifyVector = FALSE
)

stopifnot(identical(manifest$run_status, "completed"))
stopifnot(identical(manifest$mode, "single"))
stopifnot(identical(manifest$pipeline_version, "0.1.0"))
stopifnot(grepl("^[0-9a-f]{40}$", manifest$git_commit))
stopifnot(grepl("^R version 4[.]", manifest$interpreter$r))
stopifnot(grepl("Python 3[.]12", manifest$interpreter$python))
stopifnot(identical(manifest$cli$refresh_taxonomy, FALSE))
stopifnot(identical(as.integer(manifest$taxonomy$unresolved_count), 46L))
stopifnot(identical(as.integer(manifest$taxonomy$conflicts_count), 26L))

source_counts <- unlist(manifest$taxonomy$resolution_source_counts)
stopifnot(all(c("source_cache", "assignment", "ncbi_refresh", "unresolved") %in%
                names(source_counts)))
stopifnot(as.integer(source_counts[["unresolved"]]) == 46L)

required <- c(
  "resolved_config.yml",
  "session_info.txt",
  "run_manifest.json",
  "01_QC/classification_reconciliation.tsv",
  "01_QC/read_length_by_status.tsv",
  "02_Alpha_Diversity/alpha_diversity.tsv",
  "02_Alpha_Diversity/rarefaction_curve.tsv",
  "02_Alpha_Diversity/rarefaction_resamples.tsv",
  "04_Taxa_Composition/classification_fraction.tsv",
  "07_Kreport/AmbarAyunda_minimap2_16S.kreport",
  "07_Kreport/taxonomy_resolution.tsv",
  "07_Kreport/taxonomy_resolution_sources.tsv",
  "07_Kreport/unresolved_taxids.tsv",
  "07_Kreport/taxonomy_conflicts.tsv",
  "07_Kreport/taxonomy_provenance.json"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Missing release outputs: ", paste(missing, collapse = ", "))

reconciliation <- read.delim(
  file.path(root, "01_QC/classification_reconciliation.tsv"),
  check.names = FALSE
)
stopifnot(sum(reconciliation$TotalReads) == 114056L)
stopifnot(sum(reconciliation$AbundanceClassified) == 80556L)
stopifnot(sum(reconciliation$AbundanceUnclassified) == 33500L)

conflicts <- read.delim(
  file.path(root, "07_Kreport/taxonomy_conflicts.tsv"),
  check.names = FALSE
)
stopifnot(nrow(conflicts) == 26L)

resolution <- read.delim(
  file.path(root, "07_Kreport/taxonomy_resolution.tsv"),
  check.names = FALSE
)
stopifnot(identical(
  sort(unique(resolution$ResolutionSource)),
  sort(c("assignment", "source_cache", "unresolved"))
))

cat("Release integration verification passed.\n")
