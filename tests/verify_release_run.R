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
stopifnot(identical(manifest$pipeline_version, "0.2.0"))
stopifnot(grepl("^[0-9a-f]{40}$", manifest$git_commit))
stopifnot(grepl("^R version 4[.]", manifest$interpreter$r))
stopifnot(grepl("Python 3[.]12", manifest$interpreter$python))
stopifnot(identical(manifest$cli$refresh_taxonomy, FALSE))
stopifnot(identical(manifest$upstream_contract$classifier, "minimap2"))
stopifnot(identical(manifest$upstream_contract$database_set, "ncbi_16s_18s_28s_ITS"))
stopifnot(identical(manifest$upstream_contract$taxonomic_rank, "S"))
stopifnot(identical(manifest$upstream_contract$wf_agent, "epi2melabs/5.2.5"))
stopifnot(is.null(manifest$upstream_contract$workflow_version))
stopifnot(is.null(manifest$upstream_contract$workflow_revision))
stopifnot(length(manifest$warnings) == 0L)
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
  "01_QC/00_read_accounting.tsv",
  "01_QC/00_read_investigation.tsv",
  "01_QC/AmbarAyunda_minimap2_16S/00a_read_accounting_donut.png",
  "01_QC/classification_reconciliation.tsv",
  "01_QC/read_length_by_status.tsv",
  "02_Alpha_Diversity/alpha_diversity.tsv",
  "02_Alpha_Diversity/02_richness_overview.tsv",
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

accounting <- read.delim(file.path(root, "01_QC/00_read_accounting.tsv"), check.names = FALSE)
stopifnot(nrow(accounting) == 1L)
stopifnot(identical(accounting$SampleID, "AmbarAyunda_minimap2_16S"))
stopifnot(accounting$AbundanceTotal == 114056L)
stopifnot(accounting$RawC == 89809L, accounting$RawU == 24247L)
stopifnot(accounting$C_TaxID0 == 9253L, accounting$TaxID_GT0 == 80556L)
stopifnot(abs(accounting$EffectiveClassifiedPct - 70.6284632110542) < 1e-10)
stopifnot(abs(accounting$C0ShareOfEffectiveUnclassifiedPct - 27.62089552239) < 1e-10)

investigation <- read.delim(file.path(root, "01_QC/00_read_investigation.tsv"), check.names = FALSE)
stopifnot(investigation$MedianClassifiedLength == 1507)
stopifnot(investigation$MedianC0Length == 1493)
stopifnot(investigation$MedianRawULength == 1494)
stopifnot(identical(investigation$BamstatsAvailable, FALSE))
stopifnot(is.na(investigation$BamstatsC0Matched))
stopifnot(is.na(investigation$IdentityOnlyFailed))
stopifnot(is.na(investigation$RefCoverageOnlyFailed))
stopifnot(is.na(investigation$BothFailed))

richness <- read.delim(file.path(root, "02_Alpha_Diversity/02_richness_overview.tsv"),
                       check.names = FALSE)
stopifnot(richness$ClassifiedReads == 80556L)
stopifnot(richness$PositiveTaxa == 1836L, richness$SingletonTaxa == 735L)
stopifnot(richness$TaxaLeq10 == 1399L, richness$ReadsInTaxaLeq10 == 3456L)
stopifnot(abs(richness$ReadsInTaxaLeq10Pct - 4.2901832265753) < 1e-10)

composition <- read.delim(
  file.path(root, "04_Taxa_Composition/classification_fraction.tsv"), check.names = FALSE
)
stopifnot(!anyNA(composition))
stopifnot(composition$TotalReads == 114056L)
stopifnot(composition$ClassifiedReads == 80556L)
stopifnot(composition$UnclassifiedReads == 33500L)
stopifnot(composition$ClassifiedReads + composition$UnclassifiedReads == composition$TotalReads)

cat("Release integration verification passed.\n")
