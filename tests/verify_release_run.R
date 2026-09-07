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
stopifnot(is.null(manifest$upstream_contract$workflow_version))
stopifnot(is.null(manifest$upstream_contract$workflow_revision))

# Expected module set completeness and status verification
expected_modules <- if (!is.null(manifest$cli$modules) && length(manifest$cli$modules) > 0L) {
  unlist(manifest$cli$modules)
} else {
  c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport")
}
stopifnot(identical(sort(names(manifest$modules)), sort(expected_modules)))

# Module statuses, output file existence, and per-module warnings check
for (mod in names(manifest$modules)) {
  mod_rec <- manifest$modules[[mod]]
  stopifnot(mod_rec$status %in% c("completed", "skipped"))
  if (length(mod_rec$outputs) > 0L) {
    out_paths <- unlist(mod_rec$outputs, use.names = FALSE)
    stopifnot(all(file.exists(out_paths)))
  }
}

# Single-sample mode status assertions
if (identical(manifest$mode, "single")) {
  for (m in intersect(c("qc", "alpha", "composition", "kreport"), expected_modules)) {
    stopifnot(identical(manifest$modules[[m]]$status, "completed"))
  }
  for (m in intersect(c("beta", "ordination", "shared"), expected_modules)) {
    stopifnot(identical(manifest$modules[[m]]$status, "skipped"))
  }
}

if ("faprotax" %in% expected_modules) {
  stopifnot(identical(manifest$modules$faprotax$status, "completed"))
  faprotax_required <- c(
    "08_FAPROTAX/faprotax_function_abundance.tsv",
    "08_FAPROTAX/faprotax_mapping_coverage.tsv",
    "08_FAPROTAX/faprotax_taxon_function_assignments.tsv",
    "08_FAPROTAX/faprotax_top_functions.png",
    "08_FAPROTAX/faprotax_provenance.json"
  )
  faprotax_missing <- faprotax_required[!file.exists(file.path(root, faprotax_required))]
  if (length(faprotax_missing)) {
    stop("Missing FAPROTAX outputs: ", paste(faprotax_missing, collapse = ", "))
  }
}

# 1. Reconciliation table verification
reconciliation <- read.delim(
  file.path(root, "01_QC/classification_reconciliation.tsv"),
  check.names = FALSE
)
stopifnot("ReconciliationPass" %in% names(reconciliation))
stopifnot(all(reconciliation$ReconciliationPass == TRUE))

# 2. Accounting & Investigation verification
accounting <- read.delim(file.path(root, "01_QC/00_read_accounting.tsv"), check.names = FALSE)
stopifnot(nrow(accounting) == 1L)

investigation <- read.delim(file.path(root, "01_QC/00_read_investigation.tsv"), check.names = FALSE)
stopifnot(nrow(investigation) == 1L)

if (isTRUE(investigation$BamstatsAvailable[1]) && isTRUE(investigation$AssignmentAvailable[1])) {
  stopifnot(investigation$BamstatsC0Matched == accounting$C_TaxID0)
  stopifnot(investigation$IdentityOnlyFailed + investigation$RefCoverageOnlyFailed +
            investigation$BothFailed == investigation$BamstatsC0Matched)
  stopifnot(investigation$IdentityOnlyFailed > 0L)
  stopifnot(investigation$RefCoverageOnlyFailed > 0L)
  stopifnot(investigation$BothFailed > 0L)
} else {
  stopifnot(is.na(investigation$BamstatsC0Matched))
  stopifnot(is.na(investigation$IdentityOnlyFailed))
  stopifnot(is.na(investigation$RefCoverageOnlyFailed))
  stopifnot(is.na(investigation$BothFailed))
}

# 3. Composition table verification
composition <- read.delim(
  file.path(root, "04_Taxa_Composition/classification_fraction.tsv"), check.names = FALSE
)
stopifnot(!anyNA(composition))
stopifnot(composition$ClassifiedReads + composition$UnclassifiedReads == composition$TotalReads)

# 4. Dataset-specific verification
if (identical(manifest$project_name, "AmbarAyunda_16S_Amplicon")) {
  stopifnot(identical(manifest$upstream_contract$wf_agent, "epi2melabs/5.2.5"))
  stopifnot(length(manifest$warnings) == 0L)
  for (mod in names(manifest$modules)) {
    stopifnot(length(manifest$modules[[mod]]$warnings) == 0L)
  }
  stopifnot(identical(as.integer(manifest$taxonomy$unresolved_count), 46L))
  stopifnot(identical(as.integer(manifest$taxonomy$conflicts_count), 26L))

  source_counts <- unlist(manifest$taxonomy$resolution_source_counts)
  stopifnot(all(c("source_cache", "assignment", "ncbi_refresh", "unresolved") %in%
                  names(source_counts)))
  stopifnot(as.integer(source_counts[["unresolved"]]) == 46L)

  required <- c(
    "resolved_config.yml", "session_info.txt", "run_manifest.json",
    "01_QC/00_read_accounting.tsv", "01_QC/00_read_investigation.tsv",
    "01_QC/AmbarAyunda_minimap2_16S/00a_read_accounting_donut.png",
    "01_QC/classification_reconciliation.tsv", "01_QC/read_length_by_status.tsv",
    "02_Alpha_Diversity/alpha_diversity.tsv", "02_Alpha_Diversity/02_richness_overview.tsv",
    "02_Alpha_Diversity/rarefaction_curve.tsv", "02_Alpha_Diversity/rarefaction_resamples.tsv",
    "04_Taxa_Composition/classification_fraction.tsv",
    "07_Kreport/AmbarAyunda_minimap2_16S.kreport", "07_Kreport/taxonomy_resolution.tsv",
    "07_Kreport/taxonomy_resolution_sources.tsv", "07_Kreport/unresolved_taxids.tsv",
    "07_Kreport/taxonomy_conflicts.tsv", "07_Kreport/taxonomy_provenance.json"
  )
  if ("faprotax" %in% expected_modules) {
    required <- c(required, faprotax_required)
  }
  missing <- required[!file.exists(file.path(root, required))]
  if (length(missing)) stop("Missing release outputs: ", paste(missing, collapse = ", "))

  stopifnot(sum(reconciliation$TotalReads) == 114056L)
  stopifnot(sum(reconciliation$AbundanceClassified) == 80556L)
  stopifnot(sum(reconciliation$AbundanceUnclassified) == 33500L)

  conflicts <- read.delim(file.path(root, "07_Kreport/taxonomy_conflicts.tsv"), check.names = FALSE)
  stopifnot(nrow(conflicts) == 26L)

  resolution <- read.delim(file.path(root, "07_Kreport/taxonomy_resolution.tsv"), check.names = FALSE)
  stopifnot(identical(
    sort(unique(resolution$ResolutionSource)),
    sort(c("assignment", "source_cache", "unresolved"))
  ))

  stopifnot(identical(accounting$SampleID, "AmbarAyunda_minimap2_16S"))
  stopifnot(accounting$AbundanceTotal == 114056L)
  stopifnot(accounting$RawC == 89809L, accounting$RawU == 24247L)
  stopifnot(accounting$C_TaxID0 == 9253L, accounting$TaxID_GT0 == 80556L)
  stopifnot(abs(accounting$EffectiveClassifiedPct - 70.6284632110542) < 1e-10)
  stopifnot(abs(accounting$C0ShareOfEffectiveUnclassifiedPct - 27.62089552239) < 1e-10)

  stopifnot(investigation$MedianClassifiedLength == 1507)
  stopifnot(investigation$MedianC0Length == 1493)
  stopifnot(investigation$MedianRawULength == 1494)
  stopifnot(identical(investigation$BamstatsAvailable, FALSE))

  richness <- read.delim(file.path(root, "02_Alpha_Diversity/02_richness_overview.tsv"), check.names = FALSE)
  stopifnot(richness$ClassifiedReads == 80556L)
  stopifnot(richness$PositiveTaxa == 1836L, richness$SingletonTaxa == 735L)
  stopifnot(richness$TaxaLeq10 == 1399L, richness$ReadsInTaxaLeq10 == 3456L)
  stopifnot(abs(richness$ReadsInTaxaLeq10Pct - 4.2901832265753) < 1e-10)

  stopifnot(composition$TotalReads == 114056L)
  stopifnot(composition$ClassifiedReads == 80556L)
  stopifnot(composition$UnclassifiedReads == 33500L)

} else if (identical(manifest$project_name, "synthetic_minimap2_bamstats")) {
  stopifnot(identical(manifest$upstream_contract$wf_agent, "epi2melabs/wf-16s"))
  stopifnot(identical(as.integer(manifest$taxonomy$unresolved_count), 0L))
  stopifnot(identical(as.integer(manifest$taxonomy$conflicts_count), 0L))

  stopifnot(sum(reconciliation$TotalReads) == 5L)
  stopifnot(sum(reconciliation$AbundanceClassified) == 1L)
  stopifnot(sum(reconciliation$AbundanceUnclassified) == 4L)

  stopifnot(identical(accounting$SampleID, "synthetic_minimap2"))
  stopifnot(accounting$AbundanceTotal == 5L)
  stopifnot(accounting$RawC == 4L, accounting$RawU == 1L)
  stopifnot(accounting$C_TaxID0 == 3L, accounting$TaxID_GT0 == 1L)
  stopifnot(abs(accounting$EffectiveClassifiedPct - 20) < 1e-10)
  stopifnot(abs(accounting$C0ShareOfEffectiveUnclassifiedPct - 75) < 1e-10)

  stopifnot(isTRUE(investigation$BamstatsAvailable[1]))
  stopifnot(investigation$BamstatsC0Matched == 3L)
  stopifnot(investigation$IdentityOnlyFailed == 1L)
  stopifnot(investigation$RefCoverageOnlyFailed == 1L)
  stopifnot(investigation$BothFailed == 1L)

  stopifnot(composition$TotalReads == 5L)
  stopifnot(composition$ClassifiedReads == 1L)
  stopifnot(composition$UnclassifiedReads == 4L)
} else {
  stop(sprintf("Unknown project_name for release verification: '%s'", manifest$project_name))
}

cat(sprintf("Release integration verification passed for '%s'.\n", manifest$project_name))
