#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: verify_release_run.R OUTPUT_DIR")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not locate verify_release_run.R.")
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "analysis", "utils", "config.R"))
source(file.path(repo_root, "analysis", "utils", "manifest.R"))
expected_pipeline_version <- trimws(readLines(
  file.path(repo_root, "VERSION"), n = 1L, warn = FALSE
))
lockfile_path <- file.path(repo_root, "renv.lock")
if (!file.exists(lockfile_path)) stop("Expected committed renv.lock.")
expected_lockfile_sha256 <- digest::digest(file = lockfile_path, algo = "sha256")
expected_git_commit <- trimws(system2(
  "git",
  c("-c", sprintf("safe.directory=%s", repo_root), "-C", repo_root, "rev-parse", "HEAD"),
  stdout = TRUE,
  stderr = TRUE
))
if (length(expected_git_commit) != 1L || !grepl("^[0-9a-f]{40}$", expected_git_commit)) {
  stop("Could not resolve the release checkout's exact Git commit.")
}
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
manifest <- jsonlite::fromJSON(
  file.path(root, "run_manifest.json"),
  simplifyVector = FALSE
)

require_json_array <- function(value, path) {
  if (!is.list(value) || !is.null(names(value))) {
    stop(sprintf("Expected JSON array at '%s'.", path))
  }
  value
}

array_strings <- function(value, path) {
  value <- require_json_array(value, path)
  values <- vapply(value, function(item) {
    if (!is.character(item) || length(item) != 1L || is.na(item)) {
      stop(sprintf("Expected string array item at '%s'.", path))
    }
    item
  }, character(1))
  unname(values)
}

stopifnot(identical(manifest$run_status, "completed"))
stopifnot(identical(manifest$mode, "single"))
stopifnot(identical(manifest$schema_version, 2L))
stopifnot(identical(manifest$schema_revision, 2L))
stopifnot(identical(manifest$config_schema_version, 1L))
stopifnot(identical(manifest$pipeline_version, expected_pipeline_version))
validate_manifest_v2(manifest, physical_root = root)
stopifnot(isTRUE(manifest$environment$locked))
stopifnot(identical(manifest$environment$lock_status, "synchronized"))
stopifnot(identical(manifest$environment$lockfile, "renv.lock"))
stopifnot(identical(manifest$environment$lockfile_sha256, expected_lockfile_sha256))
stopifnot(identical(manifest$git_commit, expected_git_commit))
stopifnot(identical(manifest$git_dirty, FALSE))
stopifnot(identical(manifest$cli$allow_dirty, FALSE))
stopifnot(identical(manifest$cli$allow_unlocked, FALSE))
stopifnot(grepl("^R version 4[.]", manifest$interpreter$r))
stopifnot(grepl("Python 3[.]12", manifest$interpreter$python))
stopifnot(identical(manifest$cli$refresh_taxonomy, FALSE))
stopifnot(identical(manifest$upstream_contract$classifier, "minimap2"))
stopifnot(identical(manifest$upstream_contract$database_set, "ncbi_16s_18s_28s_ITS"))
stopifnot(identical(manifest$upstream_contract$taxonomic_rank, "S"))
stopifnot(is.null(manifest$upstream_contract$workflow_version))
stopifnot(is.null(manifest$upstream_contract$workflow_revision))

# Expected module set completeness and status verification.  The registry is
# fixed so a release manifest cannot silently omit a maintained module.
ALL_MODULES <- c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport", "faprotax")
expected_modules <- array_strings(manifest$cli$modules, "cli.modules")
if (length(expected_modules) == 0L) expected_modules <- ALL_MODULES
invisible(require_json_array(manifest$samples, "samples"))
invisible(require_json_array(manifest$command, "command"))
invisible(require_json_array(manifest$warnings, "warnings"))
invisible(require_json_array(manifest$package_versions, "package_versions"))

# Manifest v2 revision 2 artifacts and preserved unowned outputs verification
owned_files <- array_strings(manifest$owned_outputs, "owned_outputs")
preserved_files <- array_strings(manifest$preserved_unowned_outputs, "preserved_unowned_outputs")
artifacts <- require_json_array(manifest$artifacts, "artifacts")
artifact_paths <- vapply(artifacts, function(a) as.character(a$relative_path), character(1))
stopifnot(setequal(tolower(artifact_paths), tolower(setdiff(owned_files, "run_manifest.json"))))
for (art in artifacts) {
  stopifnot(nzchar(art$relative_path))
  stopifnot(is.numeric(art$size_bytes), art$size_bytes >= 0)
  stopifnot(grepl("^[0-9a-f]{64}$", art$sha256))
  if (art$relative_path %in% c("resolved_config.yml", "session_info.txt")) {
    stopifnot(is.null(art$producer_module))
  } else {
    stopifnot(art$producer_module %in% ALL_MODULES)
  }
}

# Physical file census check: physical files must equal (owned_outputs \ run_manifest.json) ∪ preserved_unowned_outputs
stopifnot(file.exists(file.path(root, "run_manifest.json")))
physical_files <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
physical_files <- physical_files[!dir.exists(file.path(root, physical_files))]
physical_without_manifest <- setdiff(gsub("\\\\", "/", physical_files), "run_manifest.json")
expected_physical <- sort(unique(c(setdiff(owned_files, "run_manifest.json"), preserved_files)))
stopifnot(setequal(tolower(physical_without_manifest), tolower(expected_physical)))

# Non-self-referential check: every maintained module has an explicit record.
stopifnot(identical(sort(names(manifest$modules)), sort(ALL_MODULES)))
stopifnot(setequal(expected_modules, intersect(expected_modules, ALL_MODULES)))

# Module statuses, output file existence, and per-module warnings check
for (mod in names(manifest$modules)) {
  mod_rec <- manifest$modules[[mod]]
  out_paths <- array_strings(mod_rec$outputs, paste0("modules.", mod, ".outputs"))
  warnings <- array_strings(mod_rec$warnings, paste0("modules.", mod, ".warnings"))
  if (mod %in% expected_modules) {
    stopifnot(mod_rec$status %in% c("completed", "skipped", "failed"))
  } else {
    stopifnot(identical(mod_rec$status, "not_run"))
    stopifnot(length(out_paths) == 0L, length(warnings) == 0L)
    stopifnot(is.null(mod_rec$error), !is.null(mod_rec$reason), nzchar(mod_rec$reason))
    stopifnot(is.null(mod_rec$start_time), is.null(mod_rec$end_time), is.null(mod_rec$duration_seconds))
  }
  if (length(out_paths) > 0L) {
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

faprotax_required <- c(
  "08_FAPROTAX/faprotax_function_abundance.tsv",
  "08_FAPROTAX/faprotax_mapping_coverage.tsv",
  "08_FAPROTAX/faprotax_taxon_function_assignments.tsv",
  "08_FAPROTAX/faprotax_top_functions.png",
  "08_FAPROTAX/faprotax_provenance.json"
)

if ("faprotax" %in% expected_modules) {
  stopifnot(manifest$modules$faprotax$status %in% c("completed", "skipped"))
  if (identical(manifest$modules$faprotax$status, "completed")) {
    faprotax_missing <- faprotax_required[!file.exists(file.path(root, faprotax_required))]
    if (length(faprotax_missing)) {
      stop("Missing FAPROTAX outputs: ", paste(faprotax_missing, collapse = ", "))
    }
  } else {
    skip_file <- file.path(root, "08_FAPROTAX/faprotax_skipped.tsv")
    if (!file.exists(skip_file)) {
      stop("FAPROTAX was skipped but '08_FAPROTAX/faprotax_skipped.tsv' is missing.")
    }
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

# 3a. Multi-rank composition artifacts and sidecar schemas
resolved_config <- yaml::read_yaml(file.path(root, "resolved_config.yml"))
composition_cfg <- resolved_config$composition
stacked_sidecar_columns <- c(
  "Rank", "TaxonPath", "DisplayTaxon", "SampleID", "Group",
  "RelativeAbundance", "MeanRelativeAbundance", "IsOther",
  "ValidDenominator", "StackOrder", "SampleOrder"
)
heatmap_sidecar_columns <- c(
  "Rank", "TaxonPath", "DisplayTaxon", "SampleID", "Group",
  "RelativeAbundance", "Transform", "PseudoCount", "TransformedValue",
  "IsOther", "RowOrder", "ColumnOrder", "ValidDenominator"
)
assert_composition_sidecar <- function(path, expected_columns, rank, kind) {
  if (!file.exists(path)) stop(sprintf("Missing %s sidecar: %s", kind, path))
  table <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(table), expected_columns)) {
    stop(sprintf("Unexpected %s sidecar schema for rank '%s'.", kind, rank))
  }
  stopifnot(all(table$Rank == rank), nrow(table) > 0L)
  stopifnot(!any(grepl("(^|;)Unclassified($|;)", table$TaxonPath, ignore.case = TRUE)))
  table
}
assert_composition_pair <- function(rank) {
  tsv <- file.path(root, "04_Taxa_Composition", sprintf("04_%s_stacked.tsv", rank))
  png <- file.path(root, "04_Taxa_Composition", sprintf("04_%s_stacked.png", rank))
  skipped <- sub("[.]tsv$", "_skipped.tsv", tsv)
  if (file.exists(skipped)) {
    stopifnot(!file.exists(tsv), !file.exists(png))
    return(invisible(NULL))
  }
  stopifnot(file.exists(tsv), file.exists(png))
  table <- assert_composition_sidecar(tsv, stacked_sidecar_columns, rank, "stacked")
  for (sample_id in unique(table$SampleID)) {
    sample <- table[table$SampleID == sample_id, , drop = FALSE]
    valid <- sample$ValidDenominator
    if (any(valid)) stopifnot(abs(sum(sample$RelativeAbundance[valid]) - 1) < 1e-8)
  }
  invisible(NULL)
}
assert_heatmap_pair <- function(rank) {
  tsv <- file.path(root, "04_Taxa_Composition", sprintf("04_heatmap_%s.tsv", rank))
  png <- file.path(root, "04_Taxa_Composition", sprintf("04_heatmap_%s.png", rank))
  skipped <- sub("[.]tsv$", "_skipped.tsv", tsv)
  if (file.exists(skipped)) {
    stopifnot(!file.exists(tsv), !file.exists(png))
    return(invisible(NULL))
  }
  stopifnot(file.exists(tsv), file.exists(png))
  table <- assert_composition_sidecar(tsv, heatmap_sidecar_columns, rank, "heatmap")
  stopifnot(all(table$Transform == composition_cfg$heatmap_transform))
  invisible(NULL)
}
for (rank in as.character(composition_cfg$stacked_bar_ranks)) assert_composition_pair(rank)
for (rank in as.character(composition_cfg$heatmap_ranks)) assert_heatmap_pair(rank)

if (isTRUE(manifest$cli$krona)) {
  krona_dir <- file.path(root, "07_Kreport", "krona")
  krona_provenance_file <- file.path(krona_dir, "krona_provenance.json")
  if (!file.exists(krona_provenance_file)) {
    stop("Krona was enabled but '07_Kreport/krona/krona_provenance.json' is missing.")
  }
  krona_provenance <- jsonlite::fromJSON(krona_provenance_file, simplifyVector = FALSE)
  krona_samples <- array_strings(manifest$samples, "samples")
  sanitize_release_filename <- function(sample_id) {
    gsub("[^A-Za-z0-9_.-]", "_", sample_id)
  }
  expected_tsv <- sort(paste0(vapply(krona_samples, sanitize_release_filename, character(1)), ".krona.tsv"))
  actual_tsv <- sort(list.files(krona_dir, pattern = "[.]krona[.]tsv$", full.names = FALSE))
  if (!identical(actual_tsv, expected_tsv)) {
    stop(sprintf(
      "Krona TSV artifact mismatch; expected [%s], found [%s].",
      paste(expected_tsv, collapse = ", "), paste(actual_tsv, collapse = ", ")
    ))
  }

  html_status <- if (is.null(krona_provenance$html_status)) {
    ""
  } else {
    as.character(krona_provenance$html_status)
  }
  expected_html <- if (identical(html_status, "rendered")) {
    sort(paste0(vapply(krona_samples, sanitize_release_filename, character(1)), ".krona.html"))
  } else {
    character(0)
  }
  actual_html <- sort(list.files(krona_dir, pattern = "[.]krona[.]html$", full.names = FALSE))
  if (!identical(actual_html, expected_html)) {
    stop(sprintf(
      "Krona HTML artifact mismatch for status '%s'; expected [%s], found [%s].",
      html_status, paste(expected_html, collapse = ", "), paste(actual_html, collapse = ", ")
    ))
  }
  stopifnot(html_status %in% c("not_requested", "rendered"))
  stopifnot(!is.null(krona_provenance$renderer_policy))
  stopifnot(!is.null(krona_provenance$vendor_sha256_manifest))
  invisible(require_json_array(krona_provenance$vendor_sha256_manifest,
                               "krona_provenance.vendor_sha256_manifest"))
  invisible(require_json_array(krona_provenance$samples, "krona_provenance.samples"))
  stopifnot(length(krona_provenance$samples) == length(krona_samples))

  krona_records <- krona_provenance$samples
  krona_record_ids <- vapply(krona_records, function(record) as.character(record$sample_id), character(1))
  stopifnot(setequal(krona_record_ids, krona_samples))

  read_krona_totals <- function(path) {
    lines <- readLines(path, warn = FALSE)
    stopifnot(length(lines) > 0L, all(nzchar(lines)))
    fields <- strsplit(lines, "\t", fixed = TRUE)
    stopifnot(all(vapply(fields, length, integer(1)) >= 2L))
    magnitudes <- suppressWarnings(as.numeric(vapply(fields, function(field) field[[1]], character(1))))
    stopifnot(all(is.finite(magnitudes)), all(magnitudes >= 0))
    stopifnot(all(abs(magnitudes - round(magnitudes)) <= sqrt(.Machine$double.eps)))
    is_unclassified <- vapply(fields, function(field) identical(field[[2]], "Unclassified"), logical(1))
    stopifnot(sum(is_unclassified) <= 1L)
    list(
      total = sum(magnitudes),
      classified = sum(magnitudes[!is_unclassified]),
      unclassified = sum(magnitudes[is_unclassified])
    )
  }

  for (record in krona_records) {
    sample_id <- as.character(record$sample_id)
    accounting_row <- accounting[accounting$SampleID == sample_id, , drop = FALSE]
    stopifnot(nrow(accounting_row) == 1L)
    tsv_path <- as.character(record$tsv_path)
    stopifnot(file.exists(tsv_path))
    stopifnot(identical(basename(tsv_path), paste0(sanitize_release_filename(sample_id), ".krona.tsv")))
    totals <- read_krona_totals(tsv_path)
    stopifnot(totals$total == accounting_row$TotalReads)
    stopifnot(totals$classified == accounting_row$ClassifiedReads)
    stopifnot(totals$unclassified == accounting_row$UnclassifiedReads)
    stopifnot(as.numeric(record$total_reads) == accounting_row$TotalReads)
    stopifnot(as.numeric(record$classified_reads) == accounting_row$ClassifiedReads)
    stopifnot(as.numeric(record$unclassified_reads) == accounting_row$UnclassifiedReads)
    stopifnot(as.numeric(record$emitted_magnitude_sum) == totals$total)

    if (identical(html_status, "rendered")) {
      html_path <- as.character(record$html_path)
      stopifnot(file.exists(html_path), file.info(html_path)$size > 0)
      stopifnot(identical(basename(html_path), paste0(sanitize_release_filename(sample_id), ".krona.html")))
    } else {
      stopifnot(is.null(record$html_path))
    }
  }
}

# 4. Dataset-specific verification
if (identical(manifest$project_name, "AmbarAyunda_16S_Amplicon")) {
  stopifnot(identical(manifest$upstream_contract$wf_agent, "epi2melabs/5.2.5"))
  unexpected_warnings <- unlist(manifest$warnings, use.names = FALSE)
  conflict_count <- as.integer(manifest$taxonomy$conflicts_count %||% 0L)
  expected_conflict_warning <- sprintf(
    "%d lineage-to-TaxID conflict(s) used the documented modal-count/minimum-TaxID tie-break; review taxonomy_conflicts.tsv.",
    conflict_count
  )
  if (conflict_count > 0L) {
    if (sum(unexpected_warnings == expected_conflict_warning) != 1L) {
      stop("Expected exactly one canonical taxonomy-conflict run warning.")
    }
    unexpected_warnings <- unexpected_warnings[unexpected_warnings != expected_conflict_warning]
  }
  if (length(unexpected_warnings) > 0L) {
    stop("Unexpected run warning(s): ", paste(unexpected_warnings, collapse = "; "))
  }
  for (mod in names(manifest$modules)) {
    module_warnings <- unlist(manifest$modules[[mod]]$warnings, use.names = FALSE)
    if (identical(mod, "kreport") && conflict_count > 0L) {
      if (sum(module_warnings == expected_conflict_warning) != 1L) {
        stop("Expected exactly one canonical taxonomy-conflict kreport warning.")
      }
      module_warnings <- module_warnings[module_warnings != expected_conflict_warning]
    }
    if (length(module_warnings) > 0L) {
      stop(sprintf("Unexpected module warning(s) in '%s': %s", mod, paste(module_warnings, collapse = "; ")))
    }
  }
  stopifnot(identical(as.integer(manifest$taxonomy$unresolved_count), 46L))
  stopifnot(identical(as.integer(manifest$taxonomy$conflicts_count), 26L))

  source_counts <- unlist(manifest$taxonomy$resolution_source_counts)
  expected_source_count_names <- c(
    "source_cache", "assignment", "assignment_conflict", "ncbi_refresh", "unresolved"
  )
  stopifnot(all(expected_source_count_names %in% names(source_counts)))
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
    if (identical(manifest$modules$faprotax$status, "completed")) {
      required <- c(required, faprotax_required)
    } else if (identical(manifest$modules$faprotax$status, "skipped")) {
      required <- c(required, "08_FAPROTAX/faprotax_skipped.tsv")
    }
  }
  missing <- required[!file.exists(file.path(root, required))]
  if (length(missing)) stop("Missing release outputs: ", paste(missing, collapse = ", "))

  stopifnot(sum(reconciliation$TotalReads) == 114056L)
  stopifnot(sum(reconciliation$AbundanceClassified) == 80556L)
  stopifnot(sum(reconciliation$AbundanceUnclassified) == 33500L)

  conflicts <- read.delim(file.path(root, "07_Kreport/taxonomy_conflicts.tsv"), check.names = FALSE)
  stopifnot(nrow(conflicts) == 26L)

  resolution <- read.delim(file.path(root, "07_Kreport/taxonomy_resolution.tsv"), check.names = FALSE)
  expected_sources <- c("assignment", "source_cache", "unresolved")
  if (nrow(conflicts) > 0L) {
    expected_sources <- c(expected_sources, "assignment_conflict")
  }
  stopifnot(identical(
    sort(unique(resolution$ResolutionSource)),
    sort(expected_sources)
  ))
  conflict_resolution <- resolution[resolution$ResolutionSource == "assignment_conflict", , drop = FALSE]
  to_assignment_lineage <- function(path) {
    ranks <- strsplit(path, ";", fixed = TRUE)[[1]]
    if (length(ranks) == 8L) ranks <- ranks[-2L]
    paste(ranks, collapse = "|")
  }
  used_conflict_lineages <- vapply(
    conflict_resolution$TaxonPath, to_assignment_lineage, character(1)
  )
  winner_index <- match(used_conflict_lineages, conflicts$Lineage)
  stopifnot(nrow(conflict_resolution) == as.integer(source_counts[["assignment_conflict"]]))
  stopifnot(!anyDuplicated(conflict_resolution$TaxonPath))
  stopifnot(!anyNA(winner_index))
  stopifnot(all(conflict_resolution$Status == "Conflicted"))
  stopifnot(identical(
    as.character(conflict_resolution$TaxID),
    as.character(conflicts$WinnerTaxID[winner_index])
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
