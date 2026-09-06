# Patch plan: `ont-wf16S-amplicon-postprocess` v0.2.0

**Handoff target:** Gemini

**Repository:** `WhyAdr/ont-wf16S-amplicon-postprocess`

**Audited base:** `1c0e3ba` on `main`

**State review:** `ont-wf16s-postprocess-v0.1.0-statereview.md`

**Scope:** exact read accounting, optional minimap2 failure diagnostics, a concise
species-richness overview, preservation of established unclassified-read
denominators, the P0 composition fix, and enforceable upstream provenance

## 1. Audit verdict

The original v0.2.0 draft was directionally correct but was not yet executable.
This revision locks the implementation to the real v0.1.0 functions and the
producer-native schemas observed in the tracked Ambar fixture and local wf-16s
v1.6.1 outputs.

Corrections made during this audit:

1. `run_qc()` currently returns immediately without assignments. It must instead
   always write abundance-derived accounting rows and use `NA` only for fields
   that truly require assignments or bamstats.
2. Do not select the first recursive bamstats match. The filename
   `bamstats.readstats.tsv.gz` repeats once per sample. Discover candidates by
   their internal `sample_name`, map them to `SampleID`, and fail on ambiguity.
3. The observed bamstats columns are `name`, `sample_name`, `iden`, and
   `ref_coverage`. Join assignment `read_id` to bamstats `name` one-to-one.
4. Threshold partitions must conserve every `status=C, TaxID=0` read. A C0 read
   that passes both thresholds, or a `TaxID>0` read that fails either threshold,
   is contract drift and must fail loudly when bamstats is configured.
5. `wf.agent` (for example `epi2melabs/5.2.5`) is not the wf-16s workflow
   version or Nextflow revision. Record it as `wf_agent`; leave
   `workflow_version` and `workflow_revision` null in this patch.
6. Classifier-only validation is insufficient. Require `params.json`, reject
   non-minimap2 classifiers, non-NCBI bundled database sets, non-species rank,
   and custom reference/taxonomy overrides before parsing assignments.
7. The read-accounting pie duplicates the existing three-way diagnostic bar.
   Retain both for backward compatibility, but call the new graphic a
   three-segment accounting donut and do not renumber existing `01a`-`01d`
   artifacts.
8. The state-review percentages were re-derived. Percent columns below are
   numeric percentages on a 0-100 scale and carry a `Pct` suffix.
9. The integration verifier cannot retrospectively detect plotting warnings.
   Capture warnings in the runner, attach them to each module result, and assert
   that the supported Ambar run has none.
10. The tracked Ambar fixture excludes `output_AAy/bams/` via `.gitignore`.
    Therefore `config.yml` must keep bamstats discovery disabled, CI must assert
    explicit `NA` partitions, and exact native bamstats values remain opt-in
    local-fixture checks.

No SILVA or Kraken2 adapter is introduced. The five local producer-result
directories remain untracked and are optional manual integration fixtures, not
CI dependencies.

## 2. Locked behavior and output contracts

### 2.1 Supported upstream contract

`input.params_json` is required. v0.2.0 accepts only:

- `classifier == "minimap2"`;
- `database_set` in `ncbi_16s_18s` or `ncbi_16s_18s_28s_ITS`;
- `taxonomic_rank == "S"`;
- no custom `taxonomy`, `reference`, `ref2taxid`, or `database` override; and
- finite thresholds with `min_len < max_len`, identity/coverage in `[0,100]`,
  and a non-negative abundance threshold.

This closes the abundance-only Kraken2 bypass described in state-review P0-2.
The existing exact eight-rank abundance validator remains a second, independent
schema gate.

### 2.2 `01_QC/00_read_accounting.tsv`

One row for every selected sample, even when assignments are absent:

| Column | Unit/source | Nullability |
|---|---|---|
| `SampleID` | abundance sample identifier | never null |
| `AssignmentAvailable` | logical | never null |
| `AbundanceTotal` | abundance table | never null |
| `AbundanceClassified` | abundance table excluding canonical unclassified row | never null |
| `AbundanceUnclassified` | canonical unclassified row | never null |
| `RawC` | assignment `status=C` | `NA` without assignments |
| `RawU` | assignment `status=U` | `NA` without assignments |
| `C_TaxID0` | assignment `status=C && TaxID=0` | `NA` without assignments |
| `TaxID_GT0` | assignment `TaxID>0` | `NA` without assignments |
| `EffectiveClassifiedPct` | `100 * AbundanceClassified / AbundanceTotal` | never null |
| `C0ShareOfEffectiveUnclassifiedPct` | `100 * C_TaxID0 / AbundanceUnclassified` | `NA` without assignments |

Required conservation checks:

```text
AbundanceClassified + AbundanceUnclassified == AbundanceTotal
RawC + RawU == AbundanceTotal                         [when assignments exist]
RawU + C_TaxID0 == AbundanceUnclassified              [when assignments exist]
TaxID_GT0 == AbundanceClassified                      [when assignments exist]
```

For each sample with assignments, write
`01_QC/<sanitized SampleID>/00a_read_accounting_donut.png`. Segments are
`Raw U` (`#bdbdbd`), `C + TaxID 0` (`#e6ab02`), and `TaxID > 0` (`#1b9e77`).
Labels contain count and percent; the centre contains total reads.

### 2.3 `01_QC/00_read_investigation.tsv`

One row per selected sample:

| Column | Meaning |
|---|---|
| `SampleID` | sample identifier |
| `AssignmentAvailable` | whether a configured assignment was parsed |
| `BamstatsAvailable` | whether a unique bamstats file was mapped to the sample |
| `MedianClassifiedLength` | median assignment length for `TaxID>0` |
| `MedianC0Length` | median assignment length for `status=C, TaxID=0` |
| `MedianRawULength` | median assignment length for `status=U` |
| `MinPercentIdentity` | upstream threshold from `params.json` |
| `MinRefCoverage` | upstream threshold from `params.json` |
| `BamstatsC0Matched` | number of C0 reads joined one-to-one |
| `IdentityOnlyFailed` | `iden < identity threshold`, coverage passed |
| `RefCoverageOnlyFailed` | identity passed, `ref_coverage < coverage threshold` |
| `BothFailed` | both thresholds failed |

Medians are `NA` without assignments. Bamstats-derived fields are `NA` unless
both assignments and a unique bamstats file are available. When available,
`BamstatsC0Matched == IdentityOnlyFailed + RefCoverageOnlyFailed + BothFailed ==
C_TaxID0` is mandatory.

`input.wf16s_output_root` is optional. If configured, it must exist. Zero
matching bamstats files is an informational absence; duplicate candidates for a
selected sample are an error. Bamstats paths and SHA-256 hashes enter the run
manifest.

### 2.4 `02_Alpha_Diversity/02_richness_overview.tsv`

One classified-only row per sample:

| Column | Meaning |
|---|---|
| `SampleID` | sample identifier |
| `ClassifiedReads` | classified-only denominator |
| `PositiveTaxa` | taxa with count `>0` |
| `SingletonTaxa` | taxa with count `==1` |
| `SingletonPct` | percent of positive taxa that are singletons |
| `TaxaLeq10` | positive taxa with count `<=10` |
| `ReadsInTaxaLeq10` | reads assigned to those taxa |
| `ReadsInTaxaLeq10Pct` | percent of classified reads in those taxa |

The table is descriptive and threshold-sensitive. It does not identify which
low-count taxa are artifacts.

### 2.5 Preserved unclassified-read handling

| Module | v0.2.0 behavior |
|---|---|
| QC | abundance totals always included; assignment plots use all reads when available |
| Alpha diversity | excluded from all diversity/richness denominators |
| Beta diversity | excluded |
| Composition | excluded from rank tables; included in all-read classification fraction |
| Ordination | excluded |
| Shared taxa | excluded |
| Kraken/Pavian | included as the first `U` line and in total-read arithmetic |

## 3. Exact implementation diff blocks

Apply these blocks in order. They are scoped to the audited base; do not stage
the untracked producer outputs or `ATW_Sesame_Greenhouse-Examples/`.

### 3.1 Version metadata

```diff
diff --git a/VERSION b/VERSION
--- a/VERSION
+++ b/VERSION
@@ -1,1 +1,1 @@
-0.1.0
+0.2.0
diff --git a/CITATION.cff b/CITATION.cff
--- a/CITATION.cff
+++ b/CITATION.cff
@@ -6,2 +6,1 @@
-version: 0.1.0
-date-released: 2026-09-05
+version: 0.2.0
```

Do not add a new `date-released` until the release/tag date is known.

### 3.2 Configuration contract and paths

```diff
diff --git a/analysis/utils/config.R b/analysis/utils/config.R
--- a/analysis/utils/config.R
+++ b/analysis/utils/config.R
@@ -39,3 +39,7 @@
   assert_nonempty_string(cfg$input$abundance_table, "input.abundance_table")
+  assert_nonempty_string(cfg$input$params_json, "input.params_json")
   assert_nonempty_string(cfg$input$tax_column, "input.tax_column")
   assert_nonempty_string(cfg$output$base_dir, "output.base_dir")
+  if (!is.null(cfg$input$wf16s_output_root)) {
+    assert_nonempty_string(cfg$input$wf16s_output_root, "input.wf16s_output_root")
+  }
@@ -153,4 +157,5 @@
       abundance_table = "output_AAy/abundance_table_species.tsv",
       metadata = NULL,
       params_json = "output_AAy/params.json",
+      wf16s_output_root = NULL,
       tax_column = "tax",
@@ -271,3 +276,4 @@
   cfg$input$abundance_table <- resolve_path(cfg$input$abundance_table, config_dir)
   cfg$input$metadata <- resolve_path(cfg$input$metadata, config_dir)
   cfg$input$params_json <- resolve_path(cfg$input$params_json, config_dir)
+  cfg$input$wf16s_output_root <- resolve_path(cfg$input$wf16s_output_root, config_dir)
```

```diff
diff --git a/config.example.yml b/config.example.yml
--- a/config.example.yml
+++ b/config.example.yml
@@ -24,2 +24,5 @@
-  # Upstream wf-16s params.json for reading min_len, max_len, abundance_threshold
+  # Required producer contract and threshold provenance
   params_json: "output_AAy/params.json"
+
+  # Optional root for per-sample bamstats.readstats.tsv.gz discovery (null in CI)
+  wf16s_output_root: null
diff --git a/config_helga.yml b/config_helga.yml
--- a/config_helga.yml
+++ b/config_helga.yml
@@ -8,3 +8,4 @@
   metadata: null
   params_json: "wf-16s_Helga16SrRNA/output/params.json"
+  wf16s_output_root: "wf-16s_Helga16SrRNA/output"
   tax_column: "tax"
```

### 3.3 Producer-contract and bamstats helpers

Insert the following after `sanitize_filename()` in `analysis/utils/io.R`:

```diff
diff --git a/analysis/utils/io.R b/analysis/utils/io.R
--- a/analysis/utils/io.R
+++ b/analysis/utils/io.R
@@ -19,3 +19,184 @@
 sanitize_filename <- function(s) {
   gsub("[^A-Za-z0-9_.-]", "_", s)
 }
+
+SUPPORTED_NCBI_DATABASE_SETS <- c("ncbi_16s_18s", "ncbi_16s_18s_28s_ITS")
+
+read_upstream_params <- function(path) {
+  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
+    stop("input.params_json is required and must identify an existing wf-16s params.json.",
+         call. = FALSE)
+  }
+  params <- tryCatch(
+    jsonlite::fromJSON(path, simplifyVector = FALSE),
+    error = function(e) stop(sprintf("Could not parse params.json '%s': %s", path, e$message),
+                             call. = FALSE)
+  )
+  required <- c("classifier", "database_set", "taxonomic_rank", "min_len", "max_len",
+                "min_read_qual", "min_percent_identity", "min_ref_coverage",
+                "abundance_threshold")
+  missing <- required[vapply(required, function(x) is.null(params[[x]]), logical(1))]
+  if (length(missing)) {
+    stop(sprintf("params.json is missing required contract field(s): %s",
+                 paste(missing, collapse = ", ")), call. = FALSE)
+  }
+  scalar_string <- function(field) {
+    value <- params[[field]]
+    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
+      stop(sprintf("params.json field '%s' must be one non-empty string.", field), call. = FALSE)
+    }
+    value
+  }
+  classifier <- scalar_string("classifier")
+  database_set <- scalar_string("database_set")
+  taxonomic_rank <- scalar_string("taxonomic_rank")
+  if (!identical(classifier, "minimap2")) {
+    stop(sprintf(
+      "Unsupported wf-16s classifier '%s'. v0.2.0 supports minimap2 only; Kraken2/Bracken requires a classifier-specific denominator model.",
+      classifier
+    ), call. = FALSE)
+  }
+  if (!database_set %in% SUPPORTED_NCBI_DATABASE_SETS) {
+    stop(sprintf(
+      "Unsupported wf-16s database_set '%s'. v0.2.0 supports the bundled NCBI database sets only.",
+      database_set
+    ), call. = FALSE)
+  }
+  if (!identical(taxonomic_rank, "S")) {
+    stop(sprintf("Unsupported wf-16s taxonomic_rank '%s'; expected species rank 'S'.",
+                 taxonomic_rank), call. = FALSE)
+  }
+  override_fields <- c("taxonomy", "reference", "ref2taxid", "database")
+  active_overrides <- override_fields[vapply(override_fields, function(field) {
+    value <- params[[field]]
+    !is.null(value) && length(value) > 0L && !all(is.na(value)) && any(nzchar(as.character(value)))
+  }, logical(1))]
+  if (length(active_overrides)) {
+    stop(sprintf("Custom wf-16s reference/taxonomy overrides are unsupported: %s",
+                 paste(active_overrides, collapse = ", ")), call. = FALSE)
+  }
+  numeric_fields <- c("min_len", "max_len", "min_read_qual", "min_percent_identity",
+                      "min_ref_coverage", "abundance_threshold")
+  for (field in numeric_fields) {
+    value <- params[[field]]
+    if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
+      stop(sprintf("params.json field '%s' must be one finite number.", field), call. = FALSE)
+    }
+  }
+  if (params$min_len <= 0 || params$max_len <= params$min_len || params$min_read_qual < 0 ||
+      params$min_percent_identity < 0 || params$min_percent_identity > 100 ||
+      params$min_ref_coverage < 0 || params$min_ref_coverage > 100 ||
+      params$abundance_threshold < 0) {
+    stop("params.json contains an invalid length, quality, identity, coverage, or abundance threshold.",
+         call. = FALSE)
+  }
+  params
+}
+
+extract_upstream_contract <- function(params) {
+  database_meta <- params$database_sets[[params$database_set]]
+  list(
+    workflow_name = "epi2me-labs/wf-16s",
+    workflow_version = NULL,
+    workflow_revision = NULL,
+    wf_agent = params$wf$agent %||% NULL,
+    classifier = params$classifier,
+    database_set = params$database_set,
+    taxonomy_namespace = "NCBI",
+    database_taxonomy_source = database_meta$taxonomy %||% NULL,
+    taxonomic_rank = params$taxonomic_rank,
+    min_len = params$min_len,
+    max_len = params$max_len,
+    min_read_qual = params$min_read_qual,
+    min_percent_identity = params$min_percent_identity,
+    min_ref_coverage = params$min_ref_coverage,
+    abundance_threshold = params$abundance_threshold,
+    output_unclassified = params$output_unclassified %||% NULL,
+    include_read_assignments = params$include_read_assignments %||% NULL
+  )
+}
+
+discover_bamstats <- function(root, sample_ids) {
+  mapped <- stats::setNames(rep(NA_character_, length(sample_ids)), sample_ids)
+  if (is.null(root)) return(mapped)
+  if (!dir.exists(root)) {
+    stop(sprintf("Configured input.wf16s_output_root does not exist: '%s'", root), call. = FALSE)
+  }
+  candidates <- sort(list.files(
+    root, pattern = "^bamstats[.]readstats[.]tsv[.]gz$",
+    recursive = TRUE, full.names = TRUE
+  ))
+  if (!length(candidates)) {
+    message(sprintf("[INFO] No bamstats.readstats.tsv.gz found under '%s'.", root))
+    return(mapped)
+  }
+  for (path in candidates) {
+    probe <- read.delim(gzfile(path), nrows = 1L, check.names = FALSE,
+                        stringsAsFactors = FALSE)
+    required <- c("name", "sample_name", "iden", "ref_coverage")
+    if (!all(required %in% names(probe)) || nrow(probe) != 1L) {
+      stop(sprintf("Invalid bamstats schema or empty file: '%s'", path), call. = FALSE)
+    }
+    sample_id <- as.character(probe$sample_name[[1]])
+    if (!sample_id %in% sample_ids) next
+    if (!is.na(mapped[[sample_id]])) {
+      stop(sprintf("Multiple bamstats files discovered for sample '%s'.", sample_id), call. = FALSE)
+    }
+    mapped[[sample_id]] <- normalizePath(path, winslash = "/", mustWork = TRUE)
+  }
+  missing <- names(mapped)[is.na(mapped)]
+  if (length(missing)) {
+    message(sprintf("[INFO] No bamstats file mapped for sample(s): %s",
+                    paste(missing, collapse = ", ")))
+  }
+  mapped
+}
+
+partition_minimap2_failures <- function(reads, bamstats_path, params, sample_id) {
+  stats <- read.delim(gzfile(bamstats_path), check.names = FALSE,
+                      stringsAsFactors = FALSE)
+  required <- c("name", "sample_name", "iden", "ref_coverage")
+  if (!all(required %in% names(stats))) {
+    stop(sprintf("Bamstats for '%s' lacks required columns: %s", sample_id,
+                 paste(setdiff(required, names(stats)), collapse = ", ")), call. = FALSE)
+  }
+  if (anyNA(stats$name) || any(!nzchar(stats$name)) || anyDuplicated(stats$name)) {
+    stop(sprintf("Bamstats read names for '%s' must be non-empty and unique.", sample_id),
+         call. = FALSE)
+  }
+  if (any(as.character(stats$sample_name) != sample_id)) {
+    stop(sprintf("Bamstats sample_name does not consistently equal '%s'.", sample_id),
+         call. = FALSE)
+  }
+  stats$iden <- suppressWarnings(as.numeric(stats$iden))
+  stats$ref_coverage <- suppressWarnings(as.numeric(stats$ref_coverage))
+  if (any(!is.finite(stats$iden)) || any(!is.finite(stats$ref_coverage))) {
+    stop(sprintf("Bamstats identity/coverage values for '%s' must be finite numbers.", sample_id),
+         call. = FALSE)
+  }
+  c_reads <- reads[reads$status == "C", , drop = FALSE]
+  matched <- match(c_reads$read_id, stats$name)
+  if (anyNA(matched)) {
+    stop(sprintf("Bamstats is missing %d status-C read(s) for '%s'.",
+                 sum(is.na(matched)), sample_id), call. = FALSE)
+  }
+  aligned <- stats[matched, , drop = FALSE]
+  identity_failed <- aligned$iden < params$min_percent_identity
+  coverage_failed <- aligned$ref_coverage < params$min_ref_coverage
+  positive <- c_reads$taxid > 0
+  c0 <- !positive
+  if (any(identity_failed[positive] | coverage_failed[positive])) {
+    stop(sprintf("At least one TaxID>0 read for '%s' fails the recorded thresholds.", sample_id),
+         call. = FALSE)
+  }
+  if (any(!identity_failed[c0] & !coverage_failed[c0])) {
+    stop(sprintf("At least one C+TaxID0 read for '%s' passes both recorded thresholds.", sample_id),
+         call. = FALSE)
+  }
+  list(
+    matched = sum(c0),
+    identity_only = sum(c0 & identity_failed & !coverage_failed),
+    coverage_only = sum(c0 & !identity_failed & coverage_failed),
+    both = sum(c0 & identity_failed & coverage_failed)
+  )
+}
```

Move params parsing ahead of abundance/assignment parsing and attach discovery,
contract data, and hashes in `build_context()`:

```diff
diff --git a/analysis/utils/io.R b/analysis/utils/io.R
--- a/analysis/utils/io.R
+++ b/analysis/utils/io.R
@@ -380,2 +380,6 @@
 build_context <- function(cfg) {
-  # 1. Read abundance table
+  # 1. Fail closed on the producer contract before parsing classifier-specific files.
+  params <- read_upstream_params(cfg$input$params_json)
+  upstream_contract <- extract_upstream_contract(params)
+
+  # 2. Read abundance table
@@ -403,1 +407,1 @@
-  # 2. Mode resolution
+  # 3. Mode resolution
@@ -421,1 +425,1 @@
-  # 3. Read metadata
+  # 4. Read metadata
@@ -428,1 +432,1 @@
-  # 4. Assignments mapping
+  # 5. Assignments mapping
@@ -468,7 +472,4 @@
-  # 5. Read params.json if available
-  params <- NULL
-  if (!is.null(cfg$input$params_json) && file.exists(cfg$input$params_json)) {
-    params <- suppressWarnings(tryCatch(jsonlite::fromJSON(cfg$input$params_json), error = function(e) NULL))
-  }
-
-  # 6. File hashes
+  # 6. Optional per-sample bamstats discovery
+  bamstats <- discover_bamstats(cfg$input$wf16s_output_root, selected_samples)
+
+  # 7. File hashes
@@ -481,5 +482,8 @@
   if (!is.null(assignments_map)) {
     for (s in names(assignments_map)) {
       file_hashes[[paste0("assignment_", s)]] <- compute_file_hash(assignments_map[[s]])
     }
   }
+  for (s in names(bamstats)[!is.na(bamstats)]) {
+    file_hashes[[paste0("bamstats_", s)]] <- compute_file_hash(bamstats[[s]])
+  }
@@ -496,3 +500,5 @@
     assignments = assignments_map,
     assignment_data = assignment_data,
+    bamstats = bamstats,
     params = params,
+    upstream_contract = upstream_contract,
```

### 3.4 Read accounting and investigation outputs

In `analysis/01_qc_diagnostics.R`, remove the early assignment-only skip and add
the two per-sample row collections:

```diff
diff --git a/analysis/01_qc_diagnostics.R b/analysis/01_qc_diagnostics.R
--- a/analysis/01_qc_diagnostics.R
+++ b/analysis/01_qc_diagnostics.R
@@ -16,13 +16,9 @@
-  if (is.null(assignments_map) || length(assignments_map) == 0) {
-    return(list(
-      status = "skipped",
-      reason = "No assignments mapping configured in input.assignments",
-      outputs = character(0)
-    ))
-  }
-
+
   dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
   reconciliation_rows <- list()
   length_summary_rows <- list()
+  accounting_rows <- list()
+  investigation_rows <- list()
+  median_or_na <- function(x) if (length(x)) stats::median(x) else NA_real_
@@ -45,16 +41,32 @@
   for (sample_id in context$samples) {
     asgn_path <- assignments_map[[sample_id]]
-    if (is.null(asgn_path) || !file.exists(asgn_path)) {
-      next
-    }
-
-    sample_out_dir <- file.path(qc_dir, sanitize_filename(sample_id))
-    dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
-
     # Expected counts from abundance context
     stat_row <- context$sample_stats[context$sample_stats$SampleID == sample_id, ]
     exp_total <- stat_row$TotalReads[1]
     exp_class <- stat_row$ClassifiedReads[1]
     exp_unclass <- stat_row$UnclassifiedReads[1]
+    assignment_available <- !is.null(asgn_path) && file.exists(asgn_path)
+    accounting_rows[[sample_id]] <- data.frame(
+      SampleID = sample_id, AssignmentAvailable = assignment_available,
+      AbundanceTotal = exp_total, AbundanceClassified = exp_class,
+      AbundanceUnclassified = exp_unclass, RawC = NA_integer_, RawU = NA_integer_,
+      C_TaxID0 = NA_integer_, TaxID_GT0 = NA_integer_,
+      EffectiveClassifiedPct = 100 * exp_class / exp_total,
+      C0ShareOfEffectiveUnclassifiedPct = NA_real_, stringsAsFactors = FALSE
+    )
+    investigation_rows[[sample_id]] <- data.frame(
+      SampleID = sample_id, AssignmentAvailable = assignment_available,
+      BamstatsAvailable = FALSE, MedianClassifiedLength = NA_real_,
+      MedianC0Length = NA_real_, MedianRawULength = NA_real_,
+      MinPercentIdentity = context$params$min_percent_identity,
+      MinRefCoverage = context$params$min_ref_coverage,
+      BamstatsC0Matched = NA_integer_, IdentityOnlyFailed = NA_integer_,
+      RefCoverageOnlyFailed = NA_integer_, BothFailed = NA_integer_,
+      stringsAsFactors = FALSE
+    )
+    if (!assignment_available) next
+
+    sample_out_dir <- file.path(qc_dir, sanitize_filename(sample_id))
+    dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
-
+
     reads <- context$assignment_data[[sample_id]]
@@ -73,5 +85,27 @@
     n_qc_reclass <- sum(reads$status == "C" & reads$taxid == 0)
     n_eff_class <- sum(reads$effective_classified)
     n_eff_unclass <- sum(!reads$effective_classified)
+    accounting_rows[[sample_id]][c("RawC", "RawU", "C_TaxID0", "TaxID_GT0")] <-
+      list(n_status_C, n_status_U, n_qc_reclass, n_eff_class)
+    accounting_rows[[sample_id]]$C0ShareOfEffectiveUnclassifiedPct <- if (exp_unclass > 0) {
+      100 * n_qc_reclass / exp_unclass
+    } else {
+      NA_real_
+    }
+    investigation_rows[[sample_id]]$MedianClassifiedLength <-
+      median_or_na(reads$read_length[reads$effective_classified])
+    investigation_rows[[sample_id]]$MedianC0Length <-
+      median_or_na(reads$read_length[reads$status == "C" & reads$taxid == 0])
+    investigation_rows[[sample_id]]$MedianRawULength <-
+      median_or_na(reads$read_length[reads$status == "U"])
+
+    bamstats_path <- context$bamstats[[sample_id]]
+    if (!is.na(bamstats_path)) {
+      partition <- partition_minimap2_failures(reads, bamstats_path, context$params, sample_id)
+      investigation_rows[[sample_id]]$BamstatsAvailable <- TRUE
+      investigation_rows[[sample_id]][c(
+        "BamstatsC0Matched", "IdentityOnlyFailed", "RefCoverageOnlyFailed", "BothFailed"
+      )] <- list(partition$matched, partition$identity_only, partition$coverage_only, partition$both)
+    }
-
+
     # 1. Reconciliation Table Row
@@ -98,3 +132,36 @@
     reads$effective_status <- ifelse(reads$effective_classified, "Classified", "Unclassified")
+
+    accounting_plot <- data.frame(
+      Category = factor(c("Raw U", "C + TaxID 0", "TaxID > 0"),
+                        levels = c("Raw U", "C + TaxID 0", "TaxID > 0")),
+      Count = c(n_status_U, n_qc_reclass, n_eff_class)
+    ) %>%
+      mutate(Fraction = Count / sum(Count),
+             Label = sprintf("%s\n%s (%.1f%%)", Category, scales::comma(Count), 100 * Fraction))
+    p0a <- ggplot(accounting_plot, aes(x = 3, y = Count, fill = Category)) +
+      geom_col(width = 1, color = "white", linewidth = 1.2) +
+      coord_polar(theta = "y") +
+      xlim(c(1, 4)) +
+      geom_text(data = accounting_plot[accounting_plot$Count > 0, , drop = FALSE],
+                aes(x = 3.5, label = Label),
+                position = position_stack(vjust = 0.5), size = 3.1) +
+      annotate("text", x = 1, y = 0,
+               label = sprintf("%s\nreads", scales::comma(exp_total)),
+               size = 4.2, fontface = "bold") +
+      scale_fill_manual(
+        values = c("Raw U" = "#bdbdbd", "C + TaxID 0" = "#e6ab02", "TaxID > 0" = "#1b9e77"),
+        labels = c(
+          "Raw U" = "Never aligned/classified against the reference",
+          "C + TaxID 0" = "Aligned, then failed identity and/or reference-coverage QC",
+          "TaxID > 0" = "Effectively classified; included in classified denominators"
+        )
+      ) +
+      theme_void() +
+      labs(title = sprintf("Exact Read Accounting: %s", sample_id), fill = NULL) +
+      theme(legend.position = "bottom",
+            plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
+    p0a_path <- file.path(sample_out_dir, "00a_read_accounting_donut.png")
+    save_plot(p0a_path, p0a, width = 7.5, height = 6)
+    all_outputs <- c(all_outputs, p0a_path)
-
+
     for (cat_name in c("All", "Classified", "Unclassified", "QC-filtered")) {
@@ -240,2 +307,28 @@
   # Export summary TSVs
+  accounting_df <- do.call(rbind, accounting_rows)
+  stopifnot(all(accounting_df$AbundanceClassified + accounting_df$AbundanceUnclassified ==
+                  accounting_df$AbundanceTotal))
+  available <- accounting_df$AssignmentAvailable
+  stopifnot(all(accounting_df$RawC[available] + accounting_df$RawU[available] ==
+                  accounting_df$AbundanceTotal[available]))
+  stopifnot(all(accounting_df$RawU[available] + accounting_df$C_TaxID0[available] ==
+                  accounting_df$AbundanceUnclassified[available]))
+  stopifnot(all(accounting_df$TaxID_GT0[available] == accounting_df$AbundanceClassified[available]))
+  accounting_file <- file.path(qc_dir, "00_read_accounting.tsv")
+  write.table(accounting_df, accounting_file, sep = "\t", row.names = FALSE, quote = FALSE)
+  all_outputs <- c(all_outputs, accounting_file)
+
+  investigation_df <- do.call(rbind, investigation_rows)
+  with_bamstats <- investigation_df$BamstatsAvailable
+  stopifnot(all(
+    investigation_df$BamstatsC0Matched[with_bamstats] ==
+      investigation_df$IdentityOnlyFailed[with_bamstats] +
+      investigation_df$RefCoverageOnlyFailed[with_bamstats] +
+      investigation_df$BothFailed[with_bamstats]
+  ))
+  investigation_file <- file.path(qc_dir, "00_read_investigation.tsv")
+  write.table(investigation_df, investigation_file, sep = "\t", row.names = FALSE, quote = FALSE,
+              na = "NA")
+  all_outputs <- c(all_outputs, investigation_file)
+
   reconciliation_file <- file.path(qc_dir, "classification_reconciliation.tsv")
```

Keep existing `01a_classification_donut.png` through
`01d_read_length_violin.png` unchanged.

### 3.5 Richness overview

Add a pure helper above `run_alpha()` and call it after `alpha_diversity.tsv`:

```diff
diff --git a/analysis/02_alpha_diversity.R b/analysis/02_alpha_diversity.R
--- a/analysis/02_alpha_diversity.R
+++ b/analysis/02_alpha_diversity.R
@@ -11,1 +11,21 @@
 })
+
+build_richness_overview <- function(class_matrix, samples) {
+  do.call(rbind, lapply(samples, function(sample_id) {
+    positive <- class_matrix[, sample_id][class_matrix[, sample_id] > 0]
+    classified_reads <- sum(positive)
+    singleton_taxa <- sum(positive == 1)
+    low_count <- positive[positive <= 10]
+    data.frame(
+      SampleID = sample_id,
+      ClassifiedReads = classified_reads,
+      PositiveTaxa = length(positive),
+      SingletonTaxa = singleton_taxa,
+      SingletonPct = 100 * singleton_taxa / length(positive),
+      TaxaLeq10 = length(low_count),
+      ReadsInTaxaLeq10 = sum(low_count),
+      ReadsInTaxaLeq10Pct = 100 * sum(low_count) / classified_reads,
+      stringsAsFactors = FALSE
+    )
+  }))
+}
@@ -54,3 +74,8 @@
   alpha_tsv <- file.path(alpha_dir, "alpha_diversity.tsv")
   write.table(alpha_wide, alpha_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, alpha_tsv)
+
+  richness_overview <- build_richness_overview(class_matrix, samples)
+  richness_tsv <- file.path(alpha_dir, "02_richness_overview.tsv")
+  write.table(richness_overview, richness_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
+  all_outputs <- c(all_outputs, richness_tsv)
```

### 3.6 P0 single-sample composition fix

```diff
diff --git a/analysis/04_taxa_composition.R b/analysis/04_taxa_composition.R
--- a/analysis/04_taxa_composition.R
+++ b/analysis/04_taxa_composition.R
@@ -32,4 +32,7 @@
   # Total and classified read denominators per sample
   sample_totals <- colSums(count_matrix)
-  unclass_counts <- count_matrix[unclass_idx, ]
+  unclass_counts <- stats::setNames(
+    context$sample_stats$UnclassifiedReads,
+    context$sample_stats$SampleID
+  )
   class_totals <- sample_totals - unclass_counts
```

### 3.7 Manifest provenance and warning capture

```diff
diff --git a/analysis/00_run_pipeline.R b/analysis/00_run_pipeline.R
--- a/analysis/00_run_pipeline.R
+++ b/analysis/00_run_pipeline.R
@@ -95,2 +95,6 @@
   cat(sprintf("Classified reads: %s\n", format(sum(context$sample_stats$ClassifiedReads), big.mark = ",")))
+  cat(sprintf("Upstream:         %s / %s / rank %s\n",
+              context$upstream_contract$classifier,
+              context$upstream_contract$database_set,
+              context$upstream_contract$taxonomic_rank))
   cat(sprintf("Abundance SHA256: %s\n", context$file_hashes$abundance_table))
@@ -133,6 +137,14 @@
   cat(sprintf("\n>>> Executing module [%s]...\n", mod_name))
   mod_fn <- module_registry[[mod_name]]
   mod_start <- Sys.time()
+  module_warnings <- character(0)
-
+
   mod_res <- tryCatch({
-    mod_fn(context)
+    withCallingHandlers(
+      mod_fn(context),
+      warning = function(w) {
+        module_warnings <<- c(module_warnings, conditionMessage(w))
+        cat(sprintf("WARNING in module [%s]: %s\n", mod_name, conditionMessage(w)), file = stderr())
+        invokeRestart("muffleWarning")
+      }
+    )
@@ -151,2 +163,3 @@
   mod_res$duration_seconds <- as.numeric(difftime(mod_end, mod_start, units = "secs"))
+  mod_res$warnings <- unique(module_warnings)
   module_results[[mod_name]] <- mod_res
@@ -226,2 +239,14 @@
-  } else NULL
+  } else NULL,
+  bamstats = if (any(!is.na(context$bamstats))) {
+    lapply(names(context$bamstats)[!is.na(context$bamstats)], function(sample_id) {
+      path <- context$bamstats[[sample_id]]
+      list(
+        sample_id = sample_id,
+        path = path,
+        size_bytes = file.info(path)$size,
+        mtime = as.character(file.info(path)$mtime),
+        sha256 = context$file_hashes[[paste0("bamstats_", sample_id)]]
+      )
+    })
+  } else NULL
 )
@@ -278,4 +303,8 @@
   cli = cfg$cli,
   inputs = input_meta,
+  upstream_contract = context$upstream_contract,
   modules = module_results,
-  warnings = context$warnings,
+  warnings = unique(c(
+    context$warnings,
+    unlist(lapply(module_results, function(x) x$warnings %||% character(0)), use.names = FALSE)
+  )),
```

### 3.8 Synthetic fixtures and focused tests

Add producer params and gzipped bamstats generators:

```diff
diff --git a/tests/testthat/helper-fixtures.R b/tests/testthat/helper-fixtures.R
--- a/tests/testthat/helper-fixtures.R
+++ b/tests/testthat/helper-fixtures.R
@@ -72,3 +72,30 @@
   write.table(df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
   file_path
 }
+
+create_temp_params <- function(dir, classifier = "minimap2",
+                               database_set = "ncbi_16s_18s", taxonomic_rank = "S") {
+  path <- file.path(dir, "params.json")
+  database_sets <- list()
+  database_sets[[database_set]] <- list(taxonomy = "new_taxdump_2025-01-01.zip")
+  jsonlite::write_json(list(
+    classifier = classifier, database_set = database_set,
+    taxonomic_rank = taxonomic_rank, min_len = 1300, max_len = 1700,
+    min_read_qual = 10, min_percent_identity = 90, min_ref_coverage = 90,
+    abundance_threshold = 1, taxonomy = NULL, reference = NULL,
+    ref2taxid = NULL, database = NULL, database_sets = database_sets,
+    output_unclassified = TRUE, include_read_assignments = TRUE,
+    wf = list(agent = "epi2melabs/test")
+  ), path, auto_unbox = TRUE, pretty = TRUE, null = "null")
+  path
+}
+
+create_temp_bamstats <- function(dir, sample_id, data) {
+  path <- file.path(dir, "bamstats.readstats.tsv.gz")
+  data$sample_name <- sample_id
+  con <- gzfile(path, open = "wt")
+  on.exit(close(con), add = TRUE)
+  write.table(data[, c("name", "sample_name", "iden", "ref_coverage")], con,
+              sep = "\t", row.names = FALSE, quote = FALSE)
+  path
+}
```

Add these tests to `tests/testthat/test-io.R` (and source `config.R` before
`io.R`):

```diff
diff --git a/tests/testthat/test-io.R b/tests/testthat/test-io.R
--- a/tests/testthat/test-io.R
+++ b/tests/testthat/test-io.R
@@ -5,1 +5,33 @@
-source(file.path("..", "..", "analysis", "utils", "io.R"))
+source(file.path("..", "..", "analysis", "utils", "config.R"))
+source(file.path("..", "..", "analysis", "utils", "io.R"))
+
+test_that("producer contract accepts supported minimap2 and rejects unsafe alternatives", {
+  root <- tempfile("params_contract_")
+  dir.create(root)
+  supported <- read_upstream_params(create_temp_params(root))
+  expect_equal(supported$classifier, "minimap2")
+  expect_error(read_upstream_params(create_temp_params(root, classifier = "kraken2")),
+               "supports minimap2 only")
+  expect_error(read_upstream_params(create_temp_params(root, database_set = "SILVA_138_1")),
+               "bundled NCBI database sets only")
+  expect_error(read_upstream_params(create_temp_params(root, taxonomic_rank = "G")),
+               "expected species rank")
+})
+
+test_that("bamstats partition is one-to-one and conserves all C0 reads", {
+  root <- tempfile("bamstats_contract_")
+  dir.create(root)
+  reads <- data.frame(
+    status = c("C", "C", "C", "C"),
+    read_id = c("positive", "identity", "coverage", "both"),
+    taxid = c(123, 0, 0, 0), stringsAsFactors = FALSE
+  )
+  bamstats <- create_temp_bamstats(root, "S1", data.frame(
+    name = reads$read_id, iden = c(95, 89, 95, 89),
+    ref_coverage = c(95, 95, 89, 89)
+  ))
+  params <- read_upstream_params(create_temp_params(root))
+  observed <- partition_minimap2_failures(reads, bamstats, params, "S1")
+  expect_equal(unlist(observed), c(matched = 3, identity_only = 1,
+                                   coverage_only = 1, both = 1))
+})
```

Every synthetic `build_context(cfg)` call must now point to a synthetic params
file. Make these mechanical additions:

```diff
diff --git a/tests/testthat/test-alpha-regression.R b/tests/testthat/test-alpha-regression.R
--- a/tests/testthat/test-alpha-regression.R
+++ b/tests/testthat/test-alpha-regression.R
@@ -74,2 +74,3 @@
     cfg <- get_default_config()
     cfg$input$abundance_table <- abundance
+    cfg$input$params_json <- create_temp_params(root)
diff --git a/tests/testthat/test-cohort.R b/tests/testthat/test-cohort.R
--- a/tests/testthat/test-cohort.R
+++ b/tests/testthat/test-cohort.R
@@ -18,2 +18,3 @@
   cfg <- get_default_config()
   cfg$input$abundance_table <- ab_file
+  cfg$input$params_json <- create_temp_params(tmp)
@@ -55,4 +56,5 @@
   cfg$mode <- "cohort"
   cfg$input$abundance_table <- ab_file
   cfg$input$metadata <- meta_file
+  cfg$input$params_json <- create_temp_params(tmp)
   cfg$output$base_dir <- file.path(tmp, "out_cohort")
@@ -105,4 +107,5 @@
   cfg$mode <- "cohort"
   cfg$input$abundance_table <- ab_file
   cfg$input$metadata <- meta_file
+  cfg$input$params_json <- create_temp_params(tmp)
   expect_error(build_context(cfg), "at least 2")
@@ -122,1 +125,2 @@
   cfg$output$base_dir <- file.path(tmp, "out_underrep")
+  cfg$input$params_json <- create_temp_params(tmp)
@@ -159,2 +163,3 @@
     cfg$input$abundance_table <- abundance
     cfg$input$metadata <- metadata
+    cfg$input$params_json <- create_temp_params(root)
diff --git a/tests/testthat/test-kreport.R b/tests/testthat/test-kreport.R
--- a/tests/testthat/test-kreport.R
+++ b/tests/testthat/test-kreport.R
@@ -62,2 +62,3 @@
   cfg$input$abundance_table <- ab_path
+  cfg$input$params_json <- file.path("..", "..", "output_AAy", "params.json")
   cfg$taxonomy$cache <- cache_path
@@ -136,2 +137,3 @@
   cfg$input$abundance_table <- abundance
+  cfg$input$params_json <- create_temp_params(root)
   cfg$input$assignments <- list(S1 = assignments)
```

Update the process fixture so `params.json` cannot be omitted and add the Kraken
rejection regression:

```diff
diff --git a/tests/testthat/test-release-process.R b/tests/testthat/test-release-process.R
--- a/tests/testthat/test-release-process.R
+++ b/tests/testthat/test-release-process.R
@@ -11,3 +11,3 @@
 write_process_config <- function(root, sample_names = "S1", mode = "auto",
                                  metadata = NULL, assignments = NULL,
-                                 unresolved_policy = "warn") {
+                                 unresolved_policy = "warn", classifier = "minimap2") {
@@ -21,1 +21,3 @@
-  cfg$input$params_json <- NULL
+  cfg$input$params_json <- normalizePath(
+    create_temp_params(root, classifier = classifier), winslash = "/"
+  )
@@ -90,3 +92,18 @@
   expect_gt(invalid_metadata$status, 0L)
   expect_match(invalid_metadata$stderr, "Metadata SampleID mismatch")
 })
+
+test_that("Kraken2 is rejected from params before assignment schema parsing", {
+  root <- tempfile("kraken_contract_")
+  dir.create(root)
+  bad_assignment <- file.path(root, "six_field.tsv")
+  writeLines("C\tread1\t123\t1500\tA:1\tBacteria", bad_assignment)
+  config <- write_process_config(
+    root, assignments = list(S1 = normalizePath(bad_assignment, winslash = "/")),
+    classifier = "kraken2"
+  )
+  result <- run_pipeline_process(c("--config", config, "--validate-only"), tempdir())
+  expect_gt(result$status, 0L)
+  expect_match(result$stderr, "supports minimap2 only")
+  expect_false(grepl("expected exactly 5", result$stderr))
+})
```

Add the richness regression to `tests/testthat/test-alpha-regression.R`:

```diff
diff --git a/tests/testthat/test-alpha-regression.R b/tests/testthat/test-alpha-regression.R
--- a/tests/testthat/test-alpha-regression.R
+++ b/tests/testthat/test-alpha-regression.R
@@ -64,3 +64,17 @@
   expect_equal(length(unique(tax_df$genus)), 867)
   expect_equal(length(unique(tax_df$species)), 1836)
 })
+
+test_that("Ambar richness overview exposes the exact low-count tail", {
+  ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
+  skip_if_not(file.exists(ab_path), "Real abundance table not found")
+  parsed <- read_abundance_table(ab_path)
+  class_matrix <- parsed$count_matrix[-parsed$unclass_index, , drop = FALSE]
+  observed <- build_richness_overview(class_matrix, parsed$samples)
+  expect_equal(observed$ClassifiedReads, 80556)
+  expect_equal(observed$PositiveTaxa, 1836)
+  expect_equal(observed$SingletonTaxa, 735)
+  expect_equal(observed$TaxaLeq10, 1399)
+  expect_equal(observed$ReadsInTaxaLeq10, 3456)
+  expect_equal(observed$ReadsInTaxaLeq10Pct, 4.2901832265753, tolerance = 1e-10)
+})
```

### 3.9 Full supported-fixture assertions

```diff
diff --git a/tests/verify_release_run.R b/tests/verify_release_run.R
--- a/tests/verify_release_run.R
+++ b/tests/verify_release_run.R
@@ -13,1 +13,1 @@
-stopifnot(identical(manifest$pipeline_version, "0.1.0"))
+stopifnot(identical(manifest$pipeline_version, "0.2.0"))
@@ -17,1 +17,8 @@
 stopifnot(identical(manifest$cli$refresh_taxonomy, FALSE))
+stopifnot(identical(manifest$upstream_contract$classifier, "minimap2"))
+stopifnot(identical(manifest$upstream_contract$database_set, "ncbi_16s_18s_28s_ITS"))
+stopifnot(identical(manifest$upstream_contract$taxonomic_rank, "S"))
+stopifnot(identical(manifest$upstream_contract$wf_agent, "epi2melabs/5.2.5"))
+stopifnot(is.null(manifest$upstream_contract$workflow_version))
+stopifnot(is.null(manifest$upstream_contract$workflow_revision))
+stopifnot(length(manifest$warnings) == 0L)
@@ -29,1 +36,4 @@
   "run_manifest.json",
+  "01_QC/00_read_accounting.tsv",
+  "01_QC/00_read_investigation.tsv",
+  "01_QC/AmbarAyunda_minimap2_16S/00a_read_accounting_donut.png",
@@ -32,1 +42,2 @@
   "02_Alpha_Diversity/alpha_diversity.tsv",
+  "02_Alpha_Diversity/02_richness_overview.tsv",
@@ -52,1 +63,36 @@
 stopifnot(sum(reconciliation$AbundanceUnclassified) == 33500L)
+
+accounting <- read.delim(file.path(root, "01_QC/00_read_accounting.tsv"), check.names = FALSE)
+stopifnot(nrow(accounting) == 1L)
+stopifnot(identical(accounting$SampleID, "AmbarAyunda_minimap2_16S"))
+stopifnot(accounting$AbundanceTotal == 114056L)
+stopifnot(accounting$RawC == 89809L, accounting$RawU == 24247L)
+stopifnot(accounting$C_TaxID0 == 9253L, accounting$TaxID_GT0 == 80556L)
+stopifnot(abs(accounting$EffectiveClassifiedPct - 70.62756540647) < 1e-10)
+stopifnot(abs(accounting$C0ShareOfEffectiveUnclassifiedPct - 27.62089552239) < 1e-10)
+
+investigation <- read.delim(file.path(root, "01_QC/00_read_investigation.tsv"), check.names = FALSE)
+stopifnot(investigation$MedianClassifiedLength == 1507)
+stopifnot(investigation$MedianC0Length == 1493)
+stopifnot(investigation$MedianRawULength == 1494)
+stopifnot(identical(investigation$BamstatsAvailable, FALSE))
+stopifnot(is.na(investigation$BamstatsC0Matched))
+stopifnot(is.na(investigation$IdentityOnlyFailed))
+stopifnot(is.na(investigation$RefCoverageOnlyFailed))
+stopifnot(is.na(investigation$BothFailed))
+
+richness <- read.delim(file.path(root, "02_Alpha_Diversity/02_richness_overview.tsv"),
+                       check.names = FALSE)
+stopifnot(richness$ClassifiedReads == 80556L)
+stopifnot(richness$PositiveTaxa == 1836L, richness$SingletonTaxa == 735L)
+stopifnot(richness$TaxaLeq10 == 1399L, richness$ReadsInTaxaLeq10 == 3456L)
+stopifnot(abs(richness$ReadsInTaxaLeq10Pct - 4.2901832265753) < 1e-10)
+
+composition <- read.delim(
+  file.path(root, "04_Taxa_Composition/classification_fraction.tsv"), check.names = FALSE
+)
+stopifnot(!anyNA(composition))
+stopifnot(composition$TotalReads == 114056L)
+stopifnot(composition$ClassifiedReads == 80556L)
+stopifnot(composition$UnclassifiedReads == 33500L)
+stopifnot(composition$ClassifiedReads + composition$UnclassifiedReads == composition$TotalReads)
```

### 3.10 User-facing documentation

```diff
diff --git a/README.md b/README.md
--- a/README.md
+++ b/README.md
@@ -5,1 +5,1 @@
-Current pipeline version: **0.1.0**.
+Current pipeline version: **0.2.0**.
@@ -47,1 +47,1 @@
-Version 0.1.0 is validated against the tracked `wf-16s` minimap2 abundance,
+Version 0.2.0 is validated against the tracked `wf-16s` minimap2 abundance,
@@ -69,1 +69,13 @@
 - **Effective Classification**: Defined as `taxid > 0`. Raw status `C` alone overstates classification; reads with status `C` and TaxID `0` are filtered out by upstream reference coverage/identity thresholds.
+
+### 3. Producer contract and optional bamstats
+
+- `input.params_json` is required. Version 0.2.0 accepts minimap2, species-rank
+  output from the bundled NCBI database sets without custom reference or
+  taxonomy overrides.
+- `input.wf16s_output_root` is optional. When supplied, the pipeline discovers
+  each sample's `bamstats.readstats.tsv.gz` by its internal `sample_name` and
+  partitions `C + TaxID 0` reads into identity-only, coverage-only, and
+  both-threshold failures.
+- Missing optional bamstats produces explicit `NA` fields; a configured invalid
+  or ambiguous bamstats contract fails rather than selecting an arbitrary file.
@@ -88,1 +100,18 @@
   - Subsampling iterations with `vegan::rrarefy()` represent **rarefaction resamples without replacement** to assess community sampling sensitivity. They are never presented as independent biological replicates.
+- **Richness overview**:
+  - `02_richness_overview.tsv` reports positive taxa, singletons, taxa with at
+    most 10 reads, and the reads represented by that low-count tail.
+  - These summaries remain conditional on upstream filtering, database choice,
+    taxonomy, and the abundance threshold; low count alone does not prove artifact.
+
+### Unclassified-read handling by module
+
+| Module | Handling |
+|---|---|
+| QC | Included in abundance accounting; all assignment reads are shown when assignments exist |
+| Alpha diversity | Excluded |
+| Beta diversity | Excluded |
+| Taxonomic composition | Excluded from rank tables; included in the classification-fraction output |
+| Ordination | Excluded |
+| Shared taxa | Excluded |
+| Kraken/Pavian | Included as the first `U` line and in total-read arithmetic |
@@ -173,1 +202,1 @@
-Every run records input hashes, package versions, interpreter versions, and
+Every run records input hashes, the validated upstream contract, package versions, interpreter versions, and
@@ -175,1 +204,1 @@
-an environment. Version 0.1.0 is therefore provenance-captured, not bitwise
+an environment. Version 0.2.0 is therefore provenance-captured, not bitwise
@@ -199,1 +228,1 @@
-Version 0.1.0 covers validated parsing, single-sample summaries, synthetically
+Version 0.2.0 covers validated NCBI/minimap2/species parsing, exact read accounting, single-sample summaries, synthetically
```

```diff
diff --git a/CHANGELOG.md b/CHANGELOG.md
--- a/CHANGELOG.md
+++ b/CHANGELOG.md
@@ -3,1 +3,23 @@
 All notable changes to this project are documented in this file.
+
+## [0.2.0] - Unreleased
+
+### Added
+
+- Exact per-sample read accounting with a three-way diagnostic donut.
+- Optional one-to-one bamstats joins for identity-only, reference-coverage-only,
+  and dual-threshold minimap2 failures.
+- Classified-only richness summaries for singleton and low-count tails.
+- Structured upstream producer-contract and bamstats provenance.
+
+### Fixed
+
+- Preserve the single-sample unclassified count in composition outputs instead
+  of producing `NA` and a misleading 100%-classified plot.
+- Reject unsupported Kraken2, SILVA, non-species, and custom-reference contracts
+  before classifier-specific assignment parsing.
+
+### Documentation
+
+- State each module's treatment of unclassified reads and the evidentiary limit
+  of low-count richness summaries.
```

## 4. Verification gates

Gemini should implement in three commits or clearly separated phases, stopping
on the first red gate:

1. P0 composition and producer contract.
2. Accounting/bamstats/richness outputs plus unit tests.
3. Manifest, documentation, and release metadata.

Run from the repository root:

```bash
Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
python -m compileall -q analysis tests
Rscript tests/testthat.R
python -m unittest -v tests/test_ncbi_taxonomy.py
Rscript analysis/00_run_pipeline.R --config config.yml --output-dir "output_v020_check"
Rscript tests/verify_release_run.R "output_v020_check"
git diff --check
```

Before the full run, use a fresh output directory. Do not use `--overwrite` on
an older output tree because stale artifacts can satisfy presence checks.

Manual local-fixture checks after setting `input.wf16s_output_root` to the
corresponding local producer output (not CI requirements):

| Fixture | Accounting / medians / threshold partitions | Richness tail |
|---|---|---|
| Ambar NCBI | `114056 / 89809 / 24247 / 9253 / 80556`; `1507 / 1493 / 1494`; `1608 / 5490 / 2155` | `1836`, `735`, `1399`, `3456`, `4.290183%` |
| Helga NCBI | `23099 / 23007 / 92 / 4852 / 18155`; `1510 / 1471.5 / 1468.5`; `972 / 3210 / 670` | `756`, `348`, `658`, `1496`, `8.240154%` |
| July NCBI | `17920 / 17785 / 135 / 8117 / 9668`; `1474 / 1472 / 1411`; `3536 / 2863 / 1718` | `678`, `276`, `555`, `1398`, `14.460074%` |

Also prove:

- an abundance+params run without assignments writes both `00_` tables with
  abundance fields populated and assignment/bamstats fields `NA`;
- a configured existing root with no bamstats logs an informational message;
- two bamstats candidates with the same internal `sample_name` fail;
- a missing status-C read, duplicate bamstats `name`, malformed numeric metric,
  C0 read passing both thresholds, or TaxID-positive read failing a threshold
  fails with an actionable message;
- Kraken2 fails from `params.json` before the six-field assignment parser;
- `--validate-only` remains mutation-free; and
- the final Git diff contains only approved tracked files. No commit, tag, or
  push is authorized by this plan alone.

## 5. Explicit non-goals

- SILVA taxonomy support or seven-rank/genus schemas.
- Kraken2/Bracken parsing or denominator reconciliation.
- Parsing `nextflow.log` to recover wf-16s version/revision.
- Reclassifying, discarding, or rerunning unclassified reads.
- Treating low-count taxa as contaminants or errors without independent evidence.
- Adding the local producer outputs to Git.
