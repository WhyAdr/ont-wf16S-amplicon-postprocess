# Patch plan: `ont-wf16S-amplicon-postprocess`

**Handoff target:** Gemini
**Repository:** `WhyAdr/ont-wf16S-amplicon-postprocess`
**Audited base commit:** `3638b7665db14830cc1e97cdfaa18c5660f3528b` (`main`, 2026-09-03)
**Scope:** Phase 1 correctness, fail-closed validation, taxonomy safety, statistical gating, tests, and CI
**Overall verdict:** promising architecture, but the audited commit is not yet safe to call complete or publication-ready.

## 1. Executive decision

Keep the modular architecture. The shared context, config-relative paths, explicit single/cohort modes, classified-only denominators, and separation of rarefaction resamples from biological replicates are all good design choices.

Do not release the current commit unchanged. One blocker is absolute: `%||%` is called throughout the runner/modules but is never defined, so normal module execution fails. Several other paths either fail or silently violate the written contracts, especially `--validate-only`, cohort alpha testing, taxonomy refresh, output overwrite protection, zero-distance ordination, and cohort taxonomy reporting.

The exact patch in Section 5 is based on the audited commit and was checked with `git apply --check`. It deliberately also removes trailing whitespace from the new R/Python/test files because the implementation commit fails its own `git diff --check` gate.

## 2. Prioritized findings

| ID | Severity | Finding and evidence at audited commit | Consequence | Resolution in patch |
|---|---|---|---|---|
| F01 | P0 | `%||%` is used in every analytical module and the runner, but no definition exists anywhere in the repository. | The first executed module reaches “could not find function `%||%`”; the test suite's cohort module calls fail for the same reason. | Define it once in `analysis/utils/config.R`, which is sourced before modules. |
| F02 | P0 | `build_context()` stores assignment paths but does not parse them. `--validate-only` exits before `run_qc()`, while QC silently `next`s over absent files. | The advertised zero-mutation validation can report PASS for missing, malformed, duplicated, or unreconciled assignment files. | Strictly parse and reconcile every supplied mapping in context construction; retain full rows only when QC needs them. |
| F03 | P0 | `run_alpha()` merges metadata into `alpha_wide`, then merges the same metadata again before group tests. | In cohort mode, `Group` becomes `Group.x`/`Group.y`; `table(alpha_with_group$Group)` fails. | Join metadata once, preserve canonical sample order, and test the cohort path. |
| F04 | P0 | The taxonomy module passes only the first sample's assignment file. Assignment-derived leaf TaxIDs are added only to Python memory in `cache_only` mode, then R reloads the unchanged source cache. | Other cohort samples are ignored and newly derived TaxIDs disappear before `.kreport` generation. | Accept repeated assignment inputs and always write/load a run-local resolved cache without mutating the offline source cache. |
| F05 | P0 | The Python resolver computes unresolved nodes before online refresh and never recomputes them. It also treats an unreadable cache as empty and may replace it in refresh mode. | A successful refresh can still fail `unresolved_policy: error`; a corrupt/partially readable cache can be overwritten. | Recompute after refresh, fail closed on invalid caches, and preserve the source cache when any request fails. |
| F06 | P1 | Forced `mode: cohort` has no `>=2` sample gate; metadata accepts blank/NA `SampleID` or `Group`; configured strata can silently degrade to unrestricted permutations. | Invalid experimental designs can enter inferential modules. | Add mode, metadata, configuration, and strata validation. |
| F07 | P1 | `--overwrite` checks only `run_manifest.json`, after context work; module writers overwrite known artifacts directly. | A previous partial run without a manifest can be overwritten despite the default safety promise. | Preflight all known output roots before creating/writing anything. |
| F08 | P1 | Cohort alpha ignores configured `resample_depth`; every sample resets to the same RNG seed; groups with inadequate replication can be silently discarded while remaining groups are tested. | Configuration is ignored, resampling streams are coupled, and the tested population can differ from the declared cohort. | Cap the common depth with the configured depth, derive sample-keyed seeds, and require every compared group to pass replication gates. |
| F09 | P1 | Distance matrices are written with row names as an unnamed first column; Jaccard ignores a configurable minimum count; PCoA/PERMANOVA lack zero-distance gates; beta resampling config is a silent no-op. | Output violates the tabular contract, rare detections can dominate Jaccard, and degenerate cohorts can error or yield meaningless inference. | Export explicit `SampleID`, add `beta.minimum_count`, gate zero distances, additive-correct PCoA, and implement Procrustes-aligned rarefaction stability output. |
| F10 | P1 | PCA does not explicitly respect `min(2, n-1, matrix_rank)`; identical samples are not gated; NMDS checks sample count but not unique profiles and omits tries/skip diagnostics. | Degenerate matrices can generate errors/NaNs, and failed NMDS can leave no auditable explanation. | Rank-gate PCA, require unique profiles, and always write NMDS diagnostics. |
| F11 | P1 | Composition lacks the promised abundance-derived all-read classification fraction; display labels can collide across different lineage prefixes. | Classification yield disappears when assignments are absent; duplicate labels can be merged or break factor/heatmap handling. | Add classification table/plot and make display labels unique while retaining `TaxonPath` as the analytical key. |
| F12 | P1 | `resolution_rows[[p]]` overwrites the same taxonomy path for each later sample; kreport counts are coerced to 32-bit R integers. | Cohort resolution output falsely contains only the last sample for shared nodes; very large libraries can overflow. | Append rows and keep whole-number counts as doubles formatted without decimals. |
| F13 | P1 | No config schema validation rejects typos/range errors. The exposed optional beta-resampling configuration is not honored. | Misspelled or nonsensical values can be silently accepted, undermining reproducibility. | Reject unknown keys; validate enums, ranges, booleans, and paths; implement the exposed resampling option. |
| F14 | P1 | The implementation commit produces hundreds of whitespace errors under `git diff --check`; no CI is present. | A stated acceptance gate already fails, and regressions are easy to merge. | Normalize whitespace and add network-free Python/R validation in GitHub Actions. |

Empirical audit notes from the tracked real fixture:

- 1,837 abundance rows, one sample, and 114,056 total reads were confirmed independently.
- The tracked cache has 3,220 entries: 3,130 positive TaxIDs and 90 zero entries.
- The existing resolver reports 46 unresolved nodes and 26 exact-lineage multi-TaxID conflicts on the supplied assignments. These must remain visible; “resolved enough to draw a Sankey” is not the same as taxonomy-complete.
- Current upstream wf-16S still declares `reads_assignments/{{ alias }}.*.assignments.tsv` as a per-sample output, supporting the mapping contract used here: [wf-16S output definition](https://github.com/epi2me-labs/wf-16s/blob/master/output_definition.json).

## 3. Scientific cautions beyond the mechanical patch

These should be tracked as the next milestone rather than silently folded into Phase 1:

1. **Pin taxonomy to the upstream database snapshot.** The tracked `params.json` points to a dated 2025 taxonomy archive. Live 2026 NCBI refresh can create a hybrid taxonomy. Prefer resolving against the exact upstream `taxdump`; if live refresh remains available, label the result as mixed-snapshot and record both snapshots.
2. **Normalize leaf TaxIDs to declared rank.** The 26 conflicts often reflect species-versus-strain/accession-level identifiers. Before majority voting, walk each TaxID to the species ancestor using the pinned taxdump. Otherwise a row marked Kraken rank `S` may carry a below-species TaxID.
3. **Add design formulas, not only `Group`.** `strata` constrains permutations; it does not adjust a batch effect. A future config should support an explicit, validated formula such as `~ Batch + Group`, while detecting rank deficiency/confounding.
4. **Report effect sizes and uncertainty.** Alpha tests should add group medians, Hodges-Lehmann or rank-biserial effects where appropriate, and intervals. PERMANOVA should foreground R² and dispersion results, not p-values alone.
5. **Keep read rarefaction in its lane.** Rarefaction curves/resamples are sensitivity tools, not biological replication or a universal normalization strategy. Classification efficiency and upstream abundance filtering can alter classified-only richness between samples.
6. **Validate on a genuine cohort.** The synthetic 2×3 fixture validates software control flow only. Add at least one de-identified real cohort plus negative/blank controls before calling the comparative path biologically validated.
7. **Lock the final R environment.** After the patched suite passes on the intended workstation, create and commit `renv.lock`. CI alone does not freeze numerical behavior across future package releases.
8. **Choose a repository license.** No license is present. Add one only after the owner selects the intended reuse terms.

## 4. Application and acceptance sequence

1. Confirm the checkout is exactly the audited base or rebase the diff deliberately:

   ```bash
   git rev-parse HEAD
   # expected: 3638b7665db14830cc1e97cdfaa18c5660f3528b
   ```

2. Copy only the fenced diff in Section 5 to `phase1-correctness.patch`, then run:

   ```bash
   git apply --check phase1-correctness.patch
   git apply phase1-correctness.patch
   git diff --check
   ```

3. Run the full gate matrix in an R-capable environment:

   ```bash
   Rscript analysis/install_packages.R
   Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
   python -m compileall -q analysis tests
   Rscript tests/testthat.R
   python -m unittest -v tests/test_ncbi_taxonomy.py
   Rscript analysis/00_run_pipeline.R --config config.yml --validate-only
   ```

4. Run integration into a new temporary output root first:

   ```bash
   Rscript analysis/00_run_pipeline.R --config config.yml --output-dir /tmp/ont-wf16s-phase1
   ```

5. Inspect at minimum:

   - `run_manifest.json` and `resolved_config.yml`;
   - 80,556 classified + 33,500 unclassified = 114,056 total;
   - all Section 2.3 alpha regression targets;
   - 46 currently unresolved nodes and 26 conflicts remain explicitly reported unless the pinned taxonomy strategy resolves them;
   - `.kreport` arithmetic, indentation, rank codes, and integer-formatted counts;
   - no source cache mutation in `cache_only` mode;
   - non-zero exit on malformed assignments, invalid metadata, module failure, and simulated refresh failure.

6. Do not commit generated `output/`. Review the exact diff and test logs before committing. A live NCBI smoke test must stay opt-in and separate from default CI.

### Verification already completed for this handoff

- GitHub-connected repository state and commit were verified.
- The patch applies cleanly to commit `3638b76` with `git apply --check`.
- `git diff --check` passes after application.
- Python source compilation passes.
- The added offline-cache and simulated-refresh-failure tests pass.
- The patched resolver completes against the real tracked fixture offline, reporting 46 unresolved nodes and 26 conflicts without mutating the source cache.
- Delimiter balance was checked across all R sources.

The full R suite was **not executable in the audit container because `Rscript` is absent**. Gemini must treat the R commands above—and the new CI job—as mandatory before calling the patch complete.

## 5. Exact unified diff

~~~~diff
diff --git a/.github/workflows/ci.yml b/.github/workflows/ci.yml
new file mode 100644
index 0000000..03b7de3
--- /dev/null
+++ b/.github/workflows/ci.yml
@@ -0,0 +1,28 @@
+name: validation
+
+on:
+  push:
+  pull_request:
+
+jobs:
+  test:
+    runs-on: ubuntu-latest
+    steps:
+      - uses: actions/checkout@v4
+      - uses: r-lib/actions/setup-r@v2
+        with:
+          use-public-rspm: true
+      - name: Install R dependencies
+        run: Rscript analysis/install_packages.R --install
+      - name: Parse R and Python sources
+        run: |
+          Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
+          python -m compileall -q analysis tests
+      - name: Run unit and regression tests
+        run: |
+          Rscript tests/testthat.R
+          python -m unittest -v tests/test_ncbi_taxonomy.py
+      - name: Validate the tracked real fixture without writes or network
+        run: Rscript analysis/00_run_pipeline.R --config config.yml --validate-only
+      - name: Check whitespace in the proposed change
+        run: git diff --check "${{ github.event.pull_request.base.sha || github.event.before }}" HEAD
diff --git a/README.md b/README.md
index 771c58e..cd282d3 100644
--- a/README.md
+++ b/README.md
@@ -6,7 +6,7 @@ Modular, config-driven downstream post-processing, statistical analysis, diversi

 ## Architecture & Overview

-This repository transforms primary outputs from ONT's Nextflow-based `wf-16s` pipeline into reproducible, publication-ready statistical figures, tables, and Pavian-compatible Kraken reports (`.kreport`).
+This repository transforms primary outputs from ONT's Nextflow-based `wf-16s` pipeline into reproducible statistical figures, tables, and Pavian-compatible Kraken reports (`.kreport`). Treat species-level calls and richness estimates as conditional on the upstream classifier, reference database, and abundance threshold—not as independent validation of organism presence.

 ```
 analysis/
@@ -90,7 +90,9 @@ The pipeline supports two execution modes:
   6. Indented scientific name (2 spaces per depth level)
 - **NCBI Taxonomy Resolution**:
   - Offline default (`network_mode: cache_only`): resolves TaxIDs using local `taxonomy_cache.json` without internet requests.
+  - Assignment-derived TaxIDs are written to a run-local resolved cache; an offline run does not mutate the configured source cache.
   - Opt-in refresh (`--refresh-taxonomy`): queries NCBI Entrez E-utilities with bounded exponential backoff, rate pacing, and atomic cache file replacement.
+  - Unresolved and conflicting mappings are exported explicitly. Named nodes with TaxID `0` should be resolved before treating the report as taxonomy-complete in Pavian.
 - **Interactive Sankey Visualization**:
   - Generated `.kreport` files can be uploaded to [Pavian](https://fbreitwieser.shinyapps.io/pavian/) for interactive Sankey and sunburst diagrams.

@@ -121,6 +123,18 @@ Rscript analysis/00_run_pipeline.R --config config.yml
 Rscript analysis/00_run_pipeline.R --config config.yml --output-dir output --overwrite
 ```

+### Validation
+
+```bash
+Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
+python -m compileall -q analysis tests
+Rscript tests/testthat.R
+python -m unittest -v tests/test_ncbi_taxonomy.py
+Rscript analysis/00_run_pipeline.R --config config.yml --validate-only
+```
+
+The included cohort fixture is synthetic and verifies software behavior only. It is not biological validation of cohort statistics or species-level classification accuracy. For reproducible publication work, record the wf-16S/database versions and create an `renv.lock` from the R environment used for the final analysis.
+
 ### Command-Line Options
 | Option | Description |
 |---|---|
diff --git a/analysis/00_run_pipeline.R b/analysis/00_run_pipeline.R
index 0c88e15..5586775 100644
--- a/analysis/00_run_pipeline.R
+++ b/analysis/00_run_pipeline.R
@@ -17,8 +17,11 @@ get_script_dir <- function() {
 script_dir <- get_script_dir()
 repo_root <- normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE)

-# Source all utilities
+# Bootstrap dependency checking before sourcing files that attach packages.
 source(file.path(script_dir, "utils", "dependencies.R"))
+check_dependencies()
+
+# Source all remaining utilities
 source(file.path(script_dir, "utils", "cli.R"))
 source(file.path(script_dir, "utils", "config.R"))
 source(file.path(script_dir, "utils", "io.R"))
@@ -37,21 +40,19 @@ source(file.path(script_dir, "07_kreport_pavian.R"))

 start_time <- Sys.time()

-# 1. Dependency check
-check_dependencies()
-
-# 2. Parse CLI options
+# 1. Parse CLI options
 cli_opts <- parse_cli_args()

-# 3. Load and resolve configuration
+# 2. Load and resolve configuration
 cfg <- tryCatch({
   load_config(cli_opts$config, cli_opts = cli_opts)
 }, error = function(e) {
   cat(sprintf("[FATAL] Configuration error: %s\n", e$message), file = stderr())
   quit(status = 1)
 })
+cfg$pipeline_root <- repo_root

-# 4. Build and validate shared context
+# 3. Build and validate shared context
 context <- tryCatch({
   build_context(cfg)
 }, error = function(e) {
@@ -59,7 +60,27 @@ context <- tryCatch({
   quit(status = 1)
 })

-# 5. Handle --validate-only
+# Module registry and request validation must happen before any output mutation.
+module_registry <- list(
+  qc          = run_qc,
+  alpha       = run_alpha,
+  beta        = run_beta,
+  composition = run_taxa_composition,
+  ordination  = run_ordination,
+  shared      = run_shared_taxa,
+  kreport     = run_kreport
+)
+
+requested_modules <- cfg$cli$modules
+invalid_modules <- setdiff(requested_modules, names(module_registry))
+if (length(invalid_modules) > 0) {
+  cat(sprintf("[FATAL] Unknown module(s) requested: %s\nAvailable: %s\n",
+              paste(invalid_modules, collapse = ", "),
+              paste(names(module_registry), collapse = ", ")), file = stderr())
+  quit(status = 1)
+}
+
+# 4. Handle --validate-only
 if (cfg$cli$validate_only) {
   cat("=== ONT wf-16s Pipeline Validation Check ===\n")
   cat(sprintf("Config file:      %s\n", cfg$config_file))
@@ -72,11 +93,20 @@ if (cfg$cli$validate_only) {
   quit(status = 0)
 }

-# 6. Overwrite check
-if (!cfg$cli$overwrite && file.exists(cfg$output$manifest_file)) {
+# 5. Overwrite check: protect every known module output, including partial runs
+# that failed before a manifest could be written.
+known_existing <- c(
+  c(cfg$output$manifest_file, cfg$output$resolved_config_file, cfg$output$session_info_file)[
+    file.exists(c(cfg$output$manifest_file, cfg$output$resolved_config_file, cfg$output$session_info_file))
+  ],
+  unlist(lapply(cfg$output$dirs, function(path) {
+    if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE) else character(0)
+  }), use.names = FALSE)
+)
+if (!cfg$cli$overwrite && length(known_existing) > 0L) {
   cat(sprintf(
-    "[FATAL] Output file already exists: '%s'\nUse --overwrite to allow replacing existing outputs.\n",
-    cfg$output$manifest_file
+    "[FATAL] Refusing to overwrite %d existing pipeline output(s); first path: '%s'\nUse --overwrite to allow replacement.\n",
+    length(known_existing), known_existing[1]
   ), file = stderr())
   quit(status = 1)
 }
@@ -84,26 +114,6 @@ if (!cfg$cli$overwrite && file.exists(cfg$output$manifest_file)) {
 # Create base output directory
 dir.create(cfg$output$base_dir, recursive = TRUE, showWarnings = FALSE)

-# Module registry
-module_registry <- list(
-  qc          = run_qc,
-  alpha       = run_alpha,
-  beta        = run_beta,
-  composition = run_taxa_composition,
-  ordination  = run_ordination,
-  shared      = run_shared_taxa,
-  kreport     = run_kreport
-)
-
-requested_modules <- cfg$cli$modules
-invalid_modules <- setdiff(requested_modules, names(module_registry))
-if (length(invalid_modules) > 0) {
-  cat(sprintf("[FATAL] Unknown module(s) requested: %s\nAvailable: %s\n",
-              paste(invalid_modules, collapse = ", "),
-              paste(names(module_registry), collapse = ", ")), file = stderr())
-  quit(status = 1)
-}
-
 cat("=============================================================================\n")
 cat(sprintf("ONT wf-16s Amplicon Post-Processing Pipeline\n"))
 cat(sprintf("Project: %s | Mode: %s | Samples: %d\n",
@@ -118,7 +128,7 @@ for (mod_name in requested_modules) {
   cat(sprintf("\n>>> Executing module [%s]...\n", mod_name))
   mod_fn <- module_registry[[mod_name]]
   mod_start <- Sys.time()
-
+
   mod_res <- tryCatch({
     mod_fn(context)
   }, error = function(e) {
@@ -129,11 +139,13 @@ for (mod_name in requested_modules) {
       outputs = character(0)
     )
   })
-
+
   mod_end <- Sys.time()
+  mod_res$start_time <- format(mod_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
+  mod_res$end_time <- format(mod_end, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
   mod_res$duration_seconds <- as.numeric(difftime(mod_end, mod_start, units = "secs"))
   module_results[[mod_name]] <- mod_res
-
+
   if (mod_res$status == "failed") {
     any_failed <- TRUE
     if (!cfg$cli$keep_going) {
@@ -151,20 +163,20 @@ for (mod_name in requested_modules) {
 end_time <- Sys.time()
 overall_status <- if (any_failed) "failed" else "completed"

-# 7. Write session info
-sink(cfg$output$session_info_file)
-cat("=== System & Interpreter ===\n")
-cat(sprintf("R version: %s\n", R.version.string))
-cat(sprintf("Platform:  %s\n", R.version$platform))
-cat(sprintf("Run time:  %s to %s\n\n", start_time, end_time))
-cat("=== Package Versions ===\n")
 deps <- get_dependency_versions()
-for (pkg in names(deps)) {
-  cat(sprintf("  %-15s: %s\n", pkg, deps[[pkg]]))
-}
-cat("\n=== Full sessionInfo() ===\n")
-print(sessionInfo())
-sink()
+session_lines <- c(
+  "=== System & Interpreter ===",
+  sprintf("R version: %s", R.version.string),
+  sprintf("Platform:  %s", R.version$platform),
+  sprintf("Run time:  %s to %s", start_time, end_time),
+  "",
+  "=== Package Versions ===",
+  sprintf("  %-15s: %s", names(deps), deps),
+  "",
+  "=== Full sessionInfo() ===",
+  capture.output(print(sessionInfo()))
+)
+writeLines(session_lines, cfg$output$session_info_file)

 # 8. Write resolved config
 yaml::write_yaml(cfg, cfg$output$resolved_config_file)
@@ -194,9 +206,28 @@ input_meta <- list(
     size_bytes = file.info(cfg$taxonomy$cache)$size,
     mtime = as.character(file.info(cfg$taxonomy$cache)$mtime),
     sha256 = context$file_hashes$taxonomy_cache
-  ) else NULL
+  ) else NULL,
+  assignments = if (length(context$assignments) > 0L) {
+    lapply(names(context$assignments), function(sample_id) {
+      path <- context$assignments[[sample_id]]
+      list(
+        sample_id = sample_id,
+        path = path,
+        size_bytes = file.info(path)$size,
+        mtime = as.character(file.info(path)$mtime),
+        sha256 = context$file_hashes[[paste0("assignment_", sample_id)]]
+      )
+    })
+  } else NULL
 )

+unresolved_file <- file.path(cfg$output$dirs$kreport, "unresolved_taxids.tsv")
+unresolved_count <- if (file.exists(unresolved_file)) {
+  max(0L, length(readLines(unresolved_file, warn = FALSE)) - 1L)
+} else {
+  NA_integer_
+}
+
 manifest <- list(
   pipeline = "ont-wf16s-postprocess",
   schema_version = cfg$schema_version,
@@ -208,8 +239,18 @@ manifest <- list(
   mode = context$mode,
   seed = cfg$seed,
   samples = context$samples,
+  config_file = cfg$config_file,
+  output_root = cfg$output$base_dir,
+  cli = cfg$cli,
   inputs = input_meta,
   modules = module_results,
+  warnings = context$warnings,
+  taxonomy = list(
+    network_mode = cfg$taxonomy$network_mode,
+    unresolved_policy = cfg$taxonomy$unresolved_policy,
+    unresolved_count = unresolved_count
+  ),
+  interpreter = list(r = R.version.string, platform = R.version$platform),
   package_versions = deps
 )

diff --git a/analysis/01_qc_diagnostics.R b/analysis/01_qc_diagnostics.R
index 4fd33f3..c66d95d 100644
--- a/analysis/01_qc_diagnostics.R
+++ b/analysis/01_qc_diagnostics.R
@@ -12,7 +12,7 @@ run_qc <- function(context) {
   cfg <- context$config
   qc_dir <- cfg$output$dirs$qc
   assignments_map <- context$assignments
-
+
   if (is.null(assignments_map) || length(assignments_map) == 0) {
     return(list(
       status = "skipped",
@@ -20,13 +20,13 @@ run_qc <- function(context) {
       outputs = character(0)
     ))
   }
-
+
   dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
   reconciliation_rows <- list()
   length_summary_rows <- list()
-
+
   # Target lengths from params.json or config fallbacks
   min_len_target <- if (!is.null(context$params$min_len)) {
     context$params$min_len
@@ -38,39 +38,42 @@ run_qc <- function(context) {
   } else {
     cfg$qc$target_max_length
   }
-
+
   display_min <- cfg$qc$display_min_length %||% 1200L
   display_max <- cfg$qc$display_max_length %||% 1800L
-
+
   for (sample_id in context$samples) {
     asgn_path <- assignments_map[[sample_id]]
     if (is.null(asgn_path) || !file.exists(asgn_path)) {
       next
     }
-
+
     sample_out_dir <- file.path(qc_dir, sanitize_filename(sample_id))
     dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
-
+
     # Expected counts from abundance context
     stat_row <- context$sample_stats[context$sample_stats$SampleID == sample_id, ]
     exp_total <- stat_row$TotalReads[1]
     exp_class <- stat_row$ClassifiedReads[1]
     exp_unclass <- stat_row$UnclassifiedReads[1]
-
-    reads <- read_assignments_file(
-      path = asgn_path,
-      sample_id = sample_id,
-      expected_total = exp_total,
-      expected_classified = exp_class,
-      expected_unclassified = exp_unclass
-    )
-
+
+    reads <- context$assignment_data[[sample_id]]
+    if (is.null(reads)) {
+      reads <- read_assignments_file(
+        path = asgn_path,
+        sample_id = sample_id,
+        expected_total = exp_total,
+        expected_classified = exp_class,
+        expected_unclassified = exp_unclass
+      )
+    }
+
     n_status_C <- sum(reads$status == "C")
     n_status_U <- sum(reads$status == "U")
     n_qc_reclass <- sum(reads$status == "C" & reads$taxid == 0)
     n_eff_class <- sum(reads$effective_classified)
     n_eff_unclass <- sum(!reads$effective_classified)
-
+
     # 1. Reconciliation Table Row
     reconciliation_rows[[sample_id]] <- data.frame(
       SampleID = sample_id,
@@ -86,14 +89,14 @@ run_qc <- function(context) {
       ReconciliationPass = (n_eff_class == exp_class) && (nrow(reads) == exp_total),
       stringsAsFactors = FALSE
     )
-
+
     # 2. Length summary statistics
     reads$status_category <- ifelse(
       reads$status == "C" & reads$taxid != 0, "Classified",
       ifelse(reads$status == "C" & reads$taxid == 0, "QC-filtered", "Never aligned")
     )
     reads$effective_status <- ifelse(reads$effective_classified, "Classified", "Unclassified")
-
+
     for (cat_name in c("All", "Classified", "Unclassified", "QC-filtered")) {
       sub_lens <- if (cat_name == "All") {
         reads$read_length
@@ -102,7 +105,7 @@ run_qc <- function(context) {
       } else {
         reads$read_length[reads$status_category == cat_name]
       }
-
+
       if (length(sub_lens) > 0) {
         length_summary_rows[[paste(sample_id, cat_name, sep = "_")]] <- data.frame(
           SampleID = sample_id,
@@ -119,7 +122,7 @@ run_qc <- function(context) {
         )
       }
     }
-
+
     # 3. Figure 1a: Donut Plot
     donut_df <- reads %>%
       count(effective_status) %>%
@@ -129,7 +132,7 @@ run_qc <- function(context) {
         ymin = c(0, head(ymax, -1)),
         label = sprintf("%s\n%s reads\n(%.1f%%)", effective_status, scales::comma(n), 100 * frac)
       )
-
+
     p1a <- ggplot(donut_df, aes(ymin = ymin, ymax = ymax, xmin = 3, xmax = 4, fill = effective_status)) +
       geom_rect(color = "white", linewidth = 1.2) +
       coord_polar(theta = "y") +
@@ -142,11 +145,11 @@ run_qc <- function(context) {
       theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) +
       geom_text(aes(x = 3.5, y = (ymin + ymax) / 2, label = label), inherit.aes = FALSE,
                 data = donut_df, color = "black", size = 3.3, fontface = "bold")
-
+
     p1a_path <- file.path(sample_out_dir, "01a_classification_donut.png")
     save_plot(p1a_path, p1a, width = 5.5, height = 5.5)
     all_outputs <- c(all_outputs, p1a_path)
-
+
     # 4. Figure 1b: Read Length Histogram
     p1b <- ggplot(reads, aes(x = read_length, fill = effective_status)) +
       geom_histogram(binwidth = 10, alpha = 0.85, position = "identity") +
@@ -165,18 +168,18 @@ run_qc <- function(context) {
       ) +
       theme_amplicon() +
       theme(legend.position = "bottom")
-
+
     if (!is.null(min_len_target)) {
       p1b <- p1b + geom_vline(xintercept = min_len_target, linetype = "dotted", color = "grey30")
     }
     if (!is.null(max_len_target)) {
       p1b <- p1b + geom_vline(xintercept = max_len_target, linetype = "dotted", color = "grey30")
     }
-
+
     p1b_path <- file.path(sample_out_dir, "01b_read_length_distribution.png")
     save_plot(p1b_path, p1b, width = 7.5, height = 5)
     all_outputs <- c(all_outputs, p1b_path)
-
+
     # 5. Figure 1c: Diagnostic Bar (Why raw status != effective classification)
     qc_diag_df <- data.frame(
       Category = factor(
@@ -189,7 +192,7 @@ run_qc <- function(context) {
       ),
       Count = c(n_eff_class, n_qc_reclass, n_status_U)
     )
-
+
     p1c <- ggplot(qc_diag_df, aes(x = Category, y = Count, fill = Category)) +
       geom_col(width = 0.6) +
       geom_text(aes(label = scales::comma(Count)), vjust = -0.4, size = 3.6, fontface = "bold") +
@@ -209,11 +212,11 @@ run_qc <- function(context) {
       ) +
       theme_amplicon() +
       theme(legend.position = "none")
-
+
     p1c_path <- file.path(sample_out_dir, "01c_qc_filter_diagnostic.png")
     save_plot(p1c_path, p1c, width = 7, height = 5.5)
     all_outputs <- c(all_outputs, p1c_path)
-
+
     # 6. Figure 1d: Violin & Boxplot of Read Length by Effective Class
     p1d <- ggplot(reads, aes(x = effective_status, y = read_length, fill = effective_status)) +
       geom_violin(alpha = 0.6, trim = FALSE) +
@@ -228,12 +231,12 @@ run_qc <- function(context) {
       ) +
       theme_amplicon() +
       theme(legend.position = "none")
-
+
     p1d_path <- file.path(sample_out_dir, "01d_read_length_violin.png")
     save_plot(p1d_path, p1d, width = 6, height = 5)
     all_outputs <- c(all_outputs, p1d_path)
   }
-
+
   # Export summary TSVs
   reconciliation_file <- file.path(qc_dir, "classification_reconciliation.tsv")
   if (length(reconciliation_rows) > 0) {
@@ -241,14 +244,14 @@ run_qc <- function(context) {
                 sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, reconciliation_file)
   }
-
+
   length_file <- file.path(qc_dir, "read_length_by_status.tsv")
   if (length(length_summary_rows) > 0) {
     write.table(do.call(rbind, length_summary_rows), length_file,
                 sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, length_file)
   }
-
+
   list(
     status = "completed",
     outputs = all_outputs
diff --git a/analysis/02_alpha_diversity.R b/analysis/02_alpha_diversity.R
index 61818a3..cc32435 100644
--- a/analysis/02_alpha_diversity.R
+++ b/analysis/02_alpha_diversity.R
@@ -14,23 +14,23 @@ run_alpha <- function(context) {
   cfg <- context$config
   alpha_dir <- cfg$output$dirs$alpha
   dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
-
+
   # Parameters
   seed <- cfg$seed %||% 42L
   n_points <- cfg$alpha$rarefaction_points %||% 25L
   resample_depth_cfg <- cfg$alpha$resample_depth %||% 50000L
   fraction_cap <- cfg$alpha$resample_fraction_cap %||% 0.90
   n_iterations <- cfg$alpha$resample_iterations %||% 100L
-
+
   unclass_idx <- context$unclass_index
   count_matrix <- context$count_matrix
   # Extract classified counts only
   class_matrix <- count_matrix[-unclass_idx, , drop = FALSE]
-
+
   samples <- context$samples
-
+
   # 1. Compute Alpha Diversity Indices per sample
   alpha_records <- list()
   for (s in samples) {
@@ -40,20 +40,21 @@ run_alpha <- function(context) {
     alpha_records[[s]] <- idx_df
   }
   alpha_combined <- do.call(rbind, alpha_records)
-
+
   # Pivot to wide table for export: SampleID, Observed_Richness, Chao1, etc.
   alpha_wide <- alpha_combined %>%
     tidyr::pivot_wider(names_from = Metric, values_from = Value)
-
+  alpha_wide <- alpha_wide[match(samples, alpha_wide$SampleID), , drop = FALSE]
+
   # Attach metadata if available
   if (!is.null(context$metadata)) {
-    alpha_wide <- merge(context$metadata, alpha_wide, by = "SampleID", sort = FALSE)
+    alpha_wide <- dplyr::left_join(context$metadata, alpha_wide, by = "SampleID")
   }
-
+
   alpha_tsv <- file.path(alpha_dir, "alpha_diversity.tsv")
   write.table(alpha_wide, alpha_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, alpha_tsv)
-
+
   # 2. Analytical Rarefaction Curves
   rare_curve_records <- list()
   for (s in samples) {
@@ -65,11 +66,11 @@ run_alpha <- function(context) {
     }
   }
   rare_curve_combined <- do.call(rbind, rare_curve_records)
-
+
   rare_tsv <- file.path(alpha_dir, "rarefaction_curve.tsv")
   write.table(rare_curve_combined, rare_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, rare_tsv)
-
+
   # Figure 2a: Rarefaction Curve Plot
   p_curve <- ggplot(rare_curve_combined, aes(x = depth, y = mean_richness, group = SampleID, color = SampleID)) +
     geom_ribbon(aes(ymin = mean_richness - sd_richness, ymax = mean_richness + sd_richness, fill = SampleID),
@@ -85,7 +86,7 @@ run_alpha <- function(context) {
       color = "Sample", fill = "Sample"
     ) +
     theme_amplicon()
-
+
   # For single sample, add actual classified depth dashed line
   if (length(samples) == 1) {
     single_depth <- sum(class_matrix[, samples[1]])
@@ -95,22 +96,18 @@ run_alpha <- function(context) {
                label = "actual depth ", hjust = 1, vjust = 0, color = "grey40", size = 3) +
       theme(legend.position = "none")
   }
-
+
   p_curve_path <- file.path(alpha_dir, "02a_rarefaction_curve.png")
   save_plot(p_curve_path, p_curve, width = 8, height = 5.5)
   all_outputs <- c(all_outputs, p_curve_path)
-
+
   # 3. Rarefaction Resamples
   # Determine realized depth
   sample_depths <- colSums(class_matrix)
   min_depth <- min(sample_depths)
-
-  realized_depth <- if (context$mode == "single") {
-    min(resample_depth_cfg, floor(min_depth * fraction_cap))
-  } else {
-    floor(min_depth * fraction_cap)
-  }
-
+
+  realized_depth <- min(resample_depth_cfg, floor(min_depth * fraction_cap))
+
   resample_records <- list()
   for (s in samples) {
     counts_s <- class_matrix[, s]
@@ -119,19 +116,21 @@ run_alpha <- function(context) {
         counts = counts_s,
         subsample_depth = realized_depth,
         n_iterations = n_iterations,
-        seed = seed
+        seed = as.integer((seed + strtoi(substr(
+          digest::digest(s, algo = "xxhash32", serialize = FALSE), 1L, 7L
+        ), base = 16L)) %% .Machine$integer.max)
       )
       res_df$SampleID <- s
       resample_records[[s]] <- res_df
     }
   }
-
+
   if (length(resample_records) > 0) {
     resamples_combined <- do.call(rbind, resample_records)
     resample_tsv <- file.path(alpha_dir, "rarefaction_resamples.tsv")
     write.table(resamples_combined, resample_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, resample_tsv)
-
+
     # Figure 2b: Resample Boxplots (Sensitivity of alpha metrics across resamples)
     resamples_long <- resamples_combined %>%
       select(SampleID, richness, shannon, ens, simpson, invsimpson) %>%
@@ -142,7 +141,7 @@ run_alpha <- function(context) {
         levels = c("richness", "shannon", "ens", "simpson", "invsimpson"),
         labels = c("Richness (S)", "Shannon (H)", "ENS (e^H)", "Simpson (D)", "Inv. Simpson")
       ))
-
+
     p_box <- ggplot(resamples_long, aes(x = SampleID, y = Value, fill = SampleID)) +
       geom_boxplot(alpha = 0.7, outlier.size = 1) +
       facet_wrap(~Metric, scales = "free_y", nrow = 2) +
@@ -157,25 +156,25 @@ run_alpha <- function(context) {
         axis.text.x = if (length(samples) > 5) element_text(angle = 45, hjust = 1) else element_text(),
         legend.position = if (length(samples) == 1) "none" else "right"
       )
-
+
     p_box_path <- file.path(alpha_dir, "02b_resample_boxplots.png")
     save_plot(p_box_path, p_box, width = 9, height = 6)
     all_outputs <- c(all_outputs, p_box_path)
   }
-
+
   # 4. Cohort Group Tests (Only when mode == cohort)
   if (context$mode == "cohort") {
     meta <- context$metadata
-    alpha_with_group <- merge(meta, alpha_wide, by = "SampleID")
+    alpha_with_group <- alpha_wide
     groups <- unique(alpha_with_group$Group)
-
+
     # Gating: At least 2 groups with at least 3 biological samples each
     group_counts <- table(alpha_with_group$Group)
-    valid_groups <- names(group_counts)[group_counts >= 3]
-
+    all_groups_replicated <- length(group_counts) >= 2L && all(group_counts >= 3L)
+
     diff_file <- file.path(alpha_dir, "group_differences.tsv")
-
-    if (length(valid_groups) < 2) {
+
+    if (!all_groups_replicated) {
       skip_note <- data.frame(
         Status = "Skipped",
         Reason = sprintf(
@@ -190,15 +189,18 @@ run_alpha <- function(context) {
       metrics_to_test <- c("Observed species richness (S)", "Chao1 (estimated richness)",
                            "Shannon (H)", "Effective number of species (e^H)",
                            "Simpson's D (1-sum p^2)", "Inverse Simpson", "Pielou's evenness (J)")
-
+
       test_rows <- list()
       for (m in metrics_to_test) {
         if (m %in% colnames(alpha_with_group)) {
-          df_sub <- alpha_with_group[alpha_with_group$Group %in% valid_groups, ]
+          df_sub <- alpha_with_group[is.finite(alpha_with_group[[m]]), , drop = FALSE]
           y <- df_sub[[m]]
           grp <- factor(df_sub$Group)
-
-          if (length(valid_groups) == 2) {
+
+          metric_group_counts <- table(grp)
+          if (length(metric_group_counts) < 2L || any(metric_group_counts < 3L)) next
+
+          if (nlevels(grp) == 2L) {
             wt <- suppressWarnings(wilcox.test(y ~ grp))
             test_rows[[m]] <- data.frame(
               Metric = m, Test = "Wilcoxon rank-sum",
@@ -215,13 +217,22 @@ run_alpha <- function(context) {
           }
         }
       }
-      diff_df <- do.call(rbind, test_rows)
-      diff_df$FDR_BH <- p.adjust(diff_df$P_Value, method = "BH")
-      write.table(diff_df, diff_file, sep = "\t", row.names = FALSE, quote = FALSE)
+      if (length(test_rows) == 0L) {
+        skip_note <- data.frame(
+          Status = "Skipped",
+          Reason = "No metric retained at least 3 finite observations in every group.",
+          stringsAsFactors = FALSE
+        )
+        write.table(skip_note, diff_file, sep = "\t", row.names = FALSE, quote = FALSE)
+      } else {
+        diff_df <- do.call(rbind, test_rows)
+        diff_df$FDR_BH <- p.adjust(diff_df$P_Value, method = "BH")
+        write.table(diff_df, diff_file, sep = "\t", row.names = FALSE, quote = FALSE)
+      }
     }
     all_outputs <- c(all_outputs, diff_file)
   }
-
+
   list(
     status = "completed",
     outputs = all_outputs
diff --git a/analysis/03_beta_diversity.R b/analysis/03_beta_diversity.R
index cf31284..137396c 100644
--- a/analysis/03_beta_diversity.R
+++ b/analysis/03_beta_diversity.R
@@ -12,9 +12,9 @@ run_beta <- function(context) {
   cfg <- context$config
   beta_dir <- cfg$output$dirs$beta
   dir.create(beta_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
-
+
   # Cohort gate
   if (context$mode != "cohort" || length(context$samples) < 2) {
     skip_file <- file.path(beta_dir, "beta_diversity_skipped.tsv")
@@ -31,59 +31,87 @@ run_beta <- function(context) {
       outputs = skip_file
     ))
   }
-
+
   samples <- context$samples
   unclass_idx <- context$unclass_index
   count_matrix <- context$count_matrix
   class_counts <- count_matrix[-unclass_idx, , drop = FALSE]
-
+
   meta <- context$metadata
   seed <- cfg$seed %||% 42L
   set.seed(seed)
-
+
   # Samples as rows, taxa as columns
   otu_table <- t(class_counts)
   # Calculate classified relative abundances for Bray-Curtis
   sample_sums <- rowSums(otu_table)
   rel_otu <- sweep(otu_table, 1, sample_sums, "/")
-
+
   distances_cfg <- cfg$beta$distances %||% c("bray", "jaccard")
   n_perm <- cfg$beta$permutations %||% 999L
+  minimum_count <- cfg$beta$minimum_count %||% 1L
   strata_col <- cfg$beta$strata_column
-
+
   dist_list <- list()
-
+
   for (d_name in distances_cfg) {
     dist_mat <- if (d_name == "bray") {
       vegan::vegdist(rel_otu, method = "bray")
     } else if (d_name == "jaccard") {
-      vegan::vegdist(otu_table > 0, method = "jaccard", binary = TRUE)
+      vegan::vegdist(otu_table >= minimum_count, method = "jaccard", binary = TRUE)
     } else {
-      vegan::vegdist(rel_otu, method = d_name)
+      stop(sprintf("Unsupported beta-diversity distance: '%s'", d_name), call. = FALSE)
     }
     dist_list[[d_name]] <- dist_mat
-
+
     # Save distance matrix TSV
     d_tsv <- file.path(beta_dir, sprintf("distance_%s.tsv", d_name))
-    write.table(as.matrix(dist_mat), d_tsv, sep = "\t", quote = FALSE)
+    distance_df <- data.frame(
+      SampleID = rownames(as.matrix(dist_mat)),
+      as.matrix(dist_mat),
+      check.names = FALSE,
+      stringsAsFactors = FALSE
+    )
+    write.table(distance_df, d_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, d_tsv)
-
-    # PCoA via cmdscale
+
+    distance_values <- as.vector(dist_mat)
+    if (!any(is.finite(distance_values) & distance_values > 0)) {
+      skip_path <- file.path(beta_dir, sprintf("pcoa_skipped_%s.tsv", d_name))
+      write.table(data.frame(
+        Status = "Skipped",
+        Reason = "All pairwise distances are zero; PCoA is undefined."
+      ), skip_path, sep = "\t", row.names = FALSE, quote = FALSE)
+      all_outputs <- c(all_outputs, skip_path)
+      next
+    }
+
+    # PCoA via additive-corrected classical scaling.
     max_k <- max(1, min(nrow(otu_table) - 1, 2))
-    pcoa <- stats::cmdscale(dist_mat, k = max_k, eig = TRUE)
+    pcoa <- tryCatch(
+      stats::cmdscale(dist_mat, k = max_k, eig = TRUE, add = TRUE),
+      error = function(e) structure(list(message = conditionMessage(e)), class = "pcoa_error")
+    )
+    if (inherits(pcoa, "pcoa_error")) {
+      skip_path <- file.path(beta_dir, sprintf("pcoa_skipped_%s.tsv", d_name))
+      write.table(data.frame(Status = "Skipped", Reason = pcoa$message), skip_path,
+                  sep = "\t", row.names = FALSE, quote = FALSE)
+      all_outputs <- c(all_outputs, skip_path)
+      next
+    }
     eig <- pcoa$eig
     pos_eig <- eig[eig > 0]
     total_pos <- if (length(pos_eig) > 0) sum(pos_eig) else 1
-
+
     var_exp <- c(
       if (length(eig) >= 1 && eig[1] > 0) round(100 * eig[1] / total_pos, 1) else 0,
       if (length(eig) >= 2 && eig[2] > 0) round(100 * eig[2] / total_pos, 1) else 0
     )
-
+
     pts <- as.matrix(pcoa$points)
     p1 <- if (ncol(pts) >= 1) pts[, 1] else rep(0, nrow(otu_table))
     p2 <- if (ncol(pts) >= 2) pts[, 2] else rep(0, nrow(otu_table))
-
+
     scores_df <- data.frame(
       SampleID = rownames(otu_table),
       PCoA1 = p1,
@@ -91,15 +119,26 @@ run_beta <- function(context) {
       stringsAsFactors = FALSE
     )
     if (!is.null(meta)) {
-      scores_df <- merge(meta, scores_df, by = "SampleID", sort = FALSE)
+      scores_df <- dplyr::left_join(meta, scores_df, by = "SampleID")
     }
-
+
     scores_tsv <- file.path(beta_dir, sprintf("pcoa_scores_%s.tsv", d_name))
     write.table(scores_df, scores_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, scores_tsv)
-
+
+    variance_tsv <- file.path(beta_dir, sprintf("pcoa_variance_%s.tsv", d_name))
+    variance_df <- data.frame(
+      Axis = paste0("PCoA", seq_along(eig)),
+      Eigenvalue = eig,
+      PositiveVariancePercent = ifelse(eig > 0, 100 * eig / total_pos, 0),
+      stringsAsFactors = FALSE
+    )
+    write.table(variance_df, variance_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
+    all_outputs <- c(all_outputs, variance_tsv)
+
     # PCoA 2D plot (only if at least 3 samples exist)
-    if (length(samples) >= 3 && ncol(pcoa$points) >= 2) {
+    n_unique_samples <- nrow(unique(as.data.frame(rel_otu)))
+    if (length(samples) >= 3 && n_unique_samples >= 3 && ncol(pts) >= 2) {
       p_pcoa <- ggplot(scores_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
         geom_point(size = 3.5, alpha = 0.85) +
         labs(
@@ -109,27 +148,94 @@ run_beta <- function(context) {
           y = sprintf("PCoA 2 (%.1f%%)", var_exp[2])
         ) +
         theme_amplicon()
-
+
       p_path <- file.path(beta_dir, sprintf("03_pcoa_%s.png", d_name))
       save_plot(p_path, p_pcoa, width = 7, height = 5.5)
       all_outputs <- c(all_outputs, p_path)
     }
   }
-
+
+  # Optional rarefaction stability analysis. Full-data distances/PCoA above
+  # remain primary; these coordinates are Procrustes-aligned sensitivity draws.
+  if (isTRUE(cfg$beta$resampling$enabled)) {
+    stability_file <- file.path(beta_dir, "pcoa_rarefaction_stability.tsv")
+    stability_diag <- file.path(beta_dir, "pcoa_rarefaction_diagnostics.tsv")
+    stability_depth <- floor(min(sample_sums) * cfg$beta$resampling$depth_fraction_of_minimum)
+    reference_dist <- vegan::vegdist(rel_otu, method = "bray")
+    reference_fit <- tryCatch(
+      stats::cmdscale(reference_dist, k = 2L, eig = TRUE, add = TRUE),
+      error = function(e) NULL
+    )
+
+    if (nrow(otu_table) < 3L || stability_depth < 1L || is.null(reference_fit) ||
+        ncol(as.matrix(reference_fit$points)) < 2L) {
+      write.table(data.frame(
+        Status = "Skipped",
+        Reason = "Rarefaction stability requires >=3 samples and a two-axis full-data PCoA solution."
+      ), stability_diag, sep = "\t", row.names = FALSE, quote = FALSE)
+      all_outputs <- c(all_outputs, stability_diag)
+    } else {
+      set.seed(seed)
+      stability_rows <- list()
+      failed_iterations <- integer(0)
+      reference_points <- as.matrix(reference_fit$points)[, 1:2, drop = FALSE]
+      for (iteration in seq_len(cfg$beta$resampling$iterations)) {
+        rare_counts <- vegan::rrarefy(otu_table, sample = stability_depth)
+        rare_rel <- sweep(rare_counts, 1, rowSums(rare_counts), "/")
+        rare_dist <- vegan::vegdist(rare_rel, method = "bray")
+        rare_fit <- tryCatch(
+          stats::cmdscale(rare_dist, k = 2L, eig = TRUE, add = TRUE),
+          error = function(e) NULL
+        )
+        if (is.null(rare_fit) || ncol(as.matrix(rare_fit$points)) < 2L) {
+          failed_iterations <- c(failed_iterations, iteration)
+          next
+        }
+        aligned <- vegan::procrustes(reference_points,
+                                     as.matrix(rare_fit$points)[, 1:2, drop = FALSE])$Yrot
+        stability_rows[[length(stability_rows) + 1L]] <- data.frame(
+          Iteration = iteration,
+          SampleID = rownames(aligned),
+          PCoA1 = aligned[, 1],
+          PCoA2 = aligned[, 2],
+          Depth = stability_depth,
+          stringsAsFactors = FALSE
+        )
+      }
+      if (length(stability_rows) == 0L) {
+        stop("All beta-diversity rarefaction stability iterations failed.", call. = FALSE)
+      }
+      write.table(do.call(rbind, stability_rows), stability_file,
+                  sep = "\t", row.names = FALSE, quote = FALSE)
+      write.table(data.frame(
+        Status = "Completed",
+        RequestedIterations = cfg$beta$resampling$iterations,
+        SuccessfulIterations = length(stability_rows),
+        FailedIterations = paste(failed_iterations, collapse = ","),
+        Depth = stability_depth,
+        Seed = seed
+      ), stability_diag, sep = "\t", row.names = FALSE, quote = FALSE)
+      all_outputs <- c(all_outputs, stability_file, stability_diag)
+    }
+  }
+
   # PERMANOVA & Betadisper (Primary distance: Bray-Curtis)
   primary_dist <- dist_list[["bray"]] %||% dist_list[[1]]
-
+
   permanova_file <- file.path(beta_dir, "permanova.tsv")
   betadisper_file <- file.path(beta_dir, "betadisper.tsv")
-
+
   # Gating: At least 2 groups with at least 2 samples per group
   group_counts <- table(meta$Group)
-  can_run_permanova <- length(group_counts) >= 2 && all(group_counts >= 2)
-
+  primary_values <- as.vector(primary_dist)
+  residual_df <- nrow(meta) - length(group_counts)
+  can_run_permanova <- length(group_counts) >= 2 && all(group_counts >= 2) &&
+    residual_df > 0 && any(is.finite(primary_values) & primary_values > 0)
+
   if (can_run_permanova) {
     # Check strata
-    strata_vec <- if (!is.null(strata_col) && strata_col %in% colnames(meta)) meta[[strata_col]] else NULL
-
+    strata_vec <- if (!is.null(strata_col)) meta[[strata_col]] else NULL
+
     set.seed(seed)
     perm_res <- suppressWarnings(vegan::adonis2(
       primary_dist ~ Group,
@@ -137,18 +243,18 @@ run_beta <- function(context) {
       permutations = n_perm,
       strata = strata_vec
     ))
-
+
     perm_df <- as.data.frame(perm_res)
     perm_df$Term <- rownames(perm_df)
     perm_df <- perm_df[, c("Term", setdiff(colnames(perm_df), "Term"))]
-
+
     write.table(perm_df, permanova_file, sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, permanova_file)
-
+
     # Betadisper
     disp_res <- vegan::betadisper(primary_dist, meta$Group)
     disp_perm <- vegan::permutest(disp_res, permutations = n_perm)
-
+
     disp_df <- data.frame(
       Analysis = "Betadisper (Homogeneity of Multivariate Dispersions)",
       F_Statistic = disp_perm$tab$F[1],
@@ -163,7 +269,7 @@ run_beta <- function(context) {
     skip_perm <- data.frame(
       Status = "Skipped",
       Reason = sprintf(
-        "PERMANOVA requires at least 2 groups with >= 2 samples each. Found: %s",
+        "PERMANOVA requires non-zero distances, residual degrees of freedom, and at least 2 groups with >= 2 samples each. Found: %s",
         paste(sprintf("%s (n=%d)", names(group_counts), as.integer(group_counts)), collapse = ", ")
       ),
       stringsAsFactors = FALSE
@@ -172,7 +278,7 @@ run_beta <- function(context) {
     write.table(skip_perm, betadisper_file, sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, permanova_file, betadisper_file)
   }
-
+
   list(
     status = "completed",
     outputs = all_outputs
diff --git a/analysis/04_taxa_composition.R b/analysis/04_taxa_composition.R
index e6196ec..a76d206 100644
--- a/analysis/04_taxa_composition.R
+++ b/analysis/04_taxa_composition.R
@@ -17,44 +17,57 @@ run_taxa_composition <- function(context) {
   cfg <- context$config
   comp_dir <- cfg$output$dirs$composition
   dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
-
+
   top_n <- cfg$composition$top_n_taxa %||% 15L
   heatmap_rank <- cfg$composition$heatmap_rank %||% "genus"
   heatmap_transform <- cfg$composition$heatmap_transform %||% "log10_relative"
-
+
   samples <- context$samples
   unclass_idx <- context$unclass_index
   count_matrix <- context$count_matrix
   taxonomy_df <- context$taxonomy
-
+
   # Total and classified read denominators per sample
   sample_totals <- colSums(count_matrix)
   unclass_counts <- count_matrix[unclass_idx, ]
   class_totals <- sample_totals - unclass_counts
-
+
   # Classified rows
   class_counts <- count_matrix[-unclass_idx, , drop = FALSE]
   class_tax <- taxonomy_df[-unclass_idx, , drop = FALSE]
-
+
+  # All-read classification fraction is available even without per-read assignments.
+  classification_df <- data.frame(
+    SampleID = samples,
+    TotalReads = as.numeric(sample_totals[samples]),
+    ClassifiedReads = as.numeric(class_totals[samples]),
+    UnclassifiedReads = as.numeric(unclass_counts[samples]),
+    ClassifiedFraction = as.numeric(class_totals[samples] / sample_totals[samples]),
+    stringsAsFactors = FALSE
+  )
+  classification_file <- file.path(comp_dir, "classification_fraction.tsv")
+  write.table(classification_df, classification_file, sep = "\t", row.names = FALSE, quote = FALSE)
+  all_outputs <- c(all_outputs, classification_file)
+
   # Analyze each rank
   rank_tables <- list()
-
+
   for (rk in RANKS_TO_ANALYZE) {
     rk_idx <- which(colnames(taxonomy_df) == rk)
-
+
     # Prefix-based aggregation up to rank rk
     rk_prefixes <- apply(class_tax[, 1:rk_idx, drop = FALSE], 1, paste, collapse = ";")
     leaf_names <- class_tax[[rk]]
-
+
     # Contextualize ambiguous labels (e.g. "Unknown" gets parent context)
     contextualized_names <- ifelse(
       leaf_names %in% c("Unknown", "unclassified", "uncultured"),
       sprintf("%s (%s)", leaf_names, class_tax[[rk_idx - 1]]),
       leaf_names
     )
-
+
     # Aggregate counts by prefix
     agg_df <- data.frame(
       TaxonPath = rk_prefixes,
@@ -65,35 +78,42 @@ run_taxa_composition <- function(context) {
     ) %>%
       group_by(TaxonPath, Taxon) %>%
       summarise(across(all_of(samples), sum), .groups = "drop")
-
+
+    # Display labels must remain unique even when the same name occurs under
+    # different parents. TaxonPath remains the analytical key.
+    duplicated_label <- duplicated(agg_df$Taxon) | duplicated(agg_df$Taxon, fromLast = TRUE)
+    agg_df$Taxon[duplicated_label] <- sprintf(
+      "%s [%s]", agg_df$Taxon[duplicated_label], agg_df$TaxonPath[duplicated_label]
+    )
+
     # Calculate relative abundances (classified-only denominator)
     rel_df <- agg_df
     for (s in samples) {
       denom <- class_totals[s]
       rel_df[[s]] <- if (denom > 0) rel_df[[s]] / denom else 0
     }
-
+
     # Write TSVs
     count_file <- file.path(comp_dir, sprintf("count_%s.tsv", rk))
     rel_file <- file.path(comp_dir, sprintf("rel_abundance_%s.tsv", rk))
-
+
     write.table(agg_df, count_file, sep = "\t", row.names = FALSE, quote = FALSE)
     write.table(rel_df, rel_file, sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, count_file, rel_file)
-
+
     rank_tables[[rk]] <- list(counts = agg_df, rel = rel_df)
   }
-
+
   # Single-Sample Bar Plots
   if (length(samples) == 1) {
     s_col <- samples[1]
-
+
     # 1. Phylum bar plot
     phylum_rel <- rank_tables[["phylum"]]$rel %>%
       select(Taxon, all_of(s_col)) %>%
       rename(rel = all_of(s_col)) %>%
       arrange(desc(rel))
-
+
     n_keep_phylum <- min(7L, nrow(phylum_rel))
     phylum_top <- phylum_rel %>%
       mutate(
@@ -102,9 +122,9 @@ run_taxa_composition <- function(context) {
       group_by(DisplayTaxon) %>%
       summarise(rel = sum(rel), .groups = "drop") %>%
       arrange(desc(rel))
-
+
     phylum_top$DisplayTaxon <- factor(phylum_top$DisplayTaxon, levels = rev(phylum_top$DisplayTaxon))
-
+
     p_phylum <- ggplot(phylum_top, aes(x = DisplayTaxon, y = rel, fill = DisplayTaxon)) +
       geom_col(width = 0.7) +
       geom_text(aes(label = sprintf("%.1f%%", 100 * rel)), hjust = -0.15, size = 3.3) +
@@ -118,11 +138,11 @@ run_taxa_composition <- function(context) {
       ) +
       theme_amplicon() +
       theme(legend.position = "none")
-
+
     p_phylum_path <- file.path(comp_dir, "04a_phylum_composition.png")
     save_plot(p_phylum_path, p_phylum, width = 7.5, height = 5)
     all_outputs <- c(all_outputs, p_phylum_path)
-
+
     # 2. Top Family bar plot
     family_counts <- rank_tables[["family"]]$counts %>%
       select(Taxon, all_of(s_col)) %>%
@@ -130,7 +150,7 @@ run_taxa_composition <- function(context) {
       arrange(desc(count)) %>%
       slice_head(n = top_n)
     family_counts$Taxon <- factor(family_counts$Taxon, levels = rev(family_counts$Taxon))
-
+
     p_family <- ggplot(family_counts, aes(x = Taxon, y = count)) +
       geom_col(width = 0.7, fill = "#1b9e77") +
       geom_text(aes(label = scales::comma(count)), hjust = -0.15, size = 3) +
@@ -141,11 +161,11 @@ run_taxa_composition <- function(context) {
         x = NULL, y = "Read count"
       ) +
       theme_amplicon()
-
+
     p_family_path <- file.path(comp_dir, "04b_family_composition.png")
     save_plot(p_family_path, p_family, width = 8.5, height = 6)
     all_outputs <- c(all_outputs, p_family_path)
-
+
     # 3. Top Genus bar plot
     genus_counts <- rank_tables[["genus"]]$counts %>%
       select(Taxon, all_of(s_col)) %>%
@@ -153,7 +173,7 @@ run_taxa_composition <- function(context) {
       arrange(desc(count)) %>%
       slice_head(n = top_n)
     genus_counts$Taxon <- factor(genus_counts$Taxon, levels = rev(genus_counts$Taxon))
-
+
     p_genus <- ggplot(genus_counts, aes(x = Taxon, y = count)) +
       geom_col(width = 0.7, fill = "#d95f02") +
       geom_text(aes(label = scales::comma(count)), hjust = -0.15, size = 3) +
@@ -164,11 +184,11 @@ run_taxa_composition <- function(context) {
         x = NULL, y = "Read count"
       ) +
       theme_amplicon()
-
+
     p_genus_path <- file.path(comp_dir, "04c_genus_composition.png")
     save_plot(p_genus_path, p_genus, width = 8.5, height = 6)
     all_outputs <- c(all_outputs, p_genus_path)
-
+
     # 4. Top Species bar plot
     species_rel <- rank_tables[["species"]]$rel %>%
       select(Taxon, all_of(s_col)) %>%
@@ -176,7 +196,7 @@ run_taxa_composition <- function(context) {
       arrange(desc(rel)) %>%
       slice_head(n = top_n)
     species_rel$Taxon <- factor(species_rel$Taxon, levels = rev(species_rel$Taxon))
-
+
     p_species <- ggplot(species_rel, aes(x = Taxon, y = rel)) +
       geom_col(width = 0.7, fill = "#7570b3") +
       geom_text(aes(label = sprintf("%.2f%%", 100 * rel)), hjust = -0.1, size = 3) +
@@ -188,12 +208,12 @@ run_taxa_composition <- function(context) {
       ) +
       theme_amplicon() +
       theme(axis.text.y = element_text(face = "italic"))
-
+
     p_species_path <- file.path(comp_dir, "04d_species_composition.png")
     save_plot(p_species_path, p_species, width = 9, height = 6)
     all_outputs <- c(all_outputs, p_species_path)
   }
-
+
   # Cohort Stacked Bars & Heatmap
   if (length(samples) >= 2) {
     # 1. Phylum Stacked Bar
@@ -202,21 +222,21 @@ run_taxa_composition <- function(context) {
       group_by(Taxon) %>%
       mutate(mean_rel = mean(rel)) %>%
       ungroup()
-
+
     top_phyla <- phylum_long %>%
       distinct(Taxon, mean_rel) %>%
       arrange(desc(mean_rel)) %>%
       slice_head(n = 8) %>%
       pull(Taxon)
-
+
     phylum_stacked <- phylum_long %>%
       mutate(DisplayTaxon = if_else(Taxon %in% top_phyla, Taxon, "Other")) %>%
       group_by(SampleID, DisplayTaxon) %>%
       summarise(rel = sum(rel), .groups = "drop")
-
+
     phylum_levels <- c(setdiff(unique(phylum_stacked$DisplayTaxon), "Other"), "Other")
     phylum_stacked$DisplayTaxon <- factor(phylum_stacked$DisplayTaxon, levels = rev(phylum_levels))
-
+
     p_stack <- ggplot(phylum_stacked, aes(x = SampleID, y = rel, fill = DisplayTaxon)) +
       geom_col(width = 0.6) +
       scale_fill_manual(values = get_phylum_colors(levels(phylum_stacked$DisplayTaxon)), name = "Phylum") +
@@ -227,39 +247,39 @@ run_taxa_composition <- function(context) {
       ) +
       theme_amplicon() +
       theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
-
+
     p_stack_path <- file.path(comp_dir, "04_phylum_stacked.png")
     save_plot(p_stack_path, p_stack, width = 8.5, height = 6.5)
     all_outputs <- c(all_outputs, p_stack_path)
-
+
     # 2. Heatmap
     if (heatmap_rank %in% names(rank_tables)) {
       rk_data <- rank_tables[[heatmap_rank]]$rel
-
+
       # Select top N taxa by mean relative abundance
       mat_data <- as.matrix(rk_data[, samples, drop = FALSE])
       rownames(mat_data) <- rk_data$Taxon
       mean_rel <- rowMeans(mat_data)
-
+
       top_idx <- order(mean_rel, decreasing = TRUE)[seq_len(min(top_n, nrow(mat_data)))]
       mat_top <- mat_data[top_idx, , drop = FALSE]
-
+
       # Apply transform
       mat_transformed <- if (heatmap_transform == "log10_relative") {
         log10(mat_top + 1e-4)
       } else {
         mat_top
       }
-
+
       # Annotate columns with metadata if available
       anno_col <- NA
       if (!is.null(context$metadata) && "Group" %in% colnames(context$metadata)) {
         anno_df <- data.frame(Group = context$metadata$Group, row.names = context$metadata$SampleID)
         anno_col <- anno_df[samples, , drop = FALSE]
       }
-
+
       heatmap_path <- file.path(comp_dir, sprintf("04_heatmap_%s.png", heatmap_rank))
-
+
       # Save pheatmap
       png(heatmap_path, width = 8, height = 7, units = "in", res = 150)
       pheatmap::pheatmap(
@@ -275,7 +295,24 @@ run_taxa_composition <- function(context) {
       all_outputs <- c(all_outputs, heatmap_path)
     }
   }
-
+
+  classification_long <- classification_df %>%
+    select(SampleID, ClassifiedReads, UnclassifiedReads) %>%
+    tidyr::pivot_longer(-SampleID, names_to = "Status", values_to = "Reads") %>%
+    mutate(Status = recode(Status,
+                           ClassifiedReads = "Classified",
+                           UnclassifiedReads = "Unclassified"))
+  p_class <- ggplot(classification_long, aes(x = SampleID, y = Reads, fill = Status)) +
+    geom_col(position = "fill") +
+    scale_y_continuous(labels = scales::percent) +
+    scale_fill_manual(values = c(Classified = "#1b9e77", Unclassified = "#bdbdbd")) +
+    labs(title = "All-read Classification Fraction", x = NULL, y = "Fraction of all reads") +
+    theme_amplicon() +
+    theme(axis.text.x = element_text(angle = if (length(samples) > 5) 45 else 0, hjust = 1))
+  classification_plot <- file.path(comp_dir, "04_classification_fraction.png")
+  save_plot(classification_plot, p_class, width = max(7, min(14, 0.45 * length(samples) + 5)), height = 5)
+  all_outputs <- c(all_outputs, classification_plot)
+
   list(
     status = "completed",
     outputs = all_outputs
diff --git a/analysis/05_ordination.R b/analysis/05_ordination.R
index e4987c5..a6f603e 100644
--- a/analysis/05_ordination.R
+++ b/analysis/05_ordination.R
@@ -12,9 +12,9 @@ run_ordination <- function(context) {
   cfg <- context$config
   ord_dir <- cfg$output$dirs$ordination
   dir.create(ord_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
-
+
   # Cohort gate
   if (context$mode != "cohort" || length(context$samples) < 2) {
     skip_file <- file.path(ord_dir, "ordination_skipped.tsv")
@@ -31,65 +31,76 @@ run_ordination <- function(context) {
       outputs = skip_file
     ))
   }
-
+
   samples <- context$samples
   unclass_idx <- context$unclass_index
   count_matrix <- context$count_matrix
   class_counts <- count_matrix[-unclass_idx, , drop = FALSE]
-
+
   meta <- context$metadata
   seed <- cfg$seed %||% 42L
   set.seed(seed)
-
+
   otu_table <- t(class_counts)
   sample_sums <- rowSums(otu_table)
   rel_otu <- sweep(otu_table, 1, sample_sums, "/")
-
+
   # 1. Hellinger PCA
   hel_otu <- vegan::decostand(rel_otu, method = "hellinger")
-  pca_res <- vegan::rda(hel_otu)
-
-  eig <- pca_res$CA$eig
-  var_exp <- round(100 * eig / sum(eig), 2)
-  var_df <- data.frame(
-    PC = paste0("PC", seq_along(var_exp)),
-    Eigenvalue = as.numeric(eig),
-    VarianceExplained = as.numeric(var_exp),
-    CumulativeVariance = cumsum(as.numeric(var_exp)),
-    stringsAsFactors = FALSE
-  )
-  var_tsv <- file.path(ord_dir, "pca_variance.tsv")
-  write.table(var_df, var_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
-  all_outputs <- c(all_outputs, var_tsv)
-
-  scores_mat <- as.matrix(scores(pca_res, display = "sites"))
-  scores_df <- data.frame(
-    SampleID = rownames(scores_mat),
-    PC1 = if (ncol(scores_mat) >= 1) scores_mat[, 1] else rep(0, nrow(scores_mat)),
-    PC2 = if (ncol(scores_mat) >= 2) scores_mat[, 2] else rep(0, nrow(scores_mat)),
-    stringsAsFactors = FALSE
-  )
-  if (!is.null(meta)) {
-    scores_df <- merge(meta, scores_df, by = "SampleID", sort = FALSE)
-  }
-  scores_tsv <- file.path(ord_dir, "pca_scores.tsv")
-  write.table(scores_df, scores_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
-  all_outputs <- c(all_outputs, scores_tsv)
-
-  # Species loadings
-  loadings_mat <- as.matrix(scores(pca_res, display = "species"))
-  loadings_df <- data.frame(
-    Taxon = rownames(loadings_mat),
-    PC1 = if (ncol(loadings_mat) >= 1) loadings_mat[, 1] else rep(0, nrow(loadings_mat)),
-    PC2 = if (ncol(loadings_mat) >= 2) loadings_mat[, 2] else rep(0, nrow(loadings_mat)),
-    stringsAsFactors = FALSE
-  )
-  load_tsv <- file.path(ord_dir, "pca_loadings.tsv")
-  write.table(loadings_df, load_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
-  all_outputs <- c(all_outputs, load_tsv)
-
-  # PCA plot (when at least 2 PCs available)
-  if (ncol(scores_mat) >= 2) {
+  centered_hel <- scale(hel_otu, center = TRUE, scale = FALSE)
+  matrix_rank <- qr(centered_hel)$rank
+  max_axes <- min(2L, nrow(hel_otu) - 1L, matrix_rank)
+
+  if (max_axes < 1L) {
+    pca_skip <- file.path(ord_dir, "pca_skipped.tsv")
+    write.table(data.frame(
+      Status = "Skipped",
+      Reason = "Hellinger profiles have zero between-sample variation; PCA is undefined."
+    ), pca_skip, sep = "\t", row.names = FALSE, quote = FALSE)
+    all_outputs <- c(all_outputs, pca_skip)
+  } else {
+    pca_res <- stats::prcomp(hel_otu, center = TRUE, scale. = FALSE, rank. = max_axes)
+    eig <- pca_res$sdev^2
+    var_exp <- round(100 * eig / sum(eig), 2)
+    var_df <- data.frame(
+      PC = paste0("PC", seq_along(var_exp)),
+      Eigenvalue = as.numeric(eig),
+      VarianceExplained = as.numeric(var_exp),
+      CumulativeVariance = cumsum(as.numeric(var_exp)),
+      stringsAsFactors = FALSE
+    )
+    var_tsv <- file.path(ord_dir, "pca_variance.tsv")
+    write.table(var_df, var_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
+    all_outputs <- c(all_outputs, var_tsv)
+
+    scores_mat <- as.matrix(pca_res$x[, seq_len(max_axes), drop = FALSE])
+    scores_df <- data.frame(
+      SampleID = rownames(scores_mat),
+      PC1 = scores_mat[, 1],
+      PC2 = if (max_axes >= 2L) scores_mat[, 2] else NA_real_,
+      stringsAsFactors = FALSE
+    )
+    if (!is.null(meta)) {
+      scores_df <- dplyr::left_join(meta, scores_df, by = "SampleID")
+    }
+    scores_tsv <- file.path(ord_dir, "pca_scores.tsv")
+    write.table(scores_df, scores_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
+    all_outputs <- c(all_outputs, scores_tsv)
+
+    # Species loadings
+    loadings_mat <- as.matrix(pca_res$rotation[, seq_len(max_axes), drop = FALSE])
+    loadings_df <- data.frame(
+      TaxonPath = rownames(loadings_mat),
+      PC1 = loadings_mat[, 1],
+      PC2 = if (max_axes >= 2L) loadings_mat[, 2] else NA_real_,
+      stringsAsFactors = FALSE
+    )
+    load_tsv <- file.path(ord_dir, "pca_loadings.tsv")
+    write.table(loadings_df, load_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
+    all_outputs <- c(all_outputs, load_tsv)
+
+    # PCA plot (when at least 2 PCs are identifiable)
+    if (max_axes >= 2L) {
     p_pca <- ggplot(scores_df, aes(x = PC1, y = PC2, color = Group)) +
       geom_point(size = 3.5, alpha = 0.85) +
       labs(
@@ -98,36 +109,39 @@ run_ordination <- function(context) {
         y = sprintf("PC2 (%.1f%%)", var_exp[2])
       ) +
       theme_amplicon()
-
+
     pca_plot_path <- file.path(ord_dir, "05a_pca_plot.png")
     save_plot(pca_plot_path, p_pca, width = 7, height = 5.5)
-    all_outputs <- c(all_outputs, pca_plot_path)
+      all_outputs <- c(all_outputs, pca_plot_path)
+    }
   }
-
+
   # 2. Bray-Curtis NMDS (requires at least 3 non-identical samples)
-  if (length(samples) >= 3) {
+  n_unique_samples <- nrow(unique(as.data.frame(rel_otu)))
+  nmds_diag_path <- file.path(ord_dir, "nmds_diagnostics.tsv")
+  if (length(samples) >= 3 && n_unique_samples >= 3) {
     set.seed(seed)
     nmds_res <- suppressWarnings(tryCatch({
       vegan::metaMDS(rel_otu, distance = "bray", k = 2, trymax = 50, trace = 0)
     }, error = function(e) NULL))
-
+
     if (!is.null(nmds_res) && !is.null(nmds_res$points)) {
       nmds_df <- data.frame(
         SampleID = rownames(nmds_res$points),
         NMDS1 = nmds_res$points[, 1],
         NMDS2 = nmds_res$points[, 2],
         Stress = nmds_res$stress,
-        Converged = nmds_res$converged,
+        Converged = isTRUE(nmds_res$converged),
         stringsAsFactors = FALSE
       )
       if (!is.null(meta)) {
-        nmds_df <- merge(meta, nmds_df, by = "SampleID", sort = FALSE)
+        nmds_df <- dplyr::left_join(meta, nmds_df, by = "SampleID")
       }
-
+
       nmds_tsv <- file.path(ord_dir, "nmds_scores.tsv")
       write.table(nmds_df, nmds_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
       all_outputs <- c(all_outputs, nmds_tsv)
-
+
       p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group)) +
         geom_point(size = 3.5, alpha = 0.85) +
         labs(
@@ -139,13 +153,34 @@ run_ordination <- function(context) {
           x = "NMDS1", y = "NMDS2"
         ) +
         theme_amplicon()
-
+
       nmds_plot_path <- file.path(ord_dir, "05b_nmds_plot.png")
       save_plot(nmds_plot_path, p_nmds, width = 7, height = 5.5)
       all_outputs <- c(all_outputs, nmds_plot_path)
+
+      write.table(data.frame(
+        Status = "Completed",
+        Stress = nmds_res$stress,
+        Converged = isTRUE(nmds_res$converged),
+        Tries = nmds_res$tries %||% NA_integer_,
+        Trymax = 50L,
+        Warning = if (nmds_res$stress > 0.2) "High stress (>0.2); interpret cautiously" else "None"
+      ), nmds_diag_path, sep = "\t", row.names = FALSE, quote = FALSE)
+    } else {
+      write.table(data.frame(
+        Status = "Skipped",
+        Reason = "metaMDS failed to return a valid two-dimensional solution."
+      ), nmds_diag_path, sep = "\t", row.names = FALSE, quote = FALSE)
     }
+  } else {
+    write.table(data.frame(
+      Status = "Skipped",
+      Reason = sprintf("NMDS requires at least 3 non-identical samples; found %d sample(s), %d unique profile(s).",
+                       length(samples), n_unique_samples)
+    ), nmds_diag_path, sep = "\t", row.names = FALSE, quote = FALSE)
   }
-
+  all_outputs <- c(all_outputs, nmds_diag_path)
+
   list(
     status = "completed",
     outputs = all_outputs
diff --git a/analysis/06_shared_taxa.R b/analysis/06_shared_taxa.R
index bd6e9c2..28cf412 100644
--- a/analysis/06_shared_taxa.R
+++ b/analysis/06_shared_taxa.R
@@ -12,9 +12,9 @@ run_shared_taxa <- function(context) {
   cfg <- context$config
   shared_dir <- cfg$output$dirs$shared_taxa
   dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
-
+
   # Cohort gate: Requires at least 2 samples
   if (context$mode != "cohort" || length(context$samples) < 2) {
     skip_file <- file.path(shared_dir, "shared_taxa_skipped.tsv")
@@ -31,33 +31,35 @@ run_shared_taxa <- function(context) {
       outputs = skip_file
     ))
   }
-
+
   samples <- context$samples
   unclass_idx <- context$unclass_index
   count_matrix <- context$count_matrix
   class_counts <- count_matrix[-unclass_idx, , drop = FALSE]
   tax_df <- context$taxonomy[-unclass_idx, , drop = FALSE]
-
+
   meta <- context$metadata
   min_count <- cfg$shared_taxa$minimum_count %||% 1L
   group_prev_thresh <- cfg$shared_taxa$group_prevalence %||% 0.5
   analysis_rank <- cfg$shared_taxa$rank %||% "species"
-
+
   # Group by analysis rank
   rk_idx <- which(colnames(context$taxonomy) == analysis_rank)
-  if (length(rk_idx) == 0) rk_idx <- 8L # default species
-
+  if (length(rk_idx) != 1L) {
+    stop(sprintf("Unsupported shared-taxa rank: '%s'", analysis_rank), call. = FALSE)
+  }
+
   rk_prefixes <- apply(tax_df[, 1:rk_idx, drop = FALSE], 1, paste, collapse = ";")
   tax_names <- tax_df[[rk_idx]]
-
+
   agg_counts <- data.frame(TaxonPath = rk_prefixes, Taxon = tax_names, class_counts, check.names = FALSE) %>%
     group_by(TaxonPath, Taxon) %>%
     summarise(across(all_of(samples), sum), .groups = "drop")
-
+
   tax_labels <- agg_counts$Taxon
   mat_counts <- as.matrix(agg_counts[, samples, drop = FALSE])
   rownames(mat_counts) <- agg_counts$TaxonPath
-
+
   # 1. Sample Presence/Absence Matrix
   presence_mat <- (mat_counts >= min_count) * 1L
   presence_df <- data.frame(
@@ -70,16 +72,16 @@ run_shared_taxa <- function(context) {
   presence_file <- file.path(shared_dir, "sample_presence.tsv")
   write.table(presence_df, presence_file, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, presence_file)
-
+
   # 2. Group Prevalence Table & Membership Matrix
   groups <- unique(meta$Group)
   prev_list <- list()
   membership_list <- list()
-
+
   for (grp in groups) {
     grp_samples <- meta$SampleID[meta$Group == grp]
     n_grp <- length(grp_samples)
-
+
     if (n_grp > 0) {
       grp_pres <- presence_mat[, grp_samples, drop = FALSE]
       pos_count <- rowSums(grp_pres)
@@ -88,7 +90,7 @@ run_shared_taxa <- function(context) {
       membership_list[[grp]] <- (prevalence >= group_prev_thresh) * 1L
     }
   }
-
+
   prev_mat <- do.call(cbind, prev_list)
   prev_df <- data.frame(
     TaxonPath = agg_counts$TaxonPath,
@@ -99,7 +101,7 @@ run_shared_taxa <- function(context) {
   prev_file <- file.path(shared_dir, "group_prevalence.tsv")
   write.table(prev_df, prev_file, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, prev_file)
-
+
   mem_mat <- do.call(cbind, membership_list)
   mem_df <- data.frame(
     TaxonPath = agg_counts$TaxonPath,
@@ -111,22 +113,23 @@ run_shared_taxa <- function(context) {
   mem_file <- file.path(shared_dir, "group_membership.tsv")
   write.table(mem_df, mem_file, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, mem_file)
-
+
   # 3. Core & Unique Taxa
   n_groups_present <- rowSums(mem_mat)
   is_core <- (n_groups_present == length(groups))
   is_unique <- (n_groups_present == 1)
-
+
   core_df <- data.frame(
     TaxonPath = agg_counts$TaxonPath[is_core],
     Taxon = tax_labels[is_core],
-    Threshold = sprintf("Present in all %d groups at prevalence >= %.2f", length(groups), group_prev_thresh),
+    Threshold = sprintf("count >= %d per sample; present in all %d groups at prevalence >= %.2f",
+                        min_count, length(groups), group_prev_thresh),
     stringsAsFactors = FALSE
   )
   core_file <- file.path(shared_dir, "core_taxa.tsv")
   write.table(core_df, core_file, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, core_file)
-
+
   unique_group_name <- apply(mem_mat[is_unique, , drop = FALSE], 1, function(row) {
     names(row)[which(row == 1)[1]]
   })
@@ -134,18 +137,19 @@ run_shared_taxa <- function(context) {
     TaxonPath = agg_counts$TaxonPath[is_unique],
     Taxon = tax_labels[is_unique],
     ExclusiveGroup = unique_group_name,
-    Threshold = sprintf("Present exclusively in 1 group at prevalence >= %.2f", group_prev_thresh),
+    Threshold = sprintf("count >= %d per sample; present exclusively in 1 group at prevalence >= %.2f",
+                        min_count, group_prev_thresh),
     stringsAsFactors = FALSE
   )
   unique_file <- file.path(shared_dir, "unique_taxa.tsv")
   write.table(unique_df, unique_file, sep = "\t", row.names = FALSE, quote = FALSE)
   all_outputs <- c(all_outputs, unique_file)
-
+
   # 4. UpSet Plot (When >= 2 non-empty groups exist)
   if (length(groups) >= 2 && sum(colSums(mem_mat) > 0) >= 2) {
     upset_df <- as.data.frame(mem_mat)
     upset_path <- file.path(shared_dir, "06_upset_plot.png")
-
+
     png(upset_path, width = 8, height = 5.5, units = "in", res = 150)
     suppressWarnings(print(UpSetR::upset(
       upset_df,
@@ -158,7 +162,7 @@ run_shared_taxa <- function(context) {
     dev.off()
     all_outputs <- c(all_outputs, upset_path)
   }
-
+
   list(
     status = "completed",
     outputs = all_outputs
diff --git a/analysis/07_kreport_pavian.R b/analysis/07_kreport_pavian.R
index ebe8c7b..f8ccce7 100644
--- a/analysis/07_kreport_pavian.R
+++ b/analysis/07_kreport_pavian.R
@@ -11,104 +11,99 @@ run_kreport <- function(context) {
   cfg <- context$config
   kreport_dir <- cfg$output$dirs$kreport
   dir.create(kreport_dir, recursive = TRUE, showWarnings = FALSE)
-
+
   all_outputs <- character(0)
   cache_file <- cfg$taxonomy$cache
   network_mode <- cfg$taxonomy$network_mode %||% "cache_only"
   unresolved_policy <- cfg$taxonomy$unresolved_policy %||% "warn"
-
+
   unresolved_tsv <- file.path(kreport_dir, "unresolved_taxids.tsv")
+  conflicts_tsv <- file.path(kreport_dir, "taxonomy_conflicts.tsv")
   prov_json <- file.path(kreport_dir, "taxonomy_provenance.json")
-
-  # Invoke Python resolver if cache does not exist or refresh is requested
-  need_resolver <- !file.exists(cache_file) || (network_mode == "refresh")
-
-  py_script <- if (!is.null(context$config$config_dir) && nzchar(context$config$config_dir)) {
-    file.path(context$config$config_dir, "analysis", "utils", "ncbi_taxonomy.py")
-  } else {
-    file.path("analysis", "utils", "ncbi_taxonomy.py")
-  }
-
+  resolved_cache <- file.path(kreport_dir, "resolved_taxonomy_cache.json")
+
+  py_script <- file.path(cfg$pipeline_root, "analysis", "utils", "ncbi_taxonomy.py")
   if (!file.exists(py_script)) {
-    candidates <- c(
-      file.path("analysis", "utils", "ncbi_taxonomy.py"),
-      file.path("..", "analysis", "utils", "ncbi_taxonomy.py"),
-      file.path("..", "..", "analysis", "utils", "ncbi_taxonomy.py")
-    )
-    for (cand in candidates) {
-      if (file.exists(cand)) {
-        py_script <- cand
-        break
-      }
-    }
+    stop(sprintf("Taxonomy resolver script not found: '%s'", py_script), call. = FALSE)
   }
-
-  first_sample <- context$samples[1]
-  first_asgn <- if (!is.null(context$assignments)) context$assignments[[first_sample]] else NULL
-
+
+  python_candidates <- c(Sys.which("python3"), Sys.which("python"))
+  python_cmd <- unname(python_candidates[nzchar(python_candidates)][1])
+  if (is.na(python_cmd) || !nzchar(python_cmd)) {
+    stop("Neither 'python3' nor 'python' was found on PATH.", call. = FALSE)
+  }
+
   cmd_args <- c(
     py_script,
     "--abundance", cfg$input$abundance_table,
+    "--tax-column", cfg$input$tax_column,
     "--cache", cache_file,
+    "--resolved-cache", resolved_cache,
     "--mode", network_mode,
+    "--email-env", cfg$taxonomy$email_env,
+    "--api-key-env", cfg$taxonomy$api_key_env,
     "--unresolved-policy", unresolved_policy,
     "--unresolved-tsv", unresolved_tsv,
+    "--conflicts-tsv", conflicts_tsv,
     "--provenance", prov_json
   )
-  if (!is.null(first_asgn) && file.exists(first_asgn)) {
-    cmd_args <- c(cmd_args, "--assignments", first_asgn)
+  assignment_paths <- unname(unlist(context$assignments, use.names = FALSE))
+  if (length(assignment_paths) > 0L) {
+    cmd_args <- c(cmd_args, as.vector(rbind("--assignments", assignment_paths)))
   }
-
+
   # Run Python script
-  res_code <- system2("python", args = cmd_args)
+  res_code <- system2(python_cmd, args = shQuote(cmd_args))
   if (res_code != 0) {
     stop(sprintf("NCBI taxonomy resolver failed with exit status %d", res_code), call. = FALSE)
   }
-
-  if (file.exists(unresolved_tsv)) all_outputs <- c(all_outputs, unresolved_tsv)
-  if (file.exists(prov_json)) all_outputs <- c(all_outputs, prov_json)
-
-  # Load taxonomy cache
-  tax_cache <- if (file.exists(cache_file)) {
-    jsonlite::fromJSON(cache_file)
-  } else {
-    list()
+
+  required_resolver_outputs <- c(resolved_cache, unresolved_tsv, conflicts_tsv, prov_json)
+  missing_resolver_outputs <- required_resolver_outputs[!file.exists(required_resolver_outputs)]
+  if (length(missing_resolver_outputs) > 0L) {
+    stop(sprintf("Taxonomy resolver omitted expected output(s): %s",
+                 paste(missing_resolver_outputs, collapse = ", ")), call. = FALSE)
   }
-
+  all_outputs <- c(all_outputs, required_resolver_outputs)
+
+  # Load the run-local cache so assignment-derived TaxIDs are available without
+  # mutating the configured source cache in cache_only mode.
+  tax_cache <- jsonlite::fromJSON(resolved_cache, simplifyVector = FALSE)
+
   # Generate .kreport for each sample
   samples <- context$samples
   unclass_idx <- context$unclass_index
   count_matrix <- context$count_matrix
   lineages <- context$taxonomy$TaxonPath
-
+
   resolution_rows <- list()
-
+
   for (s in samples) {
     counts_s <- count_matrix[, s]
-    total_reads <- as.integer(sum(counts_s))
-    uncl_reads <- as.integer(counts_s[unclass_idx])
-
+    total_reads <- as.numeric(sum(counts_s))
+    uncl_reads <- as.numeric(counts_s[unclass_idx])
+
     # Build DFS abundance-sorted tree
     nodes_sorted <- build_kreport_tree(lineages, counts_s)
-
+
     # Validate tree invariants
     validate_kreport_tree(nodes_sorted, total_reads, uncl_reads)
-
+
     # Format 6-column lines
     kreport_lines <- format_kreport_lines(nodes_sorted, total_reads, uncl_reads, tax_cache)
-
+
     # Write .kreport file
     out_file <- file.path(kreport_dir, sprintf("%s.kreport", sanitize_filename(s)))
     writeLines(kreport_lines, out_file)
     all_outputs <- c(all_outputs, out_file)
-
+
     # Collect resolution info
     for (i in seq_len(nrow(nodes_sorted))) {
       p <- nodes_sorted$path[i]
       tid <- tax_cache[[p]]
       if (is.null(tid)) tid <- 0L
-
-      resolution_rows[[p]] <- data.frame(
+
+      resolution_rows[[length(resolution_rows) + 1L]] <- data.frame(
         SampleID = s,
         Depth = nodes_sorted$depth[i],
         RankCode = nodes_sorted$rank_code[i],
@@ -120,7 +115,7 @@ run_kreport <- function(context) {
       )
     }
   }
-
+
   # Export taxonomy resolution summary
   res_summary_file <- file.path(kreport_dir, "taxonomy_resolution.tsv")
   if (length(resolution_rows) > 0) {
@@ -128,7 +123,7 @@ run_kreport <- function(context) {
     write.table(res_df, res_summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
     all_outputs <- c(all_outputs, res_summary_file)
   }
-
+
   list(
     status = "completed",
     outputs = all_outputs
diff --git a/analysis/install_packages.R b/analysis/install_packages.R
index 48f85e0..6ce2406 100644
--- a/analysis/install_packages.R
+++ b/analysis/install_packages.R
@@ -48,7 +48,7 @@ if (length(missing_pkgs) > 0) {
     cat(sprintf("\nAttempting installation of %d missing packages...\n", length(missing_pkgs)))
     repos <- "https://cloud.r-project.org"
     install.packages(missing_pkgs, repos = repos)
-
+
     # Re-verify
     installed_now <- rownames(installed.packages())
     still_missing <- setdiff(REQUIRED_PACKAGES, installed_now)
diff --git a/analysis/utils/cli.R b/analysis/utils/cli.R
index 17e2da2..cc9c76b 100644
--- a/analysis/utils/cli.R
+++ b/analysis/utils/cli.R
@@ -56,7 +56,7 @@ get_cli_parser <- function() {
       help = "Allow overwriting existing output files"
     )
   )
-
+
   optparse::OptionParser(
     usage = "%prog [options]",
     description = "ONT wf-16s Post-Processing Analytical Pipeline",
diff --git a/analysis/utils/config.R b/analysis/utils/config.R
index 34f10e0..cf162b2 100644
--- a/analysis/utils/config.R
+++ b/analysis/utils/config.R
@@ -4,6 +4,131 @@

 suppressMessages(library(yaml))

+`%||%` <- function(x, fallback) {
+  if (is.null(x) || length(x) == 0L) fallback else x
+}
+
+assert_scalar_number <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE,
+                                 lower_open = FALSE, upper_open = FALSE) {
+  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
+  if (valid && integer) valid <- abs(x - round(x)) <= sqrt(.Machine$double.eps)
+  if (valid) valid <- if (lower_open) x > lower else x >= lower
+  if (valid) valid <- if (upper_open) x < upper else x <= upper
+  if (!valid) stop(sprintf("Invalid configuration value '%s'.", name), call. = FALSE)
+  invisible(TRUE)
+}
+
+assert_nonempty_string <- function(x, name) {
+  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
+    stop(sprintf("'%s' must be one non-empty string.", name), call. = FALSE)
+  }
+  invisible(TRUE)
+}
+
+validate_config <- function(cfg) {
+  if (!identical(as.integer(cfg$schema_version), 1L)) {
+    stop("Unsupported schema_version; expected 1.", call. = FALSE)
+  }
+  assert_nonempty_string(cfg$project_name, "project_name")
+  if (!is.character(cfg$mode) || length(cfg$mode) != 1L ||
+      !cfg$mode %in% c("auto", "single", "cohort")) {
+    stop("'mode' must be one of: auto, single, cohort.", call. = FALSE)
+  }
+  assert_scalar_number(cfg$seed, "seed", lower = 0, integer = TRUE)
+
+  assert_nonempty_string(cfg$input$abundance_table, "input.abundance_table")
+  assert_nonempty_string(cfg$input$tax_column, "input.tax_column")
+  assert_nonempty_string(cfg$output$base_dir, "output.base_dir")
+  if (!is.null(cfg$input$aggregate_columns) &&
+      (!is.character(cfg$input$aggregate_columns) || anyNA(cfg$input$aggregate_columns) ||
+       any(!nzchar(cfg$input$aggregate_columns)) || anyDuplicated(cfg$input$aggregate_columns))) {
+    stop("'input.aggregate_columns' must contain unique, non-empty strings.", call. = FALSE)
+  }
+  if (cfg$input$tax_column %in% cfg$input$aggregate_columns) {
+    stop("'input.tax_column' cannot also be listed in 'input.aggregate_columns'.", call. = FALSE)
+  }
+  if (!is.null(cfg$input$include_samples) &&
+      (!is.character(cfg$input$include_samples) || anyNA(cfg$input$include_samples) ||
+       any(!nzchar(cfg$input$include_samples)) || anyDuplicated(cfg$input$include_samples))) {
+    stop("'input.include_samples' must contain unique, non-empty sample IDs.", call. = FALSE)
+  }
+  if (!is.null(cfg$input$assignments)) {
+    if (!is.list(cfg$input$assignments) || is.null(names(cfg$input$assignments)) ||
+        any(!nzchar(names(cfg$input$assignments))) || anyDuplicated(names(cfg$input$assignments))) {
+      stop("'input.assignments' must be a named mapping of unique SampleID -> path.", call. = FALSE)
+    }
+    valid_paths <- vapply(cfg$input$assignments, function(path) {
+      is.character(path) && length(path) == 1L && !is.na(path) && nzchar(trimws(path))
+    }, logical(1))
+    if (!all(valid_paths)) stop("Every 'input.assignments' value must be one non-empty path.", call. = FALSE)
+  }
+
+  assert_scalar_number(cfg$qc$display_min_length, "qc.display_min_length", lower = 1, integer = TRUE)
+  assert_scalar_number(cfg$qc$display_max_length, "qc.display_max_length", lower = 1, integer = TRUE)
+  if (!is.null(cfg$qc$target_min_length)) {
+    assert_scalar_number(cfg$qc$target_min_length, "qc.target_min_length", lower = 1, integer = TRUE)
+  }
+  if (!is.null(cfg$qc$target_max_length)) {
+    assert_scalar_number(cfg$qc$target_max_length, "qc.target_max_length", lower = 1, integer = TRUE)
+  }
+  if (cfg$qc$display_min_length >= cfg$qc$display_max_length) {
+    stop("'qc.display_min_length' must be smaller than 'qc.display_max_length'.", call. = FALSE)
+  }
+  assert_scalar_number(cfg$alpha$rarefaction_points, "alpha.rarefaction_points", lower = 2, integer = TRUE)
+  assert_scalar_number(cfg$alpha$resample_depth, "alpha.resample_depth", lower = 1, integer = TRUE)
+  assert_scalar_number(cfg$alpha$resample_fraction_cap, "alpha.resample_fraction_cap",
+                       lower = 0, upper = 1, lower_open = TRUE)
+  assert_scalar_number(cfg$alpha$resample_iterations, "alpha.resample_iterations", lower = 1, integer = TRUE)
+  assert_scalar_number(cfg$composition$top_n_taxa, "composition.top_n_taxa", lower = 1, integer = TRUE)
+  assert_nonempty_string(cfg$composition$heatmap_rank, "composition.heatmap_rank")
+  assert_nonempty_string(cfg$composition$heatmap_transform, "composition.heatmap_transform")
+  if (!cfg$composition$heatmap_rank %in% c("phylum", "class", "order", "family", "genus", "species")) {
+    stop("'composition.heatmap_rank' is not a supported rank.", call. = FALSE)
+  }
+  if (!cfg$composition$heatmap_transform %in% c("log10_relative", "none")) {
+    stop("'composition.heatmap_transform' must be 'log10_relative' or 'none'.", call. = FALSE)
+  }
+  if (!is.character(cfg$beta$distances) || length(cfg$beta$distances) == 0L ||
+      any(!cfg$beta$distances %in% c("bray", "jaccard")) || anyDuplicated(cfg$beta$distances)) {
+    stop("'beta.distances' must contain unique values drawn from: bray, jaccard.", call. = FALSE)
+  }
+  assert_scalar_number(cfg$beta$permutations, "beta.permutations", lower = 1, integer = TRUE)
+  assert_scalar_number(cfg$beta$minimum_count, "beta.minimum_count", lower = 1, integer = TRUE)
+  if (!is.logical(cfg$beta$resampling$enabled) || length(cfg$beta$resampling$enabled) != 1L ||
+      is.na(cfg$beta$resampling$enabled)) {
+    stop("'beta.resampling.enabled' must be true or false.", call. = FALSE)
+  }
+  assert_scalar_number(cfg$beta$resampling$iterations, "beta.resampling.iterations", lower = 1, integer = TRUE)
+  assert_scalar_number(cfg$beta$resampling$depth_fraction_of_minimum,
+                       "beta.resampling.depth_fraction_of_minimum",
+                       lower = 0, upper = 1, lower_open = TRUE)
+  if (!is.null(cfg$beta$strata_column) &&
+      (!is.character(cfg$beta$strata_column) || length(cfg$beta$strata_column) != 1L ||
+       !nzchar(cfg$beta$strata_column))) {
+    stop("'beta.strata_column' must be null or one non-empty metadata column name.", call. = FALSE)
+  }
+  assert_nonempty_string(cfg$shared_taxa$rank, "shared_taxa.rank")
+  if (!cfg$shared_taxa$rank %in% c("superkingdom", "kingdom", "phylum", "class",
+                                   "order", "family", "genus", "species")) {
+    stop("'shared_taxa.rank' is not a supported rank.", call. = FALSE)
+  }
+  assert_scalar_number(cfg$shared_taxa$minimum_count, "shared_taxa.minimum_count", lower = 1, integer = TRUE)
+  assert_scalar_number(cfg$shared_taxa$group_prevalence, "shared_taxa.group_prevalence",
+                       lower = 0, upper = 1, lower_open = TRUE)
+  assert_nonempty_string(cfg$taxonomy$cache, "taxonomy.cache")
+  assert_nonempty_string(cfg$taxonomy$network_mode, "taxonomy.network_mode")
+  assert_nonempty_string(cfg$taxonomy$unresolved_policy, "taxonomy.unresolved_policy")
+  assert_nonempty_string(cfg$taxonomy$email_env, "taxonomy.email_env")
+  assert_nonempty_string(cfg$taxonomy$api_key_env, "taxonomy.api_key_env")
+  if (!cfg$taxonomy$network_mode %in% c("cache_only", "refresh")) {
+    stop("'taxonomy.network_mode' must be 'cache_only' or 'refresh'.", call. = FALSE)
+  }
+  if (!cfg$taxonomy$unresolved_policy %in% c("warn", "error")) {
+    stop("'taxonomy.unresolved_policy' must be 'warn' or 'error'.", call. = FALSE)
+  }
+  invisible(cfg)
+}
+
 is_absolute_path <- function(p) {
   if (is.null(p) || is.na(p) || length(p) == 0 || !nzchar(p)) return(FALSE)
   grepl("^(/|\\\\|[A-Za-z]:[/\\])", p)
@@ -56,6 +181,7 @@ get_default_config <- function() {
     beta = list(
       distances = c("bray", "jaccard"),
       permutations = 999L,
+      minimum_count = 1L,
       strata_column = NULL,
       resampling = list(
         enabled = FALSE,
@@ -78,11 +204,16 @@ get_default_config <- function() {
   )
 }

-merge_config <- function(default_cfg, user_cfg) {
+merge_config <- function(default_cfg, user_cfg, path = "") {
+  unknown <- setdiff(names(user_cfg), names(default_cfg))
+  if (length(unknown) > 0L) {
+    qualified <- paste0(path, unknown)
+    stop(sprintf("Unknown configuration key(s): %s", paste(qualified, collapse = ", ")), call. = FALSE)
+  }
   merged <- default_cfg
   for (key in names(user_cfg)) {
     if (is.list(user_cfg[[key]]) && is.list(merged[[key]])) {
-      merged[[key]] <- merge_config(merged[[key]], user_cfg[[key]])
+      merged[[key]] <- merge_config(merged[[key]], user_cfg[[key]], paste0(path, key, "."))
     } else {
       merged[[key]] <- user_cfg[[key]]
     }
@@ -94,25 +225,27 @@ load_config <- function(config_path = "config.yml", cli_opts = list()) {
   if (!file.exists(config_path)) {
     stop(sprintf("Configuration file not found: '%s'", config_path), call. = FALSE)
   }
-
+
   config_file_abs <- normalizePath(config_path, winslash = "/", mustWork = TRUE)
   config_dir <- dirname(config_file_abs)
-
+
   raw_yaml <- yaml::read_yaml(config_file_abs)
   default_cfg <- get_default_config()
   cfg <- merge_config(default_cfg, raw_yaml)
-
+
   # CLI overrides
   if (!is.null(cli_opts$output_dir) && nzchar(cli_opts$output_dir)) {
     cfg$output$base_dir <- cli_opts$output_dir
   } else if (!is.null(cli_opts[["output-dir"]]) && nzchar(cli_opts[["output-dir"]])) {
     cfg$output$base_dir <- cli_opts[["output-dir"]]
   }
-
+
   if (isTRUE(cli_opts$refresh_taxonomy) || isTRUE(cli_opts[["refresh-taxonomy"]])) {
     cfg$taxonomy$network_mode <- "refresh"
   }
-
+
+  validate_config(cfg)
+
   cfg$cli <- list(
     validate_only = isTRUE(cli_opts$validate_only) || isTRUE(cli_opts[["validate-only"]]),
     keep_going = isTRUE(cli_opts$keep_going) || isTRUE(cli_opts[["keep-going"]]),
@@ -123,23 +256,23 @@ load_config <- function(config_path = "config.yml", cli_opts = list()) {
       c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport")
     }
   )
-
+
   # Path resolution against config_dir
   cfg$input$abundance_table <- resolve_path(cfg$input$abundance_table, config_dir)
   cfg$input$metadata <- resolve_path(cfg$input$metadata, config_dir)
   cfg$input$params_json <- resolve_path(cfg$input$params_json, config_dir)
-
+
   if (!is.null(cfg$input$assignments) && is.list(cfg$input$assignments)) {
     for (s in names(cfg$input$assignments)) {
       cfg$input$assignments[[s]] <- resolve_path(cfg$input$assignments[[s]], config_dir)
     }
   }
-
+
   cfg$taxonomy$cache <- resolve_path(cfg$taxonomy$cache, config_dir)
-
+
   # Resolve base output directory
   cfg$output$base_dir <- resolve_path(cfg$output$base_dir, config_dir)
-
+
   # Derive all module output directories
   base_out <- cfg$output$base_dir
   cfg$output$dirs <- list(
@@ -151,13 +284,13 @@ load_config <- function(config_path = "config.yml", cli_opts = list()) {
     shared_taxa = file.path(base_out, "06_Shared_Taxa"),
     kreport = file.path(base_out, "07_Kreport")
   )
-
+
   cfg$output$manifest_file <- file.path(base_out, "run_manifest.json")
   cfg$output$resolved_config_file <- file.path(base_out, "resolved_config.yml")
   cfg$output$session_info_file <- file.path(base_out, "session_info.txt")
-
+
   cfg$config_file <- config_file_abs
   cfg$config_dir <- config_dir
-
+
   cfg
 }
diff --git a/analysis/utils/io.R b/analysis/utils/io.R
index 0592f7f..d862d2f 100644
--- a/analysis/utils/io.R
+++ b/analysis/utils/io.R
@@ -24,28 +24,28 @@ validate_sample_ids <- function(sample_ids) {
   if (length(sample_ids) == 0) {
     stop("Abundance table validation error: No sample columns detected.", call. = FALSE)
   }
-
+
   if (any(is.na(sample_ids)) || any(sample_ids == "")) {
     stop("Sample ID validation error: Empty or NA sample ID detected.", call. = FALSE)
   }
-
+
   if (any(sample_ids %in% c(".", ".."))) {
     stop("Sample ID validation error: Sample ID cannot be '.' or '..'.", call. = FALSE)
   }
-
+
   if (any(grepl("[/\\\\]", sample_ids))) {
     stop("Sample ID validation error: Sample ID cannot contain path separators ('/' or '\\').", call. = FALSE)
   }
-
+
   if (any(grepl("[[:cntrl:]]", sample_ids))) {
     stop("Sample ID validation error: Sample ID cannot contain control characters.", call. = FALSE)
   }
-
+
   sanitized <- vapply(sample_ids, sanitize_filename, character(1))
   if (any(duplicated(sanitized))) {
     stop("Sample ID validation error: Sample IDs collide after filename sanitization.", call. = FALSE)
   }
-
+
   invisible(TRUE)
 }

@@ -53,36 +53,36 @@ read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("tota
   if (!file.exists(path)) {
     stop(sprintf("Abundance table not found: '%s'", path), call. = FALSE)
   }
-
+
   # Read header first
   raw_lines <- readLines(path, n = 5)
   if (length(raw_lines) == 0) {
     stop(sprintf("Abundance table is empty: '%s'", path), call. = FALSE)
   }
-
+
   header <- strsplit(raw_lines[1], "\t")[[1]]
   if (any(duplicated(header))) {
-    stop(sprintf("Abundance table has duplicate column names: %s",
+    stop(sprintf("Abundance table has duplicate column names: %s",
                  paste(header[duplicated(header)], collapse = ", ")), call. = FALSE)
   }
-
+
   if (!tax_col %in% header) {
     stop(sprintf("Abundance table missing configured tax column '%s'. Columns found: %s",
                  tax_col, paste(header, collapse = ", ")), call. = FALSE)
   }
-
+
   raw_df <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
-
+
   # Check rows
   if (nrow(raw_df) == 0) {
     stop("Abundance table has 0 data rows.", call. = FALSE)
   }
-
+
   all_cols <- colnames(raw_df)
   all_sample_cols <- setdiff(all_cols, c(tax_col, aggregate_cols))
-
+
   validate_sample_ids(all_sample_cols)
-
+
   # Validate counts for all sample columns
   for (sc in all_sample_cols) {
     vals <- raw_df[[sc]]
@@ -105,12 +105,17 @@ read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("tota
       stop(sprintf("Column '%s' contains non-integer count values.", sc), call. = FALSE)
     }
   }
-
+
   # Validate aggregate columns (e.g. 'total') if present
   for (ac in aggregate_cols) {
     if (ac %in% all_cols) {
       actual_sum <- rowSums(as.matrix(raw_df[, all_sample_cols, drop = FALSE]))
-      stated_total <- as.numeric(raw_df[[ac]])
+      stated_total <- suppressWarnings(as.numeric(raw_df[[ac]]))
+      if (anyNA(stated_total) || any(!is.finite(stated_total)) || any(stated_total < 0) ||
+          any(abs(stated_total - round(stated_total)) > 1e-6)) {
+        stop(sprintf("Aggregate column '%s' must contain finite, non-negative integer counts.", ac),
+             call. = FALSE)
+      }
       if (any(abs(actual_sum - stated_total) > 1e-4)) {
         diff_idx <- which(abs(actual_sum - stated_total) > 1e-4)[1]
         stop(sprintf(
@@ -120,7 +125,7 @@ read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("tota
       }
     }
   }
-
+
   # Select requested samples
   selected_samples <- if (!is.null(include_samples) && length(include_samples) > 0) {
     missing_sel <- setdiff(include_samples, all_sample_cols)
@@ -132,14 +137,17 @@ read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("tota
   } else {
     all_sample_cols
   }
-
+
   # Validate lineages
   lineages <- raw_df[[tax_col]]
+  if (anyNA(lineages) || any(!nzchar(lineages)) || any(lineages != trimws(lineages))) {
+    stop("Taxonomy lineages must be non-empty and must not have leading/trailing whitespace.", call. = FALSE)
+  }
   if (any(duplicated(lineages))) {
     dup <- lineages[duplicated(lineages)][1]
     stop(sprintf("Duplicate full lineage detected in abundance table: '%s'", dup), call. = FALSE)
   }
-
+
   parsed_lineages <- strsplit(lineages, ";")
   field_counts <- vapply(parsed_lineages, length, integer(1))
   if (any(field_counts != 8)) {
@@ -149,7 +157,7 @@ read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("tota
       bad_idx + 1, field_counts[bad_idx], lineages[bad_idx]
     ), call. = FALSE)
   }
-
+
   # Check unclassified rows
   unclass_indices <- which(vapply(parsed_lineages, function(x) x[1] == "Unclassified", logical(1)))
   if (length(unclass_indices) == 0) {
@@ -159,12 +167,16 @@ read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("tota
     stop(sprintf("Abundance table validation error: Multiple (%d) 'Unclassified' rows detected.",
                  length(unclass_indices)), call. = FALSE)
   }
-
+  expected_unclassified <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
+  if (!identical(lineages[unclass_indices], expected_unclassified)) {
+    stop(sprintf("Unclassified lineage must be exactly '%s'.", expected_unclassified), call. = FALSE)
+  }
+
   # Build count matrix (taxa x samples)
   count_mat <- as.matrix(raw_df[, selected_samples, drop = FALSE])
   rownames(count_mat) <- lineages
   mode(count_mat) <- "numeric"
-
+
   # Check sample read counts
   for (s in selected_samples) {
     tot_reads <- sum(count_mat[, s])
@@ -177,13 +189,13 @@ read_abundance_table <- function(path, tax_col = "tax", aggregate_cols = c("tota
       stop(sprintf("Sample '%s' has 0 classified reads.", s), call. = FALSE)
     }
   }
-
+
   # Parse taxonomy data frame with 8 ranks
   tax_mat <- do.call(rbind, parsed_lineages)
   colnames(tax_mat) <- RANKS_8
   tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
   tax_df$TaxonPath <- lineages
-
+
   list(
     count_matrix = count_mat,
     taxonomy = tax_df,
@@ -197,10 +209,24 @@ read_assignments_file <- function(path, sample_id, expected_total = NULL, expect
   if (!file.exists(path)) {
     stop(sprintf("Assignments file for sample '%s' not found: '%s'", sample_id, path), call. = FALSE)
   }
-
-  # Schema: 5 tab-separated fields: status, read_id, taxid, len_field, lineage
+
+  # Validate the physical five-field schema before read.delim() can normalize it.
+  lines <- readLines(path, warn = FALSE)
+  nonempty_line_numbers <- which(nzchar(lines))
+  if (length(nonempty_line_numbers) == 0L) {
+    stop(sprintf("Assignments file for sample '%s' is empty.", sample_id), call. = FALSE)
+  }
+  field_counts <- lengths(strsplit(lines[nonempty_line_numbers], "\t", fixed = TRUE))
+  if (any(field_counts != 5L)) {
+    bad <- which(field_counts != 5L)[1]
+    stop(sprintf(
+      "Assignments file for sample '%s' has %d fields at line %d; expected exactly 5.",
+      sample_id, field_counts[bad], nonempty_line_numbers[bad]
+    ), call. = FALSE)
+  }
+
   raw_reads <- read.delim(
-    path,
+    text = paste(lines[nonempty_line_numbers], collapse = "\n"),
     header = FALSE,
     sep = "\t",
     col.names = c("status", "read_id", "taxid", "len_field", "lineage"),
@@ -209,57 +235,68 @@ read_assignments_file <- function(path, sample_id, expected_total = NULL, expect
     quote = "",
     fill = FALSE
   )
-
-  if (nrow(raw_reads) == 0) {
-    stop(sprintf("Assignments file for sample '%s' is empty.", sample_id), call. = FALSE)
+
+  if (any(is.na(raw_reads$read_id) | !nzchar(raw_reads$read_id))) {
+    bad_idx <- which(is.na(raw_reads$read_id) | !nzchar(raw_reads$read_id))[1]
+    stop(sprintf("Assignments file for sample '%s' has an empty read ID at line %d",
+                 sample_id, nonempty_line_numbers[bad_idx]), call. = FALSE)
   }
-
+
   # Check unique read IDs
   if (any(duplicated(raw_reads$read_id))) {
     dup_id <- raw_reads$read_id[duplicated(raw_reads$read_id)][1]
     stop(sprintf("Assignments file for sample '%s' contains duplicate read ID: '%s'", sample_id, dup_id), call. = FALSE)
   }
-
+
   # Validate status
   valid_statuses <- raw_reads$status %in% c("C", "U")
   if (!all(valid_statuses)) {
     bad_idx <- which(!valid_statuses)[1]
     stop(sprintf("Assignments file for sample '%s' has invalid status '%s' at line %d",
-                 sample_id, raw_reads$status[bad_idx], bad_idx), call. = FALSE)
+                 sample_id, raw_reads$status[bad_idx], nonempty_line_numbers[bad_idx]), call. = FALSE)
   }
-
+
   # Parse TaxID
-  taxid_num <- suppressWarnings(as.integer(raw_reads$taxid))
-  if (any(is.na(taxid_num))) {
-    bad_idx <- which(is.na(taxid_num))[1]
+  valid_taxid <- grepl("^[0-9]+$", raw_reads$taxid)
+  taxid_num <- suppressWarnings(as.numeric(raw_reads$taxid))
+  if (any(!valid_taxid | is.na(taxid_num) | !is.finite(taxid_num))) {
+    bad_idx <- which(!valid_taxid | is.na(taxid_num) | !is.finite(taxid_num))[1]
     stop(sprintf("Assignments file for sample '%s' has non-integer TaxID '%s' at line %d",
-                 sample_id, raw_reads$taxid[bad_idx], bad_idx), call. = FALSE)
+                 sample_id, raw_reads$taxid[bad_idx], nonempty_line_numbers[bad_idx]), call. = FALSE)
   }
   raw_reads$taxid <- taxid_num
-
+
+  inconsistent <- raw_reads$status == "U" & raw_reads$taxid > 0
+  if (any(inconsistent)) {
+    bad_idx <- which(inconsistent)[1]
+    stop(sprintf("Assignments file for sample '%s' has status U with positive TaxID at line %d",
+                 sample_id, nonempty_line_numbers[bad_idx]), call. = FALSE)
+  }
+
   # Parse read length defensively: single integer or last numeric part of pipe-delimited string
   # Examples: "0|1481" -> 1481; "1500" -> 1500
   parsed_len <- vapply(raw_reads$len_field, function(lf) {
     if (is.null(lf) || is.na(lf) || lf == "") return(NA_integer_)
     parts <- strsplit(lf, "\\|")[[1]]
     last_p <- parts[length(parts)]
+    if (!grepl("^[0-9]+$", last_p)) return(NA_integer_)
     suppressWarnings(as.integer(last_p))
   }, integer(1), USE.NAMES = FALSE)
-
-  if (any(is.na(parsed_len))) {
-    bad_idx <- which(is.na(parsed_len))[1]
+
+  if (any(is.na(parsed_len) | !is.finite(parsed_len) | parsed_len <= 0)) {
+    bad_idx <- which(is.na(parsed_len) | !is.finite(parsed_len) | parsed_len <= 0)[1]
     stop(sprintf("Assignments file for sample '%s' has malformed length field '%s' at line %d",
-                 sample_id, raw_reads$len_field[bad_idx], bad_idx), call. = FALSE)
+                 sample_id, raw_reads$len_field[bad_idx], nonempty_line_numbers[bad_idx]), call. = FALSE)
   }
   raw_reads$read_length <- parsed_len
-
+
   # Effective classification: taxid > 0
   raw_reads$effective_classified <- (raw_reads$taxid > 0)
-
+
   n_total <- nrow(raw_reads)
   n_eff_class <- sum(raw_reads$effective_classified)
   n_eff_unclass <- n_total - n_eff_class
-
+
   # Reconcile against abundance expectations
   if (!is.null(expected_total) && n_total != expected_total) {
     stop(sprintf(
@@ -279,7 +316,7 @@ read_assignments_file <- function(path, sample_id, expected_total = NULL, expect
       sample_id, n_eff_unclass, expected_unclassified
     ), call. = FALSE)
   }
-
+
   raw_reads
 }

@@ -288,25 +325,41 @@ read_metadata_table <- function(path, selected_samples) {
   if (!file.exists(path)) {
     stop(sprintf("Metadata file not found: '%s'", path), call. = FALSE)
   }
-
+
   meta <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
-
+
+  if (anyDuplicated(colnames(meta))) {
+    stop("Metadata table contains duplicate column names.", call. = FALSE)
+  }
+
   if (!"SampleID" %in% colnames(meta)) {
     stop("Metadata table must contain a 'SampleID' column.", call. = FALSE)
   }
   if (!"Group" %in% colnames(meta)) {
     stop("Metadata table must contain a 'Group' column.", call. = FALSE)
   }
-
+
+  invalid_sample <- is.na(meta$SampleID) | !nzchar(trimws(meta$SampleID))
+  invalid_group <- is.na(meta$Group) | !nzchar(trimws(meta$Group))
+  if (any(invalid_sample)) {
+    stop(sprintf("Metadata contains an empty SampleID at data row %d.", which(invalid_sample)[1]), call. = FALSE)
+  }
+  if (any(invalid_group)) {
+    stop(sprintf("Metadata contains an empty Group at data row %d.", which(invalid_group)[1]), call. = FALSE)
+  }
+  if (any(meta$SampleID != trimws(meta$SampleID)) || any(meta$Group != trimws(meta$Group))) {
+    stop("Metadata SampleID and Group values must not have leading or trailing whitespace.", call. = FALSE)
+  }
+
   if (any(duplicated(meta$SampleID))) {
     dup_ids <- meta$SampleID[duplicated(meta$SampleID)]
     stop(sprintf("Metadata contains duplicate SampleID values: %s", paste(unique(dup_ids), collapse = ", ")), call. = FALSE)
   }
-
+
   meta_ids <- meta$SampleID
   missing_in_meta <- setdiff(selected_samples, meta_ids)
   extra_in_meta <- setdiff(meta_ids, selected_samples)
-
+
   if (length(missing_in_meta) > 0 || length(extra_in_meta) > 0) {
     msg <- "Metadata SampleID mismatch with selected abundance samples:"
     if (length(missing_in_meta) > 0) {
@@ -317,7 +370,7 @@ read_metadata_table <- function(path, selected_samples) {
     }
     stop(msg, call. = FALSE)
   }
-
+
   # Align metadata to exact order of selected_samples
   meta_aligned <- meta[match(selected_samples, meta$SampleID), , drop = FALSE]
   rownames(meta_aligned) <- selected_samples
@@ -332,12 +385,12 @@ build_context <- function(cfg) {
     aggregate_cols = cfg$input$aggregate_columns,
     include_samples = cfg$input$include_samples
   )
-
+
   selected_samples <- ab_res$samples
   count_matrix <- ab_res$count_matrix
   taxonomy_df <- ab_res$taxonomy
   unclass_idx <- ab_res$unclass_index
-
+
   # Calculate per-sample read stats
   sample_stats <- data.frame(
     SampleID = selected_samples,
@@ -346,7 +399,7 @@ build_context <- function(cfg) {
     ClassifiedReads = colSums(count_matrix[-unclass_idx, , drop = FALSE]),
     stringsAsFactors = FALSE
   )
-
+
   # 2. Mode resolution
   configured_mode <- cfg$mode
   resolved_mode <- if (configured_mode == "auto") {
@@ -356,30 +409,68 @@ build_context <- function(cfg) {
   } else {
     stop(sprintf("Invalid mode '%s' in configuration. Must be 'auto', 'single', or 'cohort'.", configured_mode), call. = FALSE)
   }
-
+
   if (resolved_mode == "single" && length(selected_samples) != 1) {
     stop(sprintf("Mode is 'single' but %d samples are selected.", length(selected_samples)), call. = FALSE)
   }
-
+  if (resolved_mode == "cohort" && length(selected_samples) < 2L) {
+    stop(sprintf("Mode is 'cohort' but only %d sample is selected; at least 2 are required.",
+                 length(selected_samples)), call. = FALSE)
+  }
+
   # 3. Read metadata
   metadata <- read_metadata_table(cfg$input$metadata, selected_samples)
-
+
   if (resolved_mode == "cohort" && is.null(metadata)) {
     stop("Cohort mode requires a metadata table mapping SampleID to Group.", call. = FALSE)
   }
-
+
   # 4. Assignments mapping
   assignments_map <- cfg$input$assignments
   if (!is.null(assignments_map) && !is.list(assignments_map)) {
     stop("Config 'input.assignments' must be a mapping of SampleID -> path or null.", call. = FALSE)
   }
-
+  assignment_data <- list()
+  retain_assignment_rows <- !is.null(cfg$cli) && "qc" %in% cfg$cli$modules && !isTRUE(cfg$cli$validate_only)
+  if (!is.null(assignments_map)) {
+    extra_assignment_ids <- setdiff(names(assignments_map), selected_samples)
+    if (length(extra_assignment_ids) > 0L) {
+      stop(sprintf("Assignments configured for unselected/unknown sample(s): %s",
+                   paste(extra_assignment_ids, collapse = ", ")), call. = FALSE)
+    }
+    for (s in names(assignments_map)) {
+      stat_row <- sample_stats[sample_stats$SampleID == s, , drop = FALSE]
+      parsed_assignments <- read_assignments_file(
+        assignments_map[[s]], s,
+        expected_total = stat_row$TotalReads[1],
+        expected_classified = stat_row$ClassifiedReads[1],
+        expected_unclassified = stat_row$UnclassifiedReads[1]
+      )
+      if (retain_assignment_rows) assignment_data[[s]] <- parsed_assignments
+      rm(parsed_assignments)
+    }
+  }
+
+  if (resolved_mode == "cohort" && !is.null(cfg$beta$strata_column)) {
+    strata_col <- cfg$beta$strata_column
+    if (!strata_col %in% colnames(metadata)) {
+      stop(sprintf("Configured beta.strata_column '%s' is absent from metadata.", strata_col), call. = FALSE)
+    }
+    strata <- metadata[[strata_col]]
+    if (anyNA(strata) || any(!nzchar(trimws(as.character(strata))))) {
+      stop(sprintf("Metadata strata column '%s' contains missing/empty values.", strata_col), call. = FALSE)
+    }
+    if (length(unique(strata)) < 2L) {
+      stop(sprintf("Metadata strata column '%s' must contain at least two strata.", strata_col), call. = FALSE)
+    }
+  }
+
   # 5. Read params.json if available
   params <- NULL
   if (!is.null(cfg$input$params_json) && file.exists(cfg$input$params_json)) {
     params <- suppressWarnings(tryCatch(jsonlite::fromJSON(cfg$input$params_json), error = function(e) NULL))
   }
-
+
   # 6. File hashes
   file_hashes <- list(
     abundance_table = compute_file_hash(cfg$input$abundance_table),
@@ -392,7 +483,7 @@ build_context <- function(cfg) {
       file_hashes[[paste0("assignment_", s)]] <- compute_file_hash(assignments_map[[s]])
     }
   }
-
+
   list(
     config = cfg,
     mode = resolved_mode,
@@ -403,6 +494,7 @@ build_context <- function(cfg) {
     sample_stats = sample_stats,
     metadata = metadata,
     assignments = assignments_map,
+    assignment_data = assignment_data,
     params = params,
     file_hashes = file_hashes,
     warnings = character(0)
diff --git a/analysis/utils/kreport.R b/analysis/utils/kreport.R
index bf910ec..6059c28 100644
--- a/analysis/utils/kreport.R
+++ b/analysis/utils/kreport.R
@@ -9,19 +9,19 @@ RANK_CODES_8 <- c("D", "K", "P", "C", "O", "F", "G", "S")

 build_kreport_tree <- function(lineages_str, counts) {
   lineages <- strsplit(lineages_str, ";")
-  counts <- as.integer(round(counts))
-
+  counts <- round(as.numeric(counts))
+
   node_env <- new.env(hash = TRUE, parent = emptyenv())
-
+
   for (i in seq_along(lineages)) {
     lin <- lineages[[i]]
     cnt <- counts[i]
-
+
     if (lin[1] == "Unclassified" || cnt == 0) next
-
+
     for (j in seq_along(lin)) {
       path <- paste(lin[1:j], collapse = ";")
-
+
       if (exists(path, envir = node_env)) {
         node <- get(path, envir = node_env)
         node$reads_clade <- node$reads_clade + cnt
@@ -40,7 +40,7 @@ build_kreport_tree <- function(lineages_str, counts) {
       }
     }
   }
-
+
   all_paths <- ls(node_env)
   if (length(all_paths) == 0) {
     return(data.frame(
@@ -50,10 +50,10 @@ build_kreport_tree <- function(lineages_str, counts) {
       stringsAsFactors = FALSE
     ))
   }
-
+
   nodes_list <- lapply(all_paths, function(p) get(p, envir = node_env))
   nodes_df <- do.call(rbind.data.frame, c(nodes_list, stringsAsFactors = FALSE))
-
+
   # Depth-first search sorting with abundance tie-breaking
   dfs_order <- function(parent) {
     children <- nodes_df[nodes_df$parent_path == parent, , drop = FALSE]
@@ -66,7 +66,7 @@ build_kreport_tree <- function(lineages_str, counts) {
     }
     res
   }
-
+
   ordered_paths <- dfs_order("")
   nodes_sorted <- nodes_df[match(ordered_paths, nodes_df$path), ]
   rownames(nodes_sorted) <- NULL
@@ -75,60 +75,60 @@ build_kreport_tree <- function(lineages_str, counts) {

 validate_kreport_tree <- function(nodes_df, total_reads, uncl_reads) {
   cl_reads <- total_reads - uncl_reads
-
+
   # 1. Total reads check
   root_nodes <- nodes_df[nodes_df$parent_path == "", , drop = FALSE]
   sum_root_clade <- sum(root_nodes$reads_clade)
-
+
   if (sum_root_clade != cl_reads) {
     stop(sprintf("Kreport tree validation error: root clades sum (%d) != total classified reads (%d)",
                  sum_root_clade, cl_reads), call. = FALSE)
   }
-
+
   if (uncl_reads + sum_root_clade != total_reads) {
     stop(sprintf("Kreport tree validation error: unclassified (%d) + root (%d) != total reads (%d)",
                  uncl_reads, sum_root_clade, total_reads), call. = FALSE)
   }
-
+
   # 2. Clade = direct + sum(child clades) check
   for (i in seq_len(nrow(nodes_df))) {
     p <- nodes_df$path[i]
     clade_cnt <- nodes_df$reads_clade[i]
     direct_cnt <- nodes_df$reads_taxon[i]
-
+
     children <- nodes_df[nodes_df$parent_path == p, , drop = FALSE]
     child_sum <- if (nrow(children) > 0) sum(children$reads_clade) else 0L
-
+
     if (clade_cnt != (direct_cnt + child_sum)) {
       stop(sprintf("Kreport tree validation error at '%s': clade (%d) != direct (%d) + child sum (%d)",
                    p, clade_cnt, direct_cnt, child_sum), call. = FALSE)
     }
   }
-
+
   invisible(TRUE)
 }

 format_kreport_lines <- function(nodes_sorted, total_reads, uncl_reads, taxid_cache = list()) {
   cl_reads <- total_reads - uncl_reads
-
+
   lines <- character(nrow(nodes_sorted) + 2)
-
+
   # Line 1: unclassified
-  lines[1] <- sprintf("%.2f\t%d\t%d\tU\t0\tunclassified",
+  lines[1] <- sprintf("%.2f\t%.0f\t%.0f\tU\t0\tunclassified",
                       100 * uncl_reads / total_reads, uncl_reads, uncl_reads)
-
+
   # Line 2: root
-  lines[2] <- sprintf("%.2f\t%d\t%d\tR\t1\troot",
+  lines[2] <- sprintf("%.2f\t%.0f\t%.0f\tR\t1\troot",
                       100 * cl_reads / total_reads, cl_reads, 0L)
-
+
   for (i in seq_len(nrow(nodes_sorted))) {
     indent <- strrep("  ", nodes_sorted$depth[i])
     p <- nodes_sorted$path[i]
     taxid <- taxid_cache[[p]]
     if (is.null(taxid)) taxid <- 0L
-
+
     pct <- 100 * nodes_sorted$reads_clade[i] / total_reads
-    lines[i + 2] <- sprintf("%.2f\t%d\t%d\t%s\t%s\t%s%s",
+    lines[i + 2] <- sprintf("%.2f\t%.0f\t%.0f\t%s\t%s\t%s%s",
                             pct,
                             nodes_sorted$reads_clade[i],
                             nodes_sorted$reads_taxon[i],
@@ -137,6 +137,6 @@ format_kreport_lines <- function(nodes_sorted, total_reads, uncl_reads, taxid_ca
                             indent,
                             nodes_sorted$name[i])
   }
-
+
   lines
 }
diff --git a/analysis/utils/metrics.R b/analysis/utils/metrics.R
index 89d6f26..564f098 100644
--- a/analysis/utils/metrics.R
+++ b/analysis/utils/metrics.R
@@ -11,7 +11,7 @@ calc_alpha_indices <- function(counts) {
   counts <- as.integer(round(counts[counts > 0]))
   total_classified <- sum(counts)
   S <- length(counts)
-
+
   if (S == 0 || total_classified == 0) {
     return(data.frame(
       Metric = c("Observed species richness (S)", "Chao1 (estimated richness)",
@@ -22,22 +22,22 @@ calc_alpha_indices <- function(counts) {
       stringsAsFactors = FALSE
     ))
   }
-
+
   shannon <- vegan::diversity(counts, index = "shannon")
   simpson <- vegan::diversity(counts, index = "simpson")
   invsimpson <- vegan::diversity(counts, index = "invsimpson")
   ens <- exp(shannon)
   pielou <- if (S > 1) shannon / log(S) else NA_real_
   bp_dominance <- max(counts) / total_classified
-
+
   chao1 <- tryCatch({
     as.numeric(vegan::estimateR(counts)["S.chao1"])
   }, error = function(e) NA_real_)
-
+
   fisher_alpha <- tryCatch({
     as.numeric(vegan::fisher.alpha(counts))
   }, error = function(e) NA_real_)
-
+
   data.frame(
     Metric = c("Observed species richness (S)", "Chao1 (estimated richness)",
                "Shannon (H)", "Effective number of species (e^H)",
@@ -51,16 +51,16 @@ calc_alpha_indices <- function(counts) {
 calc_analytical_rarefaction <- function(counts, n_points = 25) {
   counts <- as.integer(round(counts[counts > 0]))
   total_classified <- sum(counts)
-
-  if (total_classified < 10) {
+
+  if (total_classified < 1) {
     return(data.frame(depth = integer(0), mean_richness = numeric(0), sd_richness = numeric(0)))
   }
-
+
   start_depth <- min(100L, total_classified)
   depth_points <- unique(round(seq(start_depth, total_classified, length.out = n_points)))
-
+
   rare_res <- vegan::rarefy(counts, sample = depth_points, se = TRUE)
-
+
   data.frame(
     depth = depth_points,
     mean_richness = as.numeric(rare_res[1, ]),
@@ -70,10 +70,13 @@ calc_analytical_rarefaction <- function(counts, n_points = 25) {

 calc_rarefaction_resamples <- function(counts, subsample_depth, n_iterations = 100, seed = 42) {
   counts <- as.integer(round(counts[counts > 0]))
+  if (length(counts) == 0L || subsample_depth < 1L || subsample_depth > sum(counts)) {
+    stop("Invalid rarefaction resampling depth for the supplied counts.", call. = FALSE)
+  }
   set.seed(seed)
-
+
   count_mat <- matrix(counts, nrow = 1)
-
+
   res_list <- vector("list", n_iterations)
   for (i in seq_len(n_iterations)) {
     sub <- vegan::rrarefy(count_mat, subsample_depth)[1, ]
@@ -85,7 +88,7 @@ calc_rarefaction_resamples <- function(counts, subsample_depth, n_iterations = 1
     ens_sub <- exp(shannon_sub)
     pielou_sub <- if (S_sub > 1) shannon_sub / log(S_sub) else NA_real_
     chao1_sub <- tryCatch(as.numeric(vegan::estimateR(sub_counts)["S.chao1"]), error = function(e) NA_real_)
-
+
     res_list[[i]] <- data.frame(
       iteration = i,
       subsample_depth = subsample_depth,
@@ -98,6 +101,6 @@ calc_rarefaction_resamples <- function(counts, subsample_depth, n_iterations = 1
       pielou = pielou_sub
     )
   }
-
+
   do.call(rbind, res_list)
 }
diff --git a/analysis/utils/ncbi_taxonomy.py b/analysis/utils/ncbi_taxonomy.py
index 4b37d0a..3f0afef 100644
--- a/analysis/utils/ncbi_taxonomy.py
+++ b/analysis/utils/ncbi_taxonomy.py
@@ -1,239 +1,260 @@
 #!/usr/bin/env python3
-# =============================================================================
-# NCBI Taxonomy Resolver & Cache Manager
-# =============================================================================
+"""Resolve NCBI TaxIDs without mutating the source cache in offline mode."""

+import argparse
+import csv
+import hashlib
+import json
 import os
 import sys
-import json
-import argparse
+import tempfile
 import time
-import urllib.request
 import urllib.parse
-import xml.etree.ElementTree as ET
-import tempfile
-import hashlib
+import urllib.request
 from collections import Counter

 TOOL_NAME = "ont_wf16s_postprocess"

+
 def normalize_abundance_path_to_7(path):
-    """Normalize 8-rank abundance path (omits kingdom) to match 7-rank minimap2 lineage."""
+    """Omit the eight-rank abundance schema's kingdom field."""
     parts = path.split(";")
-    if len(parts) >= 8:
-        # omit rank index 1 (kingdom)
-        return "|".join([parts[0]] + parts[2:8])
-    return "|".join(parts)
+    return "|".join([parts[0]] + parts[2:8]) if len(parts) == 8 else "|".join(parts)
+

 def compute_sha256(filepath):
     if not filepath or not os.path.exists(filepath):
         return None
-    h = hashlib.sha256()
-    with open(filepath, "rb") as f:
-        while chunk := f.read(65536):
-            h.update(chunk)
-    return h.hexdigest()
+    digest = hashlib.sha256()
+    with open(filepath, "rb") as handle:
+        for chunk in iter(lambda: handle.read(65536), b""):
+            digest.update(chunk)
+    return digest.hexdigest()
+
+
+def atomic_write_json(path, payload):
+    directory = os.path.dirname(os.path.abspath(path))
+    os.makedirs(directory, exist_ok=True)
+    temp_path = None
+    try:
+        with tempfile.NamedTemporaryFile("w", dir=directory, delete=False, encoding="utf-8") as handle:
+            json.dump(payload, handle, indent=2, sort_keys=True)
+            handle.write("\n")
+            temp_path = handle.name
+        os.replace(temp_path, path)
+    finally:
+        if temp_path and os.path.exists(temp_path):
+            os.unlink(temp_path)
+
+
+def load_cache(path):
+    if not os.path.exists(path):
+        return {}
+    try:
+        with open(path, "r", encoding="utf-8") as handle:
+            cache = json.load(handle)
+    except (OSError, json.JSONDecodeError) as exc:
+        raise ValueError(f"Could not read taxonomy cache safely: {exc}") from exc
+    if not isinstance(cache, dict):
+        raise ValueError("Taxonomy cache root must be a JSON object.")
+    for taxon_path, taxid in cache.items():
+        if not isinstance(taxon_path, str) or not isinstance(taxid, int) or isinstance(taxid, bool) or taxid < 0:
+            raise ValueError(f"Invalid cache entry for {taxon_path!r}: expected a non-negative integer TaxID.")
+    return cache
+
+
+def read_abundance_paths(path, tax_column):
+    with open(path, "r", encoding="utf-8", newline="") as handle:
+        reader = csv.DictReader(handle, delimiter="\t")
+        if reader.fieldnames is None or tax_column not in reader.fieldnames:
+            raise ValueError(f"Abundance table does not contain tax column {tax_column!r}.")
+        paths = [row[tax_column].strip() for row in reader if row.get(tax_column, "").strip()]
+    return list(dict.fromkeys(paths))
+
+
+def read_assignment_taxids(paths):
+    lineage_to_taxids = {}
+    for path in paths:
+        with open(path, "r", encoding="utf-8") as handle:
+            for line_number, line in enumerate(handle, start=1):
+                if not line.rstrip("\r\n"):
+                    continue
+                fields = line.rstrip("\r\n").split("\t")
+                if len(fields) != 5:
+                    raise ValueError(f"{path}:{line_number}: expected exactly 5 assignment fields.")
+                taxid_text = fields[2].strip()
+                lineage = fields[4].strip()
+                if not taxid_text.isdigit():
+                    raise ValueError(f"{path}:{line_number}: invalid TaxID {taxid_text!r}.")
+                taxid = int(taxid_text)
+                if taxid > 0 and lineage:
+                    lineage_to_taxids.setdefault(lineage, []).append(taxid)
+
+    resolved = {}
+    conflicts = []
+    for lineage in sorted(lineage_to_taxids):
+        counts = Counter(lineage_to_taxids[lineage])
+        maximum = max(counts.values())
+        winner = min(taxid for taxid, count in counts.items() if count == maximum)
+        resolved[lineage] = winner
+        if len(counts) > 1:
+            conflicts.append({
+                "lineage": lineage,
+                "winner_taxid": winner,
+                "counts": {str(key): counts[key] for key in sorted(counts)},
+            })
+    return resolved, conflicts
+
+
+def find_unresolved(abundance_paths, cache):
+    unresolved = {}
+    for path in abundance_paths:
+        if path.startswith("Unclassified;"):
+            continue
+        parts = path.split(";")
+        for depth in range(1, len(parts) + 1):
+            subpath = ";".join(parts[:depth])
+            if cache.get(subpath, 0) <= 0:
+                unresolved[subpath] = {"path": subpath, "depth": depth, "name": parts[depth - 1]}
+    return [unresolved[path] for path in sorted(unresolved)]
+
+
+def query_exact_scientific_name(name, email, api_key, attempts=3):
+    params = {
+        "db": "taxonomy",
+        "term": f'"{name}"[Scientific Name]',
+        "retmode": "json",
+        "tool": TOOL_NAME,
+        "email": email,
+    }
+    if api_key:
+        params["api_key"] = api_key
+    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
+    last_error = None
+    for attempt in range(attempts):
+        try:
+            request = urllib.request.Request(url, headers={"User-Agent": f"{TOOL_NAME}/1.0"})
+            with urllib.request.urlopen(request, timeout=15) as response:
+                result = json.loads(response.read().decode("utf-8"))
+            ids = result.get("esearchresult", {}).get("idlist", [])
+            return (min(map(int, ids)) if ids else 0), None
+        except Exception as exc:
+            last_error = f"{type(exc).__name__}: {exc}"
+            if attempt + 1 < attempts:
+                time.sleep(2 ** attempt)
+    return 0, last_error
+
+
+def write_unresolved_tsv(path, unresolved):
+    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
+    with open(path, "w", encoding="utf-8", newline="") as handle:
+        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
+        writer.writerow(["Depth", "NodeName", "TaxonPath"])
+        for item in unresolved:
+            writer.writerow([item["depth"], item["name"], item["path"]])
+
+
+def write_conflicts_tsv(path, conflicts):
+    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
+    with open(path, "w", encoding="utf-8", newline="") as handle:
+        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
+        writer.writerow(["Lineage", "WinnerTaxID", "TaxIDCountsJSON"])
+        for item in conflicts:
+            writer.writerow([item["lineage"], item["winner_taxid"], json.dumps(item["counts"], sort_keys=True)])
+

 def main():
     parser = argparse.ArgumentParser(description="Resolve NCBI TaxIDs from wf-16s data")
-    parser.add_argument("--abundance", required=True, help="Path to abundance_table_species.tsv")
-    parser.add_argument("--assignments", default=None, help="Path to read assignments TSV")
-    parser.add_argument("--cache", required=True, help="Path to taxonomy_cache.json")
+    parser.add_argument("--abundance", required=True)
+    parser.add_argument("--tax-column", default="tax")
+    parser.add_argument("--assignments", action="append", default=[])
+    parser.add_argument("--cache", required=True, help="Read-only source cache in cache_only mode")
+    parser.add_argument("--resolved-cache", required=True, help="Run-local resolved cache output")
     parser.add_argument("--mode", choices=["cache_only", "refresh"], default="cache_only")
-    parser.add_argument("--email", default=None, help="NCBI Contact Email")
-    parser.add_argument("--api-key", default=None, help="NCBI API Key")
+    parser.add_argument("--email-env", default="NCBI_EMAIL")
+    parser.add_argument("--api-key-env", default="NCBI_API_KEY")
     parser.add_argument("--unresolved-policy", choices=["warn", "error"], default="warn")
-    parser.add_argument("--unresolved-tsv", default=None, help="Path to write unresolved TaxIDs")
-    parser.add_argument("--provenance", default=None, help="Path to write resolver provenance JSON")
+    parser.add_argument("--unresolved-tsv", required=True)
+    parser.add_argument("--conflicts-tsv", required=True)
+    parser.add_argument("--provenance", required=True)
     args = parser.parse_args()

-    print(f"[taxonomy] Running in '{args.mode}' mode...")
+    try:
+        if not os.path.exists(args.abundance):
+            raise ValueError(f"Abundance file not found: {args.abundance}")
+        abundance_paths = read_abundance_paths(args.abundance, args.tax_column)
+        cache_sha_before = compute_sha256(args.cache)
+        cache = load_cache(args.cache)
+        assignment_map, conflicts = read_assignment_taxids(args.assignments)

-    # 1. Load abundance table lineages
-    if not os.path.exists(args.abundance):
-        print(f"[taxonomy] ERROR: Abundance file not found: {args.abundance}", file=sys.stderr)
-        sys.exit(1)
+        for path in abundance_paths:
+            parts = path.split(";")
+            if len(parts) == 8 and cache.get(path, 0) <= 0:
+                assignment_taxid = assignment_map.get(normalize_abundance_path_to_7(path), 0)
+                if assignment_taxid > 0:
+                    cache[path] = assignment_taxid

-    ab_paths = []
-    with open(args.abundance, "r", encoding="utf-8") as f:
-        header = f.readline().strip().split("\t")
-        for line in f:
-            parts = line.strip().split("\t")
-            if parts and parts[0]:
-                ab_paths.append(parts[0])
+        unresolved = find_unresolved(abundance_paths, cache)
+        query_failures = []
+        cache_updated = False

-    print(f"[taxonomy] Loaded {len(ab_paths)} lineages from abundance table.")
+        if args.mode == "refresh" and unresolved:
+            email = os.environ.get(args.email_env, "").strip()
+            if not email:
+                raise ValueError(f"Environment variable {args.email_env!r} is required for refresh mode.")
+            api_key = os.environ.get(args.api_key_env, "").strip() or None
+            delay = 0.12 if api_key else 0.35
+            name_results = {}
+            for item in unresolved:
+                name = item["name"]
+                if name not in name_results:
+                    name_results[name] = query_exact_scientific_name(name, email, api_key)
+                    time.sleep(delay)
+                taxid, error = name_results[name]
+                if error:
+                    query_failures.append({"path": item["path"], "name": name, "error": error})
+                elif taxid > 0:
+                    cache[item["path"]] = taxid

-    # 2. Parse assignments for normalized lineage -> TaxID mapping
-    lineage_to_taxids = {}
-    if args.assignments and os.path.exists(args.assignments):
-        print(f"[taxonomy] Parsing assignments file: {args.assignments}")
-        with open(args.assignments, "r", encoding="utf-8") as f:
-            for line in f:
-                parts = line.strip().split("\t")
-                if len(parts) >= 5:
-                    taxid = parts[2].strip()
-                    lineage = parts[4].strip()
-                    if taxid and taxid != "0" and lineage:
-                        if lineage not in lineage_to_taxids:
-                            lineage_to_taxids[lineage] = []
-                        lineage_to_taxids[lineage].append(taxid)
-
-    # Resolve conflicts via majority vote
-    resolved_lineage_taxid = {}
-    conflict_records = []
-    for lin, tids in lineage_to_taxids.items():
-        counts = Counter(tids)
-        most_common = counts.most_common()
-        winner = most_common[0][0]
-        resolved_lineage_taxid[lin] = int(winner)
-        if len(most_common) > 1:
-            conflict_records.append({
-                "lineage": lin,
-                "winner_taxid": int(winner),
-                "counts": dict(counts)
-            })
+            unresolved = find_unresolved(abundance_paths, cache)
+            if not query_failures:
+                atomic_write_json(args.cache, cache)
+                cache_updated = True

-    if conflict_records:
-        print(f"[taxonomy] Resolved {len(conflict_records)} multi-TaxID conflicts via majority voting.")
+        atomic_write_json(args.resolved_cache, cache)
+        write_unresolved_tsv(args.unresolved_tsv, unresolved)
+        write_conflicts_tsv(args.conflicts_tsv, conflicts)

-    # 3. Load existing cache
-    cache = {}
-    if os.path.exists(args.cache):
-        try:
-            with open(args.cache, "r", encoding="utf-8") as f:
-                cache = json.load(f)
-            print(f"[taxonomy] Loaded existing cache with {len(cache)} entries.")
-        except Exception as e:
-            print(f"[taxonomy] WARNING: Could not read cache: {e}", file=sys.stderr)
-
-    # 4. In offline mode, identify unresolved nodes
-    unresolved_nodes = []
-    for path in ab_paths:
-        if path.startswith("Unclassified"):
-            continue
-        parts = path.split(";")
-        for j in range(1, len(parts) + 1):
-            subpath = ";".join(parts[:j])
-            current_taxid = cache.get(subpath, 0)
-            if current_taxid == 0:
-                # Try matching from assignments if leaf
-                if j == 8:
-                    norm_7 = normalize_abundance_path_to_7(subpath)
-                    asgn_taxid = resolved_lineage_taxid.get(norm_7)
-                    if asgn_taxid:
-                        cache[subpath] = asgn_taxid
-                        current_taxid = asgn_taxid
-
-            if current_taxid == 0:
-                unresolved_nodes.append({
-                    "path": subpath,
-                    "depth": j,
-                    "name": parts[j - 1]
-                })
-
-    # De-duplicate unresolved nodes
-    seen_unresolved = set()
-    dedup_unresolved = []
-    for u in unresolved_nodes:
-        if u["path"] not in seen_unresolved:
-            seen_unresolved.add(u["path"])
-            dedup_unresolved.append(u)
-
-    print(f"[taxonomy] Unresolved nodes in cache: {len(dedup_unresolved)}")
-
-    # 5. Refresh mode (if requested and authorized)
-    if args.mode == "refresh" and dedup_unresolved:
-        email = args.email or os.environ.get("NCBI_EMAIL")
-        if not email:
-            print("[taxonomy] ERROR: NCBI_EMAIL environment variable or --email is required for refresh mode.", file=sys.stderr)
-            sys.exit(1)
-        api_key = args.api_key or os.environ.get("NCBI_API_KEY")
-
-        # Collect unique species names to query
-        needed_names = list({u["name"] for u in dedup_unresolved if u["depth"] == 8})
-        print(f"[taxonomy] Querying NCBI for {len(needed_names)} species...")
-
-        # Batch querying NCBI esearch / efetch with bounded retries
-        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
-        delay = 0.35 if api_key else 1.0
-
-        for sp_name in needed_names:
-            params = {
-                "db": "taxonomy",
-                "term": sp_name,
-                "retmode": "json",
-                "tool": TOOL_NAME,
-                "email": email
-            }
-            if api_key:
-                params["api_key"] = api_key
-
-            url = f"{base_url}?{urllib.parse.urlencode(params)}"
-            retries = 3
-            found_taxid = 0
-            while retries > 0:
-                try:
-                    req = urllib.request.Request(url, headers={"User-Agent": f"{TOOL_NAME}/1.0"})
-                    with urllib.request.urlopen(req, timeout=10) as resp:
-                        res_json = json.loads(resp.read().decode("utf-8"))
-                        id_list = res_json.get("esearchresult", {}).get("idlist", [])
-                        if id_list:
-                            found_taxid = int(id_list[0])
-                    break
-                except Exception as ex:
-                    retries -= 1
-                    time.sleep(2.0)
-
-            if found_taxid > 0:
-                for path in ab_paths:
-                    if path.endswith(f";{sp_name}"):
-                        cache[path] = found_taxid
-
-            time.sleep(delay)
-
-        # Atomic write back to cache
-        cache_dir = os.path.dirname(os.path.abspath(args.cache))
-        with tempfile.NamedTemporaryFile("w", dir=cache_dir, delete=False, encoding="utf-8") as tmp_f:
-            json.dump(cache, tmp_f, indent=2)
-            tmp_path = tmp_f.name
-        os.replace(tmp_path, args.cache)
-        print(f"[taxonomy] Successfully updated cache atomically: {args.cache}")
-
-    # 6. Write unresolved TSV
-    if args.unresolved_tsv:
-        os.makedirs(os.path.dirname(os.path.abspath(args.unresolved_tsv)), exist_ok=True)
-        with open(args.unresolved_tsv, "w", encoding="utf-8") as f:
-            f.write("Depth\tNodeName\tTaxonPath\n")
-            for u in dedup_unresolved:
-                f.write(f"{u['depth']}\t{u['name']}\t{u['path']}\n")
-        print(f"[taxonomy] Wrote unresolved nodes to: {args.unresolved_tsv}")
-
-    # 7. Write provenance sidecar
-    if args.provenance:
-        os.makedirs(os.path.dirname(os.path.abspath(args.provenance)), exist_ok=True)
-        prov = {
+        provenance = {
             "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             "mode": args.mode,
             "tool": TOOL_NAME,
             "abundance_sha256": compute_sha256(args.abundance),
-            "assignments_sha256": compute_sha256(args.assignments) if args.assignments else None,
-            "cache_sha256": compute_sha256(args.cache),
-            "total_lineages": len(ab_paths),
-            "unresolved_count": len(dedup_unresolved),
-            "conflicts_count": len(conflict_records),
-            "conflicts": conflict_records[:50]
+            "assignments_sha256": {path: compute_sha256(path) for path in args.assignments},
+            "source_cache_sha256_before": cache_sha_before,
+            "source_cache_sha256_after": compute_sha256(args.cache),
+            "resolved_cache_sha256": compute_sha256(args.resolved_cache),
+            "source_cache_updated": cache_updated,
+            "total_lineages": len(abundance_paths),
+            "unresolved_count": len(unresolved),
+            "conflicts_count": len(conflicts),
+            "conflicts": conflicts,
+            "query_failures": query_failures,
         }
-        with open(args.provenance, "w", encoding="utf-8") as f:
-            json.dump(prov, f, indent=2)
+        atomic_write_json(args.provenance, provenance)

-    # 8. Check policy
-    if args.unresolved_policy == "error" and len(dedup_unresolved) > 0:
-        print(f"[taxonomy] ERROR: unresolved_policy is 'error' and {len(dedup_unresolved)} nodes unresolved.", file=sys.stderr)
-        sys.exit(1)
+        if query_failures:
+            print(f"[taxonomy] ERROR: {len(query_failures)} NCBI query failure(s); source cache preserved.", file=sys.stderr)
+            return 1
+        if args.unresolved_policy == "error" and unresolved:
+            print(f"[taxonomy] ERROR: {len(unresolved)} taxonomy nodes remain unresolved.", file=sys.stderr)
+            return 1
+        print(f"[taxonomy] Resolution complete: {len(unresolved)} unresolved, {len(conflicts)} conflicts.")
+        return 0
+    except (OSError, ValueError, json.JSONDecodeError) as exc:
+        print(f"[taxonomy] ERROR: {exc}", file=sys.stderr)
+        return 1

-    print("[taxonomy] Resolution complete.")
-    sys.exit(0)

 if __name__ == "__main__":
-    main()
+    sys.exit(main())
diff --git a/config.example.yml b/config.example.yml
index 41ba8e2..6744006 100644
--- a/config.example.yml
+++ b/config.example.yml
@@ -63,6 +63,7 @@ composition:
 beta:
   distances: ["bray", "jaccard"]
   permutations: 999
+  minimum_count: 1           # presence threshold used for binary Jaccard
   strata_column: null
   resampling:
     enabled: false
diff --git a/config.yml b/config.yml
index a7182e6..14d2058 100644
--- a/config.yml
+++ b/config.yml
@@ -36,6 +36,7 @@ composition:
 beta:
   distances: ["bray", "jaccard"]
   permutations: 999
+  minimum_count: 1           # presence threshold used for binary Jaccard
   strata_column: null
   resampling:
     enabled: false
diff --git a/tests/test_ncbi_taxonomy.py b/tests/test_ncbi_taxonomy.py
new file mode 100644
index 0000000..e9a7d88
--- /dev/null
+++ b/tests/test_ncbi_taxonomy.py
@@ -0,0 +1,77 @@
+import hashlib
+import importlib.util
+import json
+import os
+import pathlib
+import tempfile
+import unittest
+from unittest import mock
+
+ROOT = pathlib.Path(__file__).resolve().parents[1]
+MODULE_PATH = ROOT / "analysis" / "utils" / "ncbi_taxonomy.py"
+SPEC = importlib.util.spec_from_file_location("ncbi_taxonomy", MODULE_PATH)
+taxonomy = importlib.util.module_from_spec(SPEC)
+SPEC.loader.exec_module(taxonomy)
+
+
+class TaxonomyResolverTests(unittest.TestCase):
+    def setUp(self):
+        self.tempdir = tempfile.TemporaryDirectory()
+        self.work = pathlib.Path(self.tempdir.name)
+        self.lineage = "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis"
+        self.abundance = self.work / "abundance.tsv"
+        self.abundance.write_text(
+            "tax\tS1\ttotal\n"
+            "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t1\t1\n"
+            f"{self.lineage}\t2\t2\n",
+            encoding="utf-8",
+        )
+        self.assignment = self.work / "assignment.tsv"
+        self.assignment.write_text(
+            "C\tread1\t1423\t0|1500\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis\n"
+            "C\tread2\t1423\t1501\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis\n"
+            "U\tread3\t0\t1490\tUnclassified\n",
+            encoding="utf-8",
+        )
+        parts = self.lineage.split(";")
+        self.cache = self.work / "cache.json"
+        cache_payload = {";".join(parts[:depth]): depth for depth in range(1, 8)}
+        cache_payload[self.lineage] = 0
+        self.cache.write_text(json.dumps(cache_payload), encoding="utf-8")
+
+    def tearDown(self):
+        self.tempdir.cleanup()
+
+    def args(self, mode="cache_only"):
+        return [
+            "ncbi_taxonomy.py", "--abundance", str(self.abundance),
+            "--assignments", str(self.assignment), "--cache", str(self.cache),
+            "--resolved-cache", str(self.work / "resolved.json"), "--mode", mode,
+            "--unresolved-policy", "warn", "--unresolved-tsv", str(self.work / "unresolved.tsv"),
+            "--conflicts-tsv", str(self.work / "conflicts.tsv"),
+            "--provenance", str(self.work / "provenance.json"),
+        ]
+
+    def test_cache_only_uses_assignments_without_mutating_source_cache(self):
+        before = hashlib.sha256(self.cache.read_bytes()).hexdigest()
+        with mock.patch("sys.argv", self.args()):
+            self.assertEqual(taxonomy.main(), 0)
+        after = hashlib.sha256(self.cache.read_bytes()).hexdigest()
+        self.assertEqual(before, after)
+        resolved = json.loads((self.work / "resolved.json").read_text(encoding="utf-8"))
+        self.assertEqual(resolved[self.lineage], 1423)
+        self.assertEqual((self.work / "unresolved.tsv").read_text(encoding="utf-8").count("\n"), 1)
+
+    def test_refresh_failure_preserves_source_cache_and_returns_nonzero(self):
+        before = hashlib.sha256(self.cache.read_bytes()).hexdigest()
+        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
+             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
+             mock.patch.object(taxonomy, "query_exact_scientific_name", return_value=(0, "simulated failure")), \
+             mock.patch("sys.argv", self.args(mode="refresh")):
+            self.assertEqual(taxonomy.main(), 1)
+        after = hashlib.sha256(self.cache.read_bytes()).hexdigest()
+        self.assertEqual(before, after)
+
+
+if __name__ == "__main__":
+    unittest.main()
diff --git a/tests/testthat/helper-fixtures.R b/tests/testthat/helper-fixtures.R
index e7212ec..4782c17 100644
--- a/tests/testthat/helper-fixtures.R
+++ b/tests/testthat/helper-fixtures.R
@@ -6,13 +6,13 @@ create_temp_abundance <- function(dir, n_species = 10, sample_names = c("Sample1
                                   include_total = TRUE, unclass_reads = 100) {
   ranks_template <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp%02d"
   lineages <- vapply(seq_len(n_species), function(i) sprintf(ranks_template, i), character(1))
-
+
   unclass_lineage <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
   all_lineages <- c(unclass_lineage, lineages)
-
+
   set.seed(123)
   count_df <- data.frame(tax = all_lineages, stringsAsFactors = FALSE)
-
+
   for (idx in seq_along(sample_names)) {
     s <- sample_names[idx]
     base_lambda <- if (idx %% 2 == 1) 40 else 80
@@ -23,11 +23,11 @@ create_temp_abundance <- function(dir, n_species = 10, sample_names = c("Sample1
     counts <- c(unclass_reads, sp_counts)
     count_df[[s]] <- counts
   }
-
+
   if (include_total) {
     count_df$total <- rowSums(as.matrix(count_df[, sample_names, drop = FALSE]))
   }
-
+
   file_path <- file.path(dir, "synthetic_abundance.tsv")
   write.table(count_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
   file_path
@@ -36,7 +36,7 @@ create_temp_abundance <- function(dir, n_species = 10, sample_names = c("Sample1
 create_temp_assignments <- function(dir, sample_id = "Sample1", n_classified = 50, n_unclassified = 10) {
   total <- n_classified + n_unclassified
   read_ids <- sprintf("read_%05d", seq_len(total))
-
+
   status <- c(rep("C", n_classified), rep("U", n_unclassified))
   taxids <- c(rep(1386L, n_classified), rep(0L, n_unclassified))
   len_fields <- c(
@@ -47,7 +47,7 @@ create_temp_assignments <- function(dir, sample_id = "Sample1", n_classified = 5
     rep("Bacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus cereus", n_classified),
     rep("Unclassified", n_unclassified)
   )
-
+
   df <- data.frame(
     status = status,
     read_id = read_ids,
@@ -56,7 +56,7 @@ create_temp_assignments <- function(dir, sample_id = "Sample1", n_classified = 5
     lineage = lineages,
     stringsAsFactors = FALSE
   )
-
+
   file_path <- file.path(dir, sprintf("%s_assignments.tsv", sample_id))
   write.table(df, file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
   file_path
diff --git a/tests/testthat/test-alpha-regression.R b/tests/testthat/test-alpha-regression.R
index 397284c..62dfbcc 100644
--- a/tests/testthat/test-alpha-regression.R
+++ b/tests/testthat/test-alpha-regression.R
@@ -9,40 +9,40 @@ source(file.path("..", "..", "analysis", "utils", "metrics.R"))
 test_that("Real Ambar Ayunda fixture reproduces exact Section 2.3 alpha regression metrics", {
   ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
   skip_if_not(file.exists(ab_path), "Real abundance table not found")
-
+
   res <- read_abundance_table(ab_path)
   sample_col <- res$samples[1]
   counts <- res$count_matrix[, sample_col]
   class_counts <- counts[-res$unclass_index]
-
+
   alpha_df <- calc_alpha_indices(class_counts)
   vals <- setNames(alpha_df$Value, alpha_df$Metric)
-
+
   # Target values from Section 2.3:
   # Observed species richness: 1,836
   expect_equal(vals[["Observed species richness (S)"]], 1836)
-
+
   # Chao1: 2,983.851
   expect_equal(round(vals[["Chao1 (estimated richness)"]], 3), 2983.851, tolerance = 0.001)
-
+
   # Shannon: 4.603
   expect_equal(round(vals[["Shannon (H)"]], 3), 4.603, tolerance = 0.001)
-
+
   # Effective species number: 99.814
   expect_equal(round(vals[["Effective number of species (e^H)"]], 3), 99.814, tolerance = 0.001)
-
+
   # Simpson: 0.957
   expect_equal(round(vals[["Simpson's D (1-sum p^2)"]], 3), 0.957, tolerance = 0.001)
-
+
   # Inverse Simpson: 23.089
   expect_equal(round(vals[["Inverse Simpson"]], 3), 23.089, tolerance = 0.001)
-
+
   # Pielou evenness: 0.613
   expect_equal(round(vals[["Pielou's evenness (J)"]], 3), 0.613, tolerance = 0.001)
-
+
   # Fisher alpha: 334.543
   expect_equal(round(vals[["Fisher's alpha"]], 3), 334.543, tolerance = 0.001)
-
+
   # Berger-Parker dominance: 0.133
   expect_equal(round(vals[["Berger-Parker dominance"]], 3), 0.133, tolerance = 0.001)
 })
@@ -50,10 +50,10 @@ test_that("Real Ambar Ayunda fixture reproduces exact Section 2.3 alpha regressi
 test_that("Detected classified taxa count by rank matches Section 2.3 targets", {
   ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
   skip_if_not(file.exists(ab_path), "Real abundance table not found")
-
+
   res <- read_abundance_table(ab_path)
   tax_df <- res$taxonomy[-res$unclass_index, ]
-
+
   # Expected: 31 phyla, 69 classes, 128 orders, 281 families, 867 genera, 1836 species
   expect_equal(length(unique(tax_df$phylum)), 31)
   expect_equal(length(unique(tax_df$class)), 69)
diff --git a/tests/testthat/test-cohort.R b/tests/testthat/test-cohort.R
index eafbfc6..3c7a9e0 100644
--- a/tests/testthat/test-cohort.R
+++ b/tests/testthat/test-cohort.R
@@ -6,6 +6,7 @@ source(file.path("..", "..", "analysis", "utils", "config.R"))
 source(file.path("..", "..", "analysis", "utils", "io.R"))
 source(file.path("..", "..", "analysis", "utils", "metrics.R"))
 source(file.path("..", "..", "analysis", "utils", "plotting.R"))
+source(file.path("..", "..", "analysis", "02_alpha_diversity.R"))
 source(file.path("..", "..", "analysis", "03_beta_diversity.R"))
 source(file.path("..", "..", "analysis", "05_ordination.R"))
 source(file.path("..", "..", "analysis", "06_shared_taxa.R"))
@@ -13,29 +14,30 @@ source(file.path("..", "..", "analysis", "06_shared_taxa.R"))
 test_that("Single-sample context skips beta, ordination, and shared_taxa gracefully", {
   tmp <- tempdir()
   ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = c("S1"))
-
+
   cfg <- get_default_config()
   cfg$input$abundance_table <- ab_file
   cfg$output$base_dir <- file.path(tmp, "out_single")
   cfg$output$dirs <- list(
     beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"),
+    alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"),
     ordination = file.path(cfg$output$base_dir, "05_Ordination"),
     shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa")
   )
-
+
   context <- build_context(cfg)
   expect_equal(context$mode, "single")
-
+
   # Run beta
   res_beta <- run_beta(context)
   expect_equal(res_beta$status, "skipped")
   expect_true(file.exists(file.path(cfg$output$dirs$beta, "beta_diversity_skipped.tsv")))
-
+
   # Run ordination
   res_ord <- run_ordination(context)
   expect_equal(res_ord$status, "skipped")
   expect_true(file.exists(file.path(cfg$output$dirs$ordination, "ordination_skipped.tsv")))
-
+
   # Run shared taxa
   res_shared <- run_shared_taxa(context)
   expect_equal(res_shared$status, "skipped")
@@ -48,7 +50,7 @@ test_that("Synthetic cohort (2 groups x 3 replicates) passes cohort gates and pa
   ab_file <- create_temp_abundance(tmp, n_species = 20, sample_names = sample_names)
   meta_file <- create_temp_metadata(tmp, sample_names = sample_names,
                                     groups = c("Control", "Control", "Control", "Treated", "Treated", "Treated"))
-
+
   cfg <- get_default_config()
   cfg$mode <- "cohort"
   cfg$input$abundance_table <- ab_file
@@ -59,24 +61,32 @@ test_that("Synthetic cohort (2 groups x 3 replicates) passes cohort gates and pa
     ordination = file.path(cfg$output$base_dir, "05_Ordination"),
     shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa")
   )
-
+
   context <- build_context(cfg)
   expect_equal(context$mode, "cohort")
   expect_equal(length(context$samples), 6)
-
+
+  # Run cohort alpha, including the metadata-aware group test path.
+  res_alpha <- run_alpha(context)
+  expect_equal(res_alpha$status, "completed")
+  alpha_out <- read.delim(file.path(cfg$output$dirs$alpha, "alpha_diversity.tsv"), check.names = FALSE)
+  expect_equal(alpha_out$SampleID, sample_names)
+  expect_true("Group" %in% names(alpha_out))
+  expect_false(any(c("Group.x", "Group.y") %in% names(alpha_out)))
+
   # Run beta
   res_beta <- run_beta(context)
   expect_equal(res_beta$status, "completed")
   expect_true(file.exists(file.path(cfg$output$dirs$beta, "permanova.tsv")))
   expect_true(file.exists(file.path(cfg$output$dirs$beta, "betadisper.tsv")))
   expect_true(file.exists(file.path(cfg$output$dirs$beta, "pcoa_scores_bray.tsv")))
-
+
   # Run ordination
   res_ord <- run_ordination(context)
   expect_equal(res_ord$status, "completed")
   expect_true(file.exists(file.path(cfg$output$dirs$ordination, "pca_scores.tsv")))
   expect_true(file.exists(file.path(cfg$output$dirs$ordination, "pca_variance.tsv")))
-
+
   # Run shared taxa
   res_shared <- run_shared_taxa(context)
   expect_equal(res_shared$status, "completed")
@@ -85,23 +95,35 @@ test_that("Synthetic cohort (2 groups x 3 replicates) passes cohort gates and pa
   expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "group_prevalence.tsv")))
 })

+test_that("Forced cohort mode rejects a one-sample table", {
+  tmp <- tempfile("forced_cohort_")
+  dir.create(tmp)
+  ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = "Only1")
+  meta_file <- create_temp_metadata(tmp, sample_names = "Only1", groups = "Control")
+  cfg <- get_default_config()
+  cfg$mode <- "cohort"
+  cfg$input$abundance_table <- ab_file
+  cfg$input$metadata <- meta_file
+  expect_error(build_context(cfg), "at least 2")
+})
+
 test_that("Under-replicated cohort skips PERMANOVA with explicit reason", {
   tmp <- tempdir()
   # 2 groups with only 1 sample each
   sample_names <- c("Ctrl1", "Trt1")
   ab_file <- create_temp_abundance(tmp, n_species = 10, sample_names = sample_names)
   meta_file <- create_temp_metadata(tmp, sample_names = sample_names, groups = c("Control", "Treated"))
-
+
   cfg <- get_default_config()
   cfg$mode <- "cohort"
   cfg$input$abundance_table <- ab_file
   cfg$input$metadata <- meta_file
   cfg$output$base_dir <- file.path(tmp, "out_underrep")
   cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
-
+
   context <- build_context(cfg)
   res_beta <- run_beta(context)
-
+
   perm_tsv <- file.path(cfg$output$dirs$beta, "permanova.tsv")
   expect_true(file.exists(perm_tsv))
   perm_df <- read.delim(perm_tsv)
diff --git a/tests/testthat/test-config.R b/tests/testthat/test-config.R
index ed9a9c2..f61b90c 100644
--- a/tests/testthat/test-config.R
+++ b/tests/testthat/test-config.R
@@ -7,17 +7,17 @@ source(file.path("..", "..", "analysis", "utils", "config.R"))
 test_that("load_config loads default config.yml and resolves relative paths to config dir", {
   config_path <- file.path("..", "..", "config.yml")
   expect_true(file.exists(config_path))
-
+
   cfg <- load_config(config_path)
-
+
   expect_equal(cfg$schema_version, 1L)
   expect_equal(cfg$mode, "auto")
   expect_equal(cfg$seed, 42L)
-
+
   # Base output directory must be derived
   expect_true(is.character(cfg$output$base_dir))
   expect_true(nzchar(cfg$output$base_dir))
-
+
   # All 7 module directories must be present and derived from base_dir
   expect_true("qc" %in% names(cfg$output$dirs))
   expect_true("alpha" %in% names(cfg$output$dirs))
@@ -26,7 +26,7 @@ test_that("load_config loads default config.yml and resolves relative paths to c
   expect_true("ordination" %in% names(cfg$output$dirs))
   expect_true("shared_taxa" %in% names(cfg$output$dirs))
   expect_true("kreport" %in% names(cfg$output$dirs))
-
+
   expect_equal(cfg$output$dirs$qc, file.path(cfg$output$base_dir, "01_QC"))
   expect_equal(cfg$output$dirs$alpha, file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
   expect_equal(cfg$output$dirs$beta, file.path(cfg$output$base_dir, "03_Beta_Diversity"))
@@ -39,9 +39,9 @@ test_that("load_config loads default config.yml and resolves relative paths to c
 test_that("CLI --output-dir overrides base_dir and all derived paths", {
   config_path <- file.path("..", "..", "config.yml")
   override_dir <- file.path(tempdir(), "test_override_output")
-
+
   cfg <- load_config(config_path, cli_opts = list(output_dir = override_dir))
-
+
   expect_equal(normalizePath(cfg$output$base_dir, winslash = "/", mustWork = FALSE),
                normalizePath(override_dir, winslash = "/", mustWork = FALSE))
   expect_equal(cfg$output$dirs$qc, file.path(cfg$output$base_dir, "01_QC"))
@@ -51,3 +51,8 @@ test_that("CLI --output-dir overrides base_dir and all derived paths", {
 test_that("load_config errors on missing config file", {
   expect_error(load_config("non_existent_config.yml"), "Configuration file not found")
 })
+
+test_that("unknown config keys fail closed", {
+  bad <- get_default_config()
+  expect_error(merge_config(bad, list(alhpa = list())), "Unknown configuration key.*alhpa")
+})
diff --git a/tests/testthat/test-io.R b/tests/testthat/test-io.R
index 97239bf..ed84e6c 100644
--- a/tests/testthat/test-io.R
+++ b/tests/testthat/test-io.R
@@ -7,7 +7,7 @@ source(file.path("..", "..", "analysis", "utils", "io.R"))
 test_that("Synthetic abundance table parses and validates correctly", {
   tmp <- tempdir()
   ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = c("S1", "S2"))
-
+
   res <- read_abundance_table(ab_file)
   expect_equal(length(res$samples), 2)
   expect_equal(res$samples, c("S1", "S2"))
@@ -18,12 +18,12 @@ test_that("Synthetic abundance table parses and validates correctly", {
 test_that("Abundance table errors if total column does not equal row sum", {
   tmp <- tempdir()
   ab_file <- create_temp_abundance(tmp, n_species = 3, sample_names = c("S1"))
-
+
   # Corrupt total column
   df <- read.delim(ab_file, check.names = FALSE)
   df$total[1] <- df$total[1] + 999
   write.table(df, ab_file, sep = "\t", row.names = FALSE, quote = FALSE)
-
+
   expect_error(read_abundance_table(ab_file), "does not equal sample row sums")
 })

@@ -37,14 +37,14 @@ test_that("Abundance table errors if rank count is not 8", {
     total = c(10, 50)
   )
   write.table(bad_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
-
+
   expect_error(read_abundance_table(file_path), "expected 8 ranks, found 6")
 })

 test_that("Assignments parser handles pipe and plain lengths and reconciles counts", {
   tmp <- tempdir()
   asgn_file <- create_temp_assignments(tmp, sample_id = "S1", n_classified = 40, n_unclassified = 10)
-
+
   reads <- read_assignments_file(
     asgn_file,
     sample_id = "S1",
@@ -52,34 +52,42 @@ test_that("Assignments parser handles pipe and plain lengths and reconciles coun
     expected_classified = 40,
     expected_unclassified = 10
   )
-
+
   expect_equal(nrow(reads), 50)
   expect_equal(sum(reads$effective_classified), 40)
   expect_true(all(reads$read_length >= 1200 & reads$read_length <= 1600))
   expect_true(all(is.integer(reads$read_length)))
 })

+test_that("Assignments parser rejects non-five-field rows with a physical line number", {
+  tmp <- tempfile("bad_assignment_")
+  dir.create(tmp)
+  bad_file <- file.path(tmp, "bad.tsv")
+  writeLines(c("C\tread_1\t123\t0|1500\tBacteria|Example", "U\tread_2\t0\t1400"), bad_file)
+  expect_error(read_assignments_file(bad_file, "S1"), "4 fields at line 2; expected exactly 5")
+})
+
 test_that("Real Ambar Ayunda fixture satisfies all Section 2.2 invariants", {
   ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
   asgn_path <- file.path("..", "..", "output_AAy", "reads_assignments",
                          "AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv")
-
+
   skip_if_not(file.exists(ab_path), "Real abundance table not found")
   skip_if_not(file.exists(asgn_path), "Real assignments file not found")
-
+
   ab_res <- read_abundance_table(ab_path)
   expect_equal(nrow(ab_res$count_matrix), 1837)
   expect_equal(ab_res$samples, "AmbarAyunda_minimap2_16S")
-
+
   sample_col <- ab_res$samples[1]
   total_reads <- sum(ab_res$count_matrix[, sample_col])
   unclass_reads <- ab_res$count_matrix[ab_res$unclass_index, sample_col]
   classified_reads <- total_reads - unclass_reads
-
+
   expect_equal(total_reads, 114056)
   expect_equal(classified_reads, 80556)
   expect_equal(unclass_reads, 33500)
-
+
   # Check assignments
   reads <- read_assignments_file(
     asgn_path,
@@ -88,7 +96,7 @@ test_that("Real Ambar Ayunda fixture satisfies all Section 2.2 invariants", {
     expected_classified = 80556,
     expected_unclassified = 33500
   )
-
+
   expect_equal(nrow(reads), 114056)
   expect_equal(sum(reads$status == "C"), 89809)
   expect_equal(sum(reads$status == "U"), 24247)
@@ -99,14 +107,22 @@ test_that("Real Ambar Ayunda fixture satisfies all Section 2.2 invariants", {
 test_that("Metadata validation aligns samples and detects discrepancies", {
   tmp <- tempdir()
   meta_file <- create_temp_metadata(tmp, sample_names = c("S2", "S1"), groups = c("GroupB", "GroupA"))
-
+
   meta_aligned <- read_metadata_table(meta_file, selected_samples = c("S1", "S2"))
   expect_equal(meta_aligned$SampleID, c("S1", "S2"))
   expect_equal(meta_aligned$Group, c("GroupA", "GroupB"))
-
+
   # Test missing sample
   expect_error(
     read_metadata_table(meta_file, selected_samples = c("S1", "S2", "S3")),
     "Missing from metadata: S3"
   )
 })
+
+test_that("Metadata rejects empty groups", {
+  tmp <- tempfile("bad_metadata_")
+  dir.create(tmp)
+  meta_file <- file.path(tmp, "metadata.tsv")
+  writeLines(c("SampleID\tGroup", "S1\t"), meta_file)
+  expect_error(read_metadata_table(meta_file, "S1"), "empty Group")
+})
diff --git a/tests/testthat/test-kreport.R b/tests/testthat/test-kreport.R
index cac04f2..0fa9aa8 100644
--- a/tests/testthat/test-kreport.R
+++ b/tests/testthat/test-kreport.R
@@ -10,12 +10,12 @@ source(file.path("..", "..", "analysis", "07_kreport_pavian.R"))
 test_that("kreport tree builder uses standard rank codes D, K, P, C, O, F, G, S", {
   ranks_template <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_subtilis"
   unclass_lineage <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
-
+
   lineages <- c(unclass_lineage, ranks_template)
   counts <- c(10, 50)
-
+
   nodes <- build_kreport_tree(lineages, counts)
-
+
   # Standard rank codes: D, K, P, C, O, F, G, S
   expect_equal(nodes$rank_code, c("D", "K", "P", "C", "O", "F", "G", "S"))
   expect_false("D1" %in% nodes$rank_code) # Defect fixed!
@@ -26,16 +26,16 @@ test_that("kreport tree validates clade arithmetic and total read invariants", {
   lin1 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp1"
   lin2 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp2"
   uncl <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
-
+
   lineages <- c(uncl, lin1, lin2)
   counts <- c(20, 30, 50)
   total <- sum(counts)
-
+
   nodes <- build_kreport_tree(lineages, counts)
-
+
   # Should validate without error
   expect_silent(validate_kreport_tree(nodes, total_reads = total, uncl_reads = 20))
-
+
   # Check failure when counts corrupted
   expect_error(
     validate_kreport_tree(nodes, total_reads = total + 10, uncl_reads = 20),
@@ -48,37 +48,38 @@ test_that("Real Ambar Ayunda fixture builds valid .kreport and runs offline", {
   cache_path <- file.path("..", "..", "output_AAy", "taxonomy_cache.json")
   asgn_path <- file.path("..", "..", "output_AAy", "reads_assignments",
                          "AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv")
-
+
   skip_if_not(file.exists(ab_path), "Abundance table not found")
   skip_if_not(file.exists(cache_path), "Cache not found")
-
+
   tmp <- tempdir()
   out_dir <- file.path(tmp, "test_kreport_out")
-
+
   cfg <- get_default_config()
   cfg$config_dir <- normalizePath(file.path("..", ".."), winslash = "/")
+  cfg$pipeline_root <- cfg$config_dir
   cfg$input$abundance_table <- ab_path
   cfg$taxonomy$cache <- cache_path
   cfg$input$assignments <- list(AmbarAyunda_minimap2_16S = asgn_path)
   cfg$output$base_dir <- out_dir
   cfg$output$dirs <- list(kreport = file.path(out_dir, "07_Kreport"))
-
+
   context <- build_context(cfg)
-
+
   res <- run_kreport(context)
   expect_equal(res$status, "completed")
-
+
   kreport_file <- file.path(cfg$output$dirs$kreport, "AmbarAyunda_minimap2_16S.kreport")
   expect_true(file.exists(kreport_file))
-
+
   lines <- readLines(kreport_file)
   expect_gt(length(lines), 100)
-
+
   # Line 1: unclassified
   expect_match(lines[1], "^[0-9.]+\t33500\t33500\tU\t0\tunclassified")
   # Line 2: root
   expect_match(lines[2], "^[0-9.]+\t80556\t0\tR\t1\troot")
-
+
   # Check kingdom row uses K
   expect_true(any(grepl("\tK\t", lines)))
   expect_false(any(grepl("\tD1\t", lines)))
~~~~
