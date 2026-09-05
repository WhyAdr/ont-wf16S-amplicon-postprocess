# `ont-wf16S-amplicon-postprocess` v0.1.0 release-readiness audit

**Repository:** `WhyAdr/ont-wf16S-amplicon-postprocess`
**Audited head:** `4810fa67ad8cd10608c6e3d910741343203e351b` (`main`, 2026-09-04)
**Blueprint:** `ont-wf16s-amplicon-postprocess-xpand1.md`
**Recommended first tag:** `v0.1.0` after all P0 gates below are green
**Verdict:** **NO-GO for tagging at the audited head; close to release after one focused hardening pass.**

## 1. Executive assessment

Phase 1 has substantially improved the repository. The shared validated context, strict abundance and assignment parsing, explicit single/cohort gates, classified-only denominators, deterministic seeds, offline taxonomy mode, run-local resolved cache, atomic refresh behavior, Kraken arithmetic checks, manifest generation, and modular architecture are all real implementations rather than placeholders. The audited fixture still reconciles to **114,056 total = 80,556 classified + 33,500 unclassified**, and the resolver independently reproduces **46 unresolved nodes and 26 assignment-lineage conflicts** in cache-only mode.

However, the repository is not yet release-ready for four decisive reasons:

1. **The required CI gate is red.** Both workflow runs to date failed. At head `4810fa6`, package installation failed through `fs -> pkgload -> testthat`; all later CI steps were skipped, including R parsing, the 105 R assertions, Python tests, real-fixture validation, and whitespace checking. See [GitHub Actions run 33828788325](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/actions/runs/33828788325).
2. **The latest Windows-path fix contradicts the R contract.** Commit `4810fa6` removed argument quoting on the claim that `system2()` quotes arguments internally. [R's own `system2()` documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/help/system2.html) says arguments containing spaces or special characters must be quoted. The current call can split paths and expose shell metacharacters. Use a vector-safe process API and test a real path containing spaces on Windows.
3. **Explicit network opt-in is not actually enforced.** YAML can set `taxonomy.network_mode: refresh`, allowing network access without `--refresh-taxonomy`, contrary to the blueprint's “only CLI action that permits NCBI network access” rule.
4. **The release surface and reproducibility story are incomplete.** `config.example.yml` combines a one-sample abundance file with six unmatched synthetic metadata rows; the README says an included cohort fixture exists when it is created only at test runtime; the manifest has no pipeline version or Python version; per-node taxonomy provenance is absent; and the repository has no license, changelog, citation metadata, or environment lock.

The right final pass is narrow: repair CI, replace the fragile subprocess call, enforce offline-by-default at the boundary, close the missing regression tests, make the example/documentation truthful, and add minimal release metadata. Do not expand Phase 1 into differential abundance, UniFrac, functional prediction, or a redesign of the statistical modules.

## 2. Evidence collected

| Check | Result |
|---|---|
| Head commit | `4810fa6` |
| Tags | None |
| Git working tree after clone | Clean |
| Current GitHub Actions run | **Failed**: run `33828788325` |
| Failed CI step | `Install R dependencies` |
| CI consequence | All parse/test/validation steps skipped |
| Local Python compilation | Passed |
| Local Python unit tests | **2/2 passed** |
| Local cache-only resolver on tracked fixture | Passed; 46 unresolved, 26 conflicts |
| Source-cache mutation in that run | None; SHA-256 remained `a2bdbeca...17e2cf` |
| Local R suite in this audit environment | Not run because `Rscript` is unavailable |
| Claimed developer-side R result in commit message | 105 assertions passed; not independently confirmed by CI |
| License detected by GitHub | None |
| Repository description/topics | Empty |

The developer-side 105-assertion result is useful evidence, but it cannot substitute for a clean-room CI run. A release tag should point only at a commit whose required status checks are green.

## 3. Blueprint gate matrix

| Blueprint gate | Status | Assessment |
|---|---|---|
| Gate A — contracts and scaffold | **Partial** | Strict parsers and real-fixture invariants exist. Tests do not yet prove outside-CWD runner invocation, process-level mutation-free validation, all sample-ID rejection cases, or CLI exit behavior. |
| Gate B — single-sample modules | **Partial/strong** | Alpha regression targets and module outputs are tested. No byte-stability assertion exists for seeded resampling, and no clean full-run integration test validates every promised artifact from the runner. |
| Gate C — cohort modules | **Partial/strong** | A synthetic 2×3 cohort exercises alpha, beta, ordination, and shared taxa. Metadata alignment and under-replication skipping are covered. Input-order invariance, deterministic keyed results, and stronger PERMANOVA/betadisper content checks are missing. |
| Gate D — kreport/taxonomy | **Partial/strong** | Offline resolution, rank `K`, arithmetic, conflicts, and atomic refresh failure are implemented. Exact 46/26 regression counts are not asserted, per-node source provenance is missing, and subprocess path handling is currently unsafe. |
| Gate E — integration/docs | **Fail** | CI is red; only Ubuntu is configured; no successful full runner execution occurs in CI; examples/docs contain mismatches; release metadata is absent. |

## 4. Findings by priority

### P0-1 — CI is red and therefore the release gate is objectively closed

The workflow uses a bespoke `install.packages(missing_pkgs)` call. On the clean Ubuntu runner, a failed `fs` installation caused `pkgload` and then `testthat` to fail. This is exactly the class of dependency-resolution problem CI is supposed to expose. The workflow currently tests only `ubuntu-latest`, despite a specific Windows compatibility claim. The maintained [`r-lib/actions`](https://github.com/r-lib/actions) suite provides dependency-aware setup actions intended for CI.

**Required outcome:** install dependencies with a dependency-aware action, pin the interpreter versions for the release branch, run Ubuntu and Windows, and require both jobs to pass.

### P0-2 — `system2()` paths with spaces are still broken or unsafe

`analysis/07_kreport_pavian.R` passes raw arguments to `system2()`. The R manual explicitly requires quoting arguments containing spaces or special characters. The newest commit's rationale is therefore incorrect. Re-applying `shQuote()` can remain platform-fragile and easy to misuse, so the highest-confidence fix is `processx::run()`, whose API accepts an argument vector without shell re-parsing.

**Impact:** the taxonomy/kreport module may fail on `D:/...` or other paths with spaces, and unquoted metacharacters create avoidable command-injection risk from user-controlled paths.

### P0-3 — Network access can be enabled without the required CLI flag

`load_config()` accepts `taxonomy.network_mode: refresh` from YAML. The blueprint states that `--refresh-taxonomy` is the only action that permits NCBI access. A reproducible offline run must not become online merely because a copied YAML says `refresh`.

**Required outcome:** reject YAML `refresh`; allow it only when the parsed CLI flag is true. Record that opt-in in the manifest.

### P0-4 — Critical release behavior is not exercised end-to-end

The suite tests many functions, but the release gate is process-level. Missing checks include:

- full runner success from a non-repository current directory;
- a full output path containing spaces;
- `--validate-only` creating no output path;
- non-zero exit on an invalid module;
- manifest status and required fields after a full run;
- exact `46` unresolved and `26` conflict counts;
- stable, complete required output inventory;
- Windows execution in CI.

### P1-1 — `config.example.yml` is internally inconsistent

It points to the tracked one-sample Ambar Ayunda abundance table but also sets `metadata.example.tsv`, which contains six unrelated sample IDs. `mode: auto` resolves to single, but metadata is still validated exactly and must fail with missing/extra IDs. An example copied by a user should not be guaranteed to fail.

### P1-2 — README overclaims an included cohort fixture

The README says “The included cohort fixture is synthetic.” No cohort data fixture is tracked; tests create temporary synthetic data at runtime. The scientific caveat is good, but the noun is inaccurate.

### P1-3 — Dependency ownership is duplicated and runtime includes a test-only package

The authoritative package vector appears independently in `analysis/utils/dependencies.R` and `analysis/install_packages.R`. This violates the blueprint and invites drift. `testthat` is also required by every pipeline execution even though it is not a runtime dependency.

### P1-4 — Manifest provenance is incomplete

The manifest records R but not Python, and contains no pipeline semantic version, git commit, or original command. `taxonomy_resolution.tsv` labels nodes only `Resolved`/`Unresolved`; it cannot distinguish source-cache, assignment-derived, or NCBI-refreshed TaxIDs. Global resolver provenance does not satisfy evidence-level source provenance.

### P1-5 — Environment reproduction is documented but not delivered

The output records installed versions and `sessionInfo()`, which is provenance. It is not an environment lock. A clean machine installing “latest CRAN” can produce a different environment next week. This is acceptable for the original Phase 1 definition only if the release is explicitly described as provenance-captured but not bitwise environment-locked. Prefer an `renv.lock` generated from the final green reference environment.

### P1-6 — Public release metadata is absent

There is no `LICENSE`, `CHANGELOG.md`, `CITATION.cff`, or machine-readable package/pipeline version. Without a license, public visibility does not grant reuse rights. Confirm the desired license; MIT is a reasonable default for this small analysis pipeline.

### P2 — Repository presentation and fixture hygiene

- GitHub description and topics are blank.
- Tracked `params.json` exposes machine-specific absolute paths and an EPI2ME instance identifier. These are not credentials, but they are needless environment disclosure in a public release fixture. Redact them only if doing so does not invalidate the intended regression contract, and document the redaction.
- Large legacy outputs, plots, and the Pavian HTML remain tracked. The blueprint intentionally preserved them, so do not remove them in this patch. Label them clearly as historical/example artifacts and consider a later `examples/` or release-asset migration.
- `--overwrite` permits stale artifacts from an older run to coexist with new outputs. The manifest lists current outputs, which partly mitigates this, but users should be told to use a fresh output root for publication runs.

## 5. Exact patch plan

Apply the patches in order. The snippets are intended as implementation instructions; re-run formatting and adjust hunk line numbers if intervening commits move the code.

### Patch 1 — Use a vector-safe subprocess API

```diff
diff --git a/analysis/utils/dependencies.R b/analysis/utils/dependencies.R
--- a/analysis/utils/dependencies.R
+++ b/analysis/utils/dependencies.R
@@
-REQUIRED_PACKAGES <- c(
+RUNTIME_PACKAGES <- c(
   "yaml",
   "optparse",
   "dplyr",
@@
   "UpSetR",
   "ggrepel",
-  "digest",
-  "testthat"
+  "digest",
+  "processx"
 )

-get_required_packages <- function() {
-  REQUIRED_PACKAGES
+TEST_PACKAGES <- c("testthat")
+REQUIRED_PACKAGES <- unique(c(RUNTIME_PACKAGES, TEST_PACKAGES))
+
+get_required_packages <- function(include_tests = FALSE) {
+  if (isTRUE(include_tests)) REQUIRED_PACKAGES else RUNTIME_PACKAGES
 }

-check_dependencies <- function(pkgs = REQUIRED_PACKAGES) {
+check_dependencies <- function(pkgs = RUNTIME_PACKAGES) {
@@
-get_dependency_versions <- function(pkgs = REQUIRED_PACKAGES) {
+get_dependency_versions <- function(pkgs = RUNTIME_PACKAGES) {
```

```diff
diff --git a/analysis/07_kreport_pavian.R b/analysis/07_kreport_pavian.R
--- a/analysis/07_kreport_pavian.R
+++ b/analysis/07_kreport_pavian.R
@@
 suppressMessages({
   library(jsonlite)
   library(dplyr)
+  library(processx)
 })
@@
-  # Run Python script
-  res_code <- system2(python_cmd, args = cmd_args)
-  if (res_code != 0) {
-    stop(sprintf("NCBI taxonomy resolver failed with exit status %d", res_code), call. = FALSE)
+  # processx passes a true argument vector on Windows and Unix; do not shell-quote.
+  resolver <- processx::run(
+    command = python_cmd,
+    args = cmd_args,
+    echo = TRUE,
+    error_on_status = FALSE
+  )
+  if (!identical(resolver$status, 0L)) {
+    stop(sprintf("NCBI taxonomy resolver failed with exit status %d", resolver$status),
+         call. = FALSE)
   }
```

Add a regression test that uses spaces in every relevant path:

```diff
diff --git a/tests/testthat/test-kreport.R b/tests/testthat/test-kreport.R
--- a/tests/testthat/test-kreport.R
+++ b/tests/testthat/test-kreport.R
@@
+test_that("kreport resolver handles input and output paths containing spaces", {
+  root <- tempfile("wf16s path with spaces ")
+  dir.create(root, recursive = TRUE, showWarnings = FALSE)
+
+  lineage <- paste(c("Bacteria", "Bacillati", "Bacillota", "Bacilli",
+                     "Bacillales", "Bacillaceae", "Bacillus",
+                     "Bacillus subtilis"), collapse = ";")
+  abundance <- file.path(root, "abundance table.tsv")
+  writeLines(c(
+    "tax\tS1\ttotal",
+    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t1\t1",
+    paste(lineage, "2", "2", sep = "\t")
+  ), abundance)
+
+  assignments <- file.path(root, "read assignments.tsv")
+  writeLines(c(
+    "C\tread1\t1423\t0|1500\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis",
+    "C\tread2\t1423\t1501\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis",
+    "U\tread3\t0\t1490\tUnclassified"
+  ), assignments)
+
+  parts <- strsplit(lineage, ";", fixed = TRUE)[[1]]
+  cache <- setNames(as.list(seq_len(7L)),
+                    vapply(seq_len(7L), function(i) paste(parts[seq_len(i)], collapse = ";"),
+                           character(1)))
+  cache[[lineage]] <- 0L
+  cache_file <- file.path(root, "taxonomy cache.json")
+  jsonlite::write_json(cache, cache_file, auto_unbox = TRUE)
+
+  cfg <- get_default_config()
+  cfg$pipeline_root <- normalizePath(file.path("..", ".."), winslash = "/")
+  cfg$input$abundance_table <- abundance
+  cfg$input$assignments <- list(S1 = assignments)
+  cfg$taxonomy$cache <- cache_file
+  cfg$output$base_dir <- file.path(root, "output directory")
+  cfg$output$dirs <- list(kreport = file.path(cfg$output$base_dir, "07_Kreport"))
+  cfg$cli <- list(modules = "kreport", validate_only = FALSE)
+
+  context <- build_context(cfg)
+  expect_equal(run_kreport(context)$status, "completed")
+  expect_true(file.exists(file.path(cfg$output$dirs$kreport, "S1.kreport")))
+})
```

### Patch 2 — Make network refresh CLI-only

```diff
diff --git a/analysis/utils/config.R b/analysis/utils/config.R
--- a/analysis/utils/config.R
+++ b/analysis/utils/config.R
@@
   default_cfg <- get_default_config()
   cfg <- merge_config(default_cfg, raw_yaml)
+  cli_refresh <- isTRUE(cli_opts$refresh_taxonomy) ||
+    isTRUE(cli_opts[["refresh-taxonomy"]])
+
+  if (identical(cfg$taxonomy$network_mode, "refresh") && !cli_refresh) {
+    stop(paste(
+      "YAML cannot enable taxonomy refresh.",
+      "Keep taxonomy.network_mode: cache_only and pass --refresh-taxonomy explicitly."
+    ), call. = FALSE)
+  }

   # CLI overrides
@@
-  if (isTRUE(cli_opts$refresh_taxonomy) || isTRUE(cli_opts[["refresh-taxonomy"]])) {
+  if (cli_refresh) {
     cfg$taxonomy$network_mode <- "refresh"
   }
@@
     overwrite = isTRUE(cli_opts$overwrite),
+    refresh_taxonomy = cli_refresh,
```

```diff
diff --git a/config.example.yml b/config.example.yml
--- a/config.example.yml
+++ b/config.example.yml
@@
-  metadata: "metadata.example.tsv"
+  # Leave null for the tracked one-sample demonstration. For a real cohort,
+  # provide metadata whose SampleID set exactly matches the abundance columns.
+  metadata: null
@@
-  network_mode: "cache_only"  # 'cache_only' (offline default) or 'refresh'
+  # YAML is deliberately offline-only. Use --refresh-taxonomy for an online run.
+  network_mode: "cache_only"
```

```diff
diff --git a/tests/testthat/test-config.R b/tests/testthat/test-config.R
--- a/tests/testthat/test-config.R
+++ b/tests/testthat/test-config.R
@@
+test_that("taxonomy refresh requires explicit CLI opt-in", {
+  tmp <- tempfile("refresh_config_")
+  dir.create(tmp)
+  cfg <- get_default_config()
+  cfg$taxonomy$network_mode <- "refresh"
+  path <- file.path(tmp, "config.yml")
+  yaml::write_yaml(cfg, path)
+
+  expect_error(load_config(path), "YAML cannot enable taxonomy refresh")
+  resolved <- load_config(path, cli_opts = list(refresh_taxonomy = TRUE))
+  expect_equal(resolved$taxonomy$network_mode, "refresh")
+  expect_true(resolved$cli$refresh_taxonomy)
+})
```

Also change the blueprint's example comment if the blueprint is retained as live documentation; it should not advertise `refresh` as a YAML choice.

### Patch 3 — Centralize dependency installation and repair CI

Replace `analysis/install_packages.R`'s duplicate `REQUIRED_PACKAGES` declaration with the canonical file:

```diff
diff --git a/analysis/install_packages.R b/analysis/install_packages.R
--- a/analysis/install_packages.R
+++ b/analysis/install_packages.R
@@
 args <- commandArgs(trailingOnly = TRUE)
 do_install <- "--install" %in% args

-REQUIRED_PACKAGES <- c(
-  "yaml",
-  "optparse",
-  "dplyr",
-  "tidyr",
-  "stringr",
-  "ggplot2",
-  "scales",
-  "vegan",
-  "RColorBrewer",
-  "jsonlite",
-  "pheatmap",
-  "UpSetR",
-  "ggrepel",
-  "digest",
-  "testthat"
-)
+all_args <- commandArgs(trailingOnly = FALSE)
+file_arg <- grep("^--file=", all_args, value = TRUE)
+script_dir <- if (length(file_arg)) {
+  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
+} else {
+  normalizePath("analysis", winslash = "/", mustWork = TRUE)
+}
+source(file.path(script_dir, "utils", "dependencies.R"))
+REQUIRED_PACKAGES <- get_required_packages(include_tests = TRUE)
```

Replace `.github/workflows/ci.yml` with:

```yaml
name: validation

on:
  push:
  pull_request:

permissions:
  contents: read

jobs:
  test:
    name: ${{ matrix.os }} / R ${{ matrix.r }} / Python ${{ matrix.python }}
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, windows-latest]
        r: ['4.5']
        python: ['3.12']

    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python }}

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.r }}
          use-public-rspm: true

      - name: Install R dependencies with pak
        uses: r-lib/actions/setup-r-dependencies@v2
        with:
          packages: |
            any::yaml
            any::optparse
            any::dplyr
            any::tidyr
            any::stringr
            any::ggplot2
            any::scales
            any::vegan
            any::RColorBrewer
            any::jsonlite
            any::pheatmap
            any::UpSetR
            any::ggrepel
            any::digest
            any::processx
            any::testthat
          cache-version: 1

      - name: Parse R and Python sources
        shell: bash
        run: |
          Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
          python -m compileall -q analysis tests

      - name: Run unit and regression tests
        shell: bash
        run: |
          Rscript tests/testthat.R
          python -m unittest -v tests/test_ncbi_taxonomy.py

      - name: Validate without filesystem mutation
        shell: bash
        run: |
          test ! -e "ci validation output"
          Rscript analysis/00_run_pipeline.R --config config.yml --output-dir "ci validation output" --validate-only
          test ! -e "ci validation output"

      - name: Run full release integration
        shell: bash
        run: |
          Rscript analysis/00_run_pipeline.R --config config.yml --output-dir "ci integration output"
          Rscript tests/verify_release_run.R "ci integration output"

      - name: Check committed whitespace
        if: runner.os == 'Linux'
        shell: bash
        run: git show --check --oneline --no-renames HEAD
```

Add `tests/verify_release_run.R`:

```r
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: verify_release_run.R OUTPUT_DIR")
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
manifest <- jsonlite::fromJSON(file.path(root, "run_manifest.json"),
                               simplifyVector = FALSE)
stopifnot(identical(manifest$run_status, "completed"))
stopifnot(identical(manifest$mode, "single"))
stopifnot(identical(as.integer(manifest$taxonomy$unresolved_count), 46L))

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
  "07_Kreport/unresolved_taxids.tsv",
  "07_Kreport/taxonomy_conflicts.tsv",
  "07_Kreport/taxonomy_provenance.json"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Missing release outputs: ", paste(missing, collapse = ", "))

conflicts <- read.delim(file.path(root, "07_Kreport/taxonomy_conflicts.tsv"),
                        check.names = FALSE)
stopifnot(nrow(conflicts) == 26L)
cat("Release integration verification passed.\n")
```

### Patch 4 — Assert the known taxonomy incompleteness and conflicts

Extend the real kreport regression test:

```diff
diff --git a/tests/testthat/test-kreport.R b/tests/testthat/test-kreport.R
--- a/tests/testthat/test-kreport.R
+++ b/tests/testthat/test-kreport.R
@@
   expect_true(any(grepl("\tK\t", lines)))
   expect_false(any(grepl("\tD1\t", lines)))
+
+  unresolved <- read.delim(file.path(cfg$output$dirs$kreport,
+                                     "unresolved_taxids.tsv"),
+                           check.names = FALSE)
+  conflicts <- read.delim(file.path(cfg$output$dirs$kreport,
+                                    "taxonomy_conflicts.tsv"),
+                          check.names = FALSE)
+  expect_equal(nrow(unresolved), 46L)
+  expect_equal(nrow(conflicts), 26L)
 })
```

Add tests for deterministic tie-breaking and ambiguous NCBI results to `tests/test_ncbi_taxonomy.py`. Do not silently select `min(idlist)` when an exact-name query returns multiple TaxIDs; retain the node as unresolved and record ambiguity. This avoids manufacturing a confident but contextually wrong parent TaxID.

### Patch 5 — Add per-node taxonomy source provenance

Extend the Python resolver with a required `--resolution-sources-tsv` output. Track sources as follows:

- positive entry present at load: `source_cache`;
- leaf filled from exact seven-rank assignment lineage: `assignment`;
- node filled by a successful unambiguous NCBI refresh: `ncbi_refresh`;
- no positive ID: `unresolved`.

Minimum TSV contract:

```text
TaxonPath<TAB>TaxID<TAB>ResolutionSource
```

Then pass the new path from `07_kreport_pavian.R`, require it alongside the other resolver outputs, and join it into `taxonomy_resolution.tsv` so every row has `ResolutionSource`. Add source counts to `taxonomy_provenance.json` and `run_manifest.json`.

This patch is intentionally specified as a contract rather than a partial code hunk: it touches resolver state across load, assignment merge, refresh, TSV writing, R joining, manifest assembly, and both test suites. Gemini should implement it atomically and add tests for all four source labels.

### Patch 6 — Add semantic version and complete manifest provenance

Create `VERSION`:

```text
0.1.0
```

At runner startup, read and validate it:

```diff
diff --git a/analysis/00_run_pipeline.R b/analysis/00_run_pipeline.R
--- a/analysis/00_run_pipeline.R
+++ b/analysis/00_run_pipeline.R
@@
 repo_root <- normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE)
+version_file <- file.path(repo_root, "VERSION")
+pipeline_version <- trimws(readLines(version_file, n = 1L, warn = FALSE))
+if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", pipeline_version)) {
+  stop("VERSION must contain a semantic version such as 0.1.0", call. = FALSE)
+}
@@
 manifest <- list(
   pipeline = "ont-wf16s-postprocess",
+  pipeline_version = pipeline_version,
   schema_version = cfg$schema_version,
@@
-  interpreter = list(r = R.version.string, platform = R.version$platform),
+  interpreter = list(
+    r = R.version.string,
+    platform = R.version$platform,
+    python = tryCatch(
+      processx::run(python_cmd, "--version", error_on_status = FALSE)$stdout,
+      error = function(e) NA_character_
+    )
+  ),
```

Do not reference the `python_cmd` local variable from `run_kreport()`; factor Python discovery into a shared utility such as `find_python()` and call it from both the runner and kreport module. Also record:

- `git_commit` when the repository is a Git checkout, else `null`;
- `commandArgs(trailingOnly = FALSE)` after redacting nothing secret (the API key value is never passed as an argument);
- `refresh_taxonomy` explicit boolean;
- taxonomy resolution-source counts.

Add manifest schema assertions to `tests/verify_release_run.R`.

### Patch 7 — Make documentation and examples truthful

```diff
diff --git a/README.md b/README.md
--- a/README.md
+++ b/README.md
@@
-The included cohort fixture is synthetic and verifies software behavior only. It is not biological validation of cohort statistics or species-level classification accuracy. For reproducible publication work, record the wf-16S/database versions and create an `renv.lock` from the R environment used for the final analysis.
+The test suite creates a temporary synthetic cohort at runtime to verify software
+behavior only; no biologically validated cohort dataset is included. Every run
+records input hashes, package versions, and session information. These capture
+provenance but do not by themselves recreate an environment. For publication work,
+use the committed `renv.lock` from the tagged release (or create and archive a
+project-specific lock from the final analysis environment).
```

Add sections covering:

1. **Supported upstream contract:** explicitly state that Phase 1 is validated against the tracked minimap2 schema and wf-16S metadata represented by the fixture, not arbitrary historical/future wf-16S schemas.
2. **Quickstart distinction:** label `config.yml` as the tracked Ambar Ayunda regression demonstration and `config.example.yml` as the copy-and-edit template.
3. **Cohort setup:** show a short mapping where abundance sample columns, metadata `SampleID`, and assignment-map keys are identical.
4. **Taxonomy incompleteness:** state that the reference offline run intentionally reports 46 unresolved nodes and 26 conflicting lineage-to-TaxID mappings; these are surfaced, not hidden.
5. **Fresh output roots:** recommend a new output directory for publication runs to prevent stale files from older executions.
6. **Release scope:** state that cohort behavior is synthetically tested but not biologically validated.

### Patch 8 — Add release metadata

After the repository owner confirms the license, add:

- `LICENSE` — recommended MIT, copyright `2026 WhyAdr`;
- `CHANGELOG.md` — `0.1.0` entry summarizing Phase 1 contracts, modules, offline taxonomy, provenance, and known limitation that cohort validation is synthetic;
- `CITATION.cff` — title, version `0.1.0`, release date, repository URL, preferred citation author identity, and license;
- repository description and topics such as `nanopore`, `16s-rrna`, `amplicon-sequencing`, `microbiome`, `r`, and `bioinformatics`.

Do not invent an ORCID or legal author name. Use the owner-confirmed identity.

### Patch 9 — Generate and validate an environment lock

Once the fixed Ubuntu job is green, generate the lock from that exact environment:

```bash
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")'
Rscript -e 'renv::init(bare = TRUE); renv::snapshot(type = "all", prompt = FALSE)'
```

Inspect `renv.lock` to ensure all runtime packages plus `testthat` and `processx` are present. Then change CI from free-floating package installation to `r-lib/actions/setup-renv@v2`. Keep one scheduled or manual “latest dependency compatibility” workflow separate from the release-blocking locked workflow.

If locking is deliberately deferred, document that `v0.1.0` records provenance but does not guarantee environment reconstruction. Do not use the word “fully reproducible.”

## 6. Additional missing tests to add before tagging

These tests are smaller than the failure classes they prevent:

1. **Config/CWD:** invoke the runner from a different working directory with an absolute config path.
2. **CLI precedence:** prove CLI output and refresh flags override YAML, and YAML cannot independently enable refresh.
3. **Mutation-free validation:** snapshot the non-existent output path before and after `--validate-only`.
4. **Exit status:** invalid module, malformed assignments, invalid metadata, simulated module failure, and simulated refresh failure must all return non-zero.
5. **Keep-going:** inject one module failure; prove later modules run, manifest records the failure, and final status is non-zero.
6. **Seed determinism:** run alpha resampling twice in the same environment and compare the TSV bytes or SHA-256.
7. **Cohort order invariance:** permute abundance columns and metadata rows; compare keyed numeric tables after sorting by identifiers.
8. **Sample IDs:** reject separators, control characters, `.`/`..`, empty IDs, and post-sanitization collisions.
9. **Taxonomy:** assert 46 unresolved, 26 conflicts, deterministic tie-break, no source-cache mutation, and no URL opener call in cache-only mode.
10. **Manifest contract:** validate semantic version, interpreters, hashes, module statuses, skip reasons, outputs, and taxonomy source counts.

## 7. Final release gate

Tag `v0.1.0` only when all boxes are true:

- [ ] Ubuntu R/Python CI is green from a clean cache.
- [ ] Windows R/Python CI is green, including a path containing spaces.
- [ ] The R manual contradiction is resolved by a vector-safe subprocess API.
- [ ] The 105 existing assertions plus the new release tests pass.
- [ ] Full tracked-fixture execution succeeds through the runner.
- [ ] `--validate-only` is proven mutation-free.
- [ ] Cache-only is proven network-free; refresh requires the explicit CLI flag.
- [ ] Tracked fixture yields 114,056 total, 80,556 classified, 33,500 unclassified.
- [ ] Offline taxonomy yields exactly 46 unresolved nodes and 26 conflicts.
- [ ] Every requested module is completed or skipped with an explicit reason.
- [ ] Manifest contains pipeline version, R/Python versions, input hashes, module results, and taxonomy provenance.
- [ ] `config.example.yml` validates for its stated scenario.
- [ ] README claims match actual outputs and limitations.
- [ ] License choice and citation identity are confirmed.
- [ ] `LICENSE`, `CHANGELOG.md`, `CITATION.cff`, and `VERSION` are present.
- [ ] `git show --check HEAD` is clean.
- [ ] The tag target commit has no failing required checks.

Suggested commands after implementation:

```bash
Rscript analysis/install_packages.R --install
Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
python -m compileall -q analysis tests
Rscript tests/testthat.R
python -m unittest -v tests/test_ncbi_taxonomy.py
Rscript analysis/00_run_pipeline.R --config config.yml --output-dir "release candidate output"
Rscript tests/verify_release_run.R "release candidate output"
git show --check --oneline --no-renames HEAD
git status --short
```

On Windows PowerShell, run the same `Rscript`/`python` commands with an output path containing spaces. Require the GitHub Windows job as the authoritative clean-machine check.

## 8. Recommended tag and release notes posture

Use `v0.1.0`, not `v1.0.0`. The software is a credible first public Phase 1 release once the gates above pass, but the cohort path is still validated only with synthetic data and the taxonomy export intentionally exposes unresolved/conflicting mappings. Those are normal early-release boundaries, not embarrassments; hiding them would be the embarrassing part.

The release notes should state:

- strict support for the documented wf-16S abundance and minimap2 assignment contracts;
- validated single-sample regression dataset and expected counts;
- synthetic-only cohort software validation;
- offline-by-default taxonomy resolution;
- expected unresolved/conflict diagnostics;
- no automatic Pavian HTML export;
- environment-lock status;
- known scope exclusions from the blueprint.

## 9. Bottom line

The core Phase 1 implementation is substantially sound and worth preserving. The main risk is not that the analytical architecture is fundamentally wrong; it is that the repository currently advertises release confidence that its clean CI has never demonstrated, while the most recent Windows-path patch likely reintroduced the exact class of failure it meant to fix. One disciplined P0 pass should convert this from “promising implementation with local evidence” into a defensible `v0.1.0` release candidate.
