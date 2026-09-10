# `ont-wf16S-amplicon-postprocess` v0.3.1 patch plan

## Decision

**Do not tag v0.3.1 from the audited tree as-is.** The current implementation at
`5fbd583192190322e895a4aa48c63d7a98bef2fa` is healthy on its tracked fixtures,
but the version bump is not internally coherent and several valid edge inputs can
still crash or be reported incorrectly.

The target state is **GO after the required patches below and four green CI jobs**
(Ubuntu/Windows core plus Ubuntu/Windows FAPROTAX). Keep this release a narrow
correctness and hardening patch: do not redesign denominators, taxonomy resolution,
FAPROTAX semantics, or the Krona count model.

## Audit basis

- Repository: `WhyAdr/ont-wf16S-amplicon-postprocess`
- Audited commit: `5fbd583192190322e895a4aa48c63d7a98bef2fa`
- Declared version: `0.3.0`
- Tree state during review: clean
- Current upstream CI run for this SHA: all four jobs passed:
  <https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/actions/runs/34090914956>
- Locally verified in the audit environment:
  - `python3 -m unittest -v tests/test_ncbi_taxonomy.py`: 6/6 passed
  - `python3 -m compileall -q analysis tests`: passed
  - `python3 tests/check_committed_whitespace.py`: passed
- R was not installed in the audit container; the existing R baseline is therefore
  supported by the green commit-pinned GitHub Actions run, while every new R change
  below must pass the required CI matrix before tagging.

## Validation of Z's findings

| ID | Verdict | Evidence and correction |
|---|---|---|
| Z1 release-bump landmine | **Confirmed, P1** | `tests/verify_release_run.R:13` requires literal `0.3.0`. Active version strings also occur in `VERSION`, `CITATION.cff`, **five** README locations (not four: lines 5, 51, 77, 307, 332), and two runtime messages in `analysis/utils/io.R`. Add a new changelog entry; do not rewrite the historical `0.3.0` heading. Make release verification derive the expected value from `VERSION` so this trap does not recur. |
| Z2 fractional upstream thresholds | **Confirmed, P2** | `read_upstream_params()` accepts any finite numeric `min_len`, `max_len`, and `abundance_threshold`; `01_qc_diagnostics.R:228-229` later passes lengths to `%d`. A value such as `1400.5` therefore survives contract validation and fails during plot construction with an unrelated formatting error. Enforce whole-number semantics at the producer boundary. |
| Z3 empty lineage ranks | **Confirmed, P2** | Only the complete lineage string and its field count are checked. An internal empty or whitespace-only field therefore reaches composition, kreport, Krona, and the Python resolver. Reject empty, whitespace-only, and boundary-whitespace rank cells before constructing `tax_df`. The two tracked abundance fixtures contain no such rows. |
| Z4a derived-seed overflow | **Confirmed, P3** | `seed + strtoi(...)` can use overflowing integer arithmetic and yield `NA`. Derive the sample seed in double precision, modulo the accepted seed range, then cast once. Also cap configured seeds at `.Machine$integer.max`, because the same seed is passed directly to `set.seed()` in other modules. |
| Z4b graphics-device leaks | **Confirmed, P3** | The raw `png()` calls in modules 04 and 06 have no cleanup path when `pheatmap()` or `UpSetR::upset()` errors. Introduce one tested graphics-device wrapper that closes the exact device and removes an incomplete file on error. |
| Z4c bamstats whitespace | **Confirmed, P3** | Discovery checks `nzchar(trimws(sample_name))` but uses the untrimmed value for mapping. A file containing only `"S1 "` is silently treated as belonging to an unknown sample. Reject boundary whitespace explicitly; do not silently trim identifiers. |
| Z4d duplicate modules | **Confirmed, P3** | `--modules qc,qc` executes QC twice and overwrites the first module record. Reject duplicates when parsing the module list. |
| Z4e late module validation | **Confirmed, P3** | `build_context(cfg)` runs before the module registry is checked, so an invalid module can trigger large assignment parsing and unrelated input errors first. Validate module names and the Krona→kreport dependency immediately after config loading and before `build_context()`. |

The three ruled-out claims in Z's note remain ruled out: the displayed `stats[...]`
and `meta[...]` expressions are not defects; `prcomp(rank. = k)` retains the full
`sdev` vector used for variance accounting; and the empty unique-taxa branch in
module 06 retains its table schema.

## Additional findings from the independent review

| ID | Severity | Finding | Required resolution |
|---|---|---|---|
| A1 | **P2** | Numeric-looking `Group` values are type-converted by `read.delim()`. `06_shared_taxa.R` then evaluates `prev_list[[grp]]` and `membership_list[[grp]]`; numeric groups are interpreted as positional indices, so group `0` errors and other values can create unnamed/sparse lists. Numeric-looking `SampleID` values can also lose leading zeroes before matching. | Force `SampleID` and `Group` to character **during import**, use character group keys defensively in module 06, and preserve exact group column names in outputs. |
| A2 | **P2** | `composition.top_n_taxa = 1` is valid, but the cohort heatmap then contains one row. `pheatmap()` defaults to row clustering and reaches `hclust()` with fewer than two objects. The same occurs when the valid input has only one classified taxon. | Set `cluster_rows = nrow(mat_transformed) >= 2L` and the analogous column gate. Add a one-row heatmap regression. |
| A3 | **P2** | NMDS convergence is misreported. vegan's `metaMDSiter()` initializes `converged <- 0`, increments it when the best solution is repeated, and returns that numeric count. `isTRUE(nmds_res$converged)` is therefore false for positive integer counts. | Convert the repetition count explicitly with `> 0`, retain the count in diagnostics, and test `0`, `1`, and `2`. Upstream semantics: <https://github.com/vegandevs/vegan/blob/1895f7c23107e73dab743456df2da5a7349a875d/R/metaMDSiter.R#L10-L10> and <https://github.com/vegandevs/vegan/blob/1895f7c23107e73dab743456df2da5a7349a875d/R/metaMDSiter.R#L137-L168>. |
| A4 | **P3** | Sample filename validation is not fully portable: leading/trailing whitespace is accepted, comparisons are case-sensitive, trailing dots are accepted, and Windows device basenames such as `CON`, `NUL`, `COM1`, or `LPT1.txt` survive `sanitize_filename()`. These can fail or collide on the supported Windows runner. | Reject boundary whitespace, trailing dots, case-folded post-sanitization collisions, and reserved Windows device basenames. Keep the original SampleID in analytical tables and manifests. |
| A5 | **P3** | A syntactically valid YAML scalar at the document root has no names, so `merge_config(default_cfg, raw_yaml)` silently returns the full defaults instead of failing closed. | Require the YAML root to be a named mapping before merging. |
| A6 | **P3** | Beta rarefaction is designed to tolerate failed iterations, but only `cmdscale()` is caught. A degenerate `procrustes()` call escapes the iteration and fails the entire module. | Put rarefaction, PCoA validation, and Procrustes alignment inside the same per-iteration `tryCatch()` and record the iteration as failed. |

## Required patch sequence

Apply the following phases in order. Do not bump release metadata until the code and
tests are ready in the same commit or tightly ordered commit series.

### Patch 1 — Fail fast on the module request

**Files:** `analysis/utils/config.R`, `analysis/00_run_pipeline.R`,
`tests/testthat/test-config.R`, `tests/testthat/test-release-process.R`

Add a single parser in `config.R`:

```r
parse_requested_modules <- function(value) {
  defaults <- c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport")
  if (is.null(value)) return(defaults)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop("'--modules' must contain at least one module name.", call. = FALSE)
  }

  modules <- strsplit(trimws(value), "[,[:space:]]+", perl = TRUE)[[1]]
  duplicates <- unique(modules[duplicated(modules)])
  if (length(duplicates) > 0L) {
    stop(sprintf("Duplicate module name(s): %s", paste(duplicates, collapse = ", ")),
         call. = FALSE)
  }
  modules
}
```

In `load_config()`, compute `requested_modules <- parse_requested_modules(cli_opts$modules)`
before constructing `cfg$cli`, then set `modules = requested_modules`. Remove the
inline `strsplit()` expression.

In `00_run_pipeline.R`, move the module registry, unknown-module check, and
Krona→kreport check to immediately after `cfg$pipeline_root <- repo_root` and before
`build_context(cfg)`. Leave module dependency/FAPROTAX runtime checks after context
construction, but before `--validate-only` and output mutation. The required order is:

```text
parse CLI → load/validate config → validate requested modules/dependencies between
modules → build input context → validate optional package/runtime dependencies →
validate-only or output mutation
```

Add regressions that prove:

1. `load_config(..., cli_opts = list(modules = "qc,qc"))` fails with
   `Duplicate module name`;
2. whitespace-separated and comma-separated unique names preserve user order;
3. an unknown module paired with a deliberately malformed assignment file reports
   `Unknown module` first and does not create the output root;
4. `--krona` without `kreport` remains an early, mutation-free failure.

### Patch 2 — Close producer, lineage, metadata, and identifier contracts

**Files:** `analysis/utils/config.R`, `analysis/utils/io.R`,
`analysis/06_shared_taxa.R`, `tests/testthat/test-io.R`,
`tests/testthat/test-cohort.R`, `README.md`

#### 2.1 Whole-number upstream fields

After the current finite-number loop in `read_upstream_params()`, add:

```r
integer_fields <- c("min_len", "max_len", "abundance_threshold")
for (field in integer_fields) {
  value <- params[[field]]
  if (abs(value - round(value)) > sqrt(.Machine$double.eps) ||
      value > .Machine$integer.max) {
    stop(sprintf(
      "params.json field '%s' must be a whole number no greater than %d.",
      field, .Machine$integer.max
    ), call. = FALSE)
  }
}
```

Keep fractional `min_read_qual`, `min_percent_identity`, and `min_ref_coverage`
legal. The current range checks remain authoritative.

Add a table-driven test that mutates each of `min_len`, `max_len`, and
`abundance_threshold` to a fractional JSON number and expects the contract error.
Also prove `1.0` remains accepted.

#### 2.2 Non-empty, trimmed rank cells

Immediately after the eight-field-count check in `read_abundance_table()` add a
per-cell gate. Report the physical table row (`data row + header`) and rank name:

```r
for (row_index in seq_along(parsed_lineages)) {
  fields <- parsed_lineages[[row_index]]
  bad_rank <- which(!nzchar(trimws(fields)) | fields != trimws(fields))
  if (length(bad_rank) > 0L) {
    rank_index <- bad_rank[1]
    stop(sprintf(
      "Lineage schema violation at row %d: rank '%s' must be non-empty and have no leading/trailing whitespace ('%s').",
      row_index + 1L, RANKS_8[rank_index], lineages[row_index]
    ), call. = FALSE)
  }
}
```

Tests must cover `Bacteria;;Bacillota;...`, `Bacteria; ;Bacillota;...`, and an
otherwise valid field with boundary whitespace. Keep the lower-level Krona
normalization test: defensive normalization of a manually supplied `nodes_df` is
still useful even though production abundance input now fails earlier.

Update the README abundance contract to state that all eight rank fields must be
non-empty after trimming and must not contain boundary whitespace. `Unknown` remains
an allowed explicit label.

#### 2.3 Import metadata identity columns as character

In `read_metadata_table()`, read and validate the header first, then force the two
identity columns to character at import time. Do not read them numerically and cast
afterward, because that destroys a value such as `01`.

```r
header_line <- readLines(path, n = 1L, warn = FALSE)
if (length(header_line) == 0L) {
  stop("Metadata table is empty.", call. = FALSE)
}
header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
if (anyDuplicated(header)) {
  stop("Metadata table contains duplicate column names.", call. = FALSE)
}
required_identity <- c("SampleID", "Group")
missing_identity <- setdiff(required_identity, header)
if (length(missing_identity) > 0L) {
  stop(sprintf("Metadata table must contain column(s): %s.",
               paste(missing_identity, collapse = ", ")), call. = FALSE)
}
col_classes <- rep(NA_character_, length(header))
col_classes[match(required_identity, header)] <- "character"
meta <- read.delim(
  path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  check.names = FALSE, colClasses = col_classes
)
```

Remove the now-duplicated post-read header checks, but retain empty values,
boundary-whitespace, duplicate SampleID, exact set matching, and order alignment.

In `run_shared_taxa()` use explicit character group values and exact column names:

```r
group_values <- as.character(meta$Group)
groups <- unique(group_values)
...
grp_samples <- meta$SampleID[group_values == grp]
...
prev_df <- data.frame(..., prev_mat, check.names = FALSE, stringsAsFactors = FALSE)
mem_df  <- data.frame(..., mem_mat,  check.names = FALSE, stringsAsFactors = FALSE)
```

Add a cohort regression with sample IDs `01`/`02` and groups `0`/`1`; require exact
identity preservation and successful shared-taxa output. Include a group label with
a space or hyphen and verify the emitted group column is not silently renamed.

#### 2.4 Portable SampleIDs and fail-closed YAML root

Extend `validate_sample_ids()`:

```r
if (any(sample_ids != trimws(sample_ids))) {
  stop("Sample ID validation error: Sample IDs must not have leading or trailing whitespace.",
       call. = FALSE)
}

sanitized <- vapply(sample_ids, sanitize_filename, character(1))
if (any(grepl("[.]$", sanitized))) {
  stop("Sample ID validation error: Portable output basenames must not end in a dot.",
       call. = FALSE)
}
if (anyDuplicated(tolower(sanitized))) {
  stop("Sample ID validation error: Sample IDs collide after portable filename normalization.",
       call. = FALSE)
}
windows_base <- toupper(sub("[.].*$", "", sanitized))
reserved <- windows_base %in% c("CON", "PRN", "AUX", "NUL",
                                paste0("COM", 1:9), paste0("LPT", 1:9))
if (any(reserved)) {
  stop(sprintf("Sample ID validation error: '%s' is a reserved Windows device basename.",
               sample_ids[which(reserved)[1]]), call. = FALSE)
}
```

Replace the old case-sensitive post-sanitization duplicate check with the
case-folded check above. Test boundary whitespace, `Sample` versus `sample`, a
trailing dot, `CON`, `NUL.txt`, `COM1`, and `LPT9.tsv`.

Before `merge_config(default_cfg, raw_yaml)` in `load_config()` add:

```r
if (!is.list(raw_yaml) || is.null(names(raw_yaml))) {
  stop("Configuration YAML root must be a named mapping.", call. = FALSE)
}
```

An entirely empty YAML document may either be rejected by the same check or handled
with a separate explicit error; it must not silently select the repository's
dataset-specific defaults. Add scalar-root and empty-document tests.

#### 2.5 Bamstats identifiers must be exact

In `discover_bamstats()`, after the missing/empty check and before selecting
`file_sample`, reject any `sample_name` whose value differs from `trimws(value)`:

```r
if (any(samples_in_file != trimws(samples_in_file))) {
  stop(sprintf("Bamstats file '%s' contains sample_name values with leading/trailing whitespace.",
               path), call. = FALSE)
}
```

Add a discovery regression using one file with `S1 ` and require a non-zero,
specific failure rather than an informational “not mapped” result.

### Patch 3 — Make seed derivation safe without changing ordinary results

**Files:** `analysis/utils/config.R`, `analysis/utils/metrics.R`,
`analysis/02_alpha_diversity.R`, `tests/testthat/test-alpha-regression.R`,
`tests/testthat/test-config.R`

Cap the configured base seed:

```r
assert_scalar_number(cfg$seed, "seed", lower = 0,
                     upper = .Machine$integer.max, integer = TRUE)
```

Add to `metrics.R`:

```r
derive_sample_seed <- function(seed, sample_id) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max ||
      abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
    stop("Base seed must be an integer in the set.seed() range.", call. = FALSE)
  }
  if (!is.character(sample_id) || length(sample_id) != 1L || is.na(sample_id) ||
      !nzchar(sample_id)) {
    stop("Sample ID for seed derivation must be one non-empty string.", call. = FALSE)
  }
  hash_part <- strtoi(substr(
    digest::digest(sample_id, algo = "xxhash32", serialize = FALSE), 1L, 7L
  ), base = 16L)
  as.integer((as.double(seed) + as.double(hash_part)) %%
               as.double(.Machine$integer.max))
}
```

Replace the inline expression in module 02 with
`seed = derive_sample_seed(seed, s)`.

Tests must prove that the maximum accepted base seed returns a finite scalar integer
for multiple sample IDs, repeated calls are identical, `seed + hash` does not
produce `NA`, over-range configuration fails, and the existing ordinary-seed
byte-stability regression remains unchanged.

### Patch 4 — Close graphics devices and support one-row heatmaps

**Files:** `analysis/utils/plotting.R`, `analysis/04_taxa_composition.R`,
`analysis/06_shared_taxa.R`, new `tests/testthat/test-plotting.R`, and
`tests/testthat/test-cohort.R`

Add a shared wrapper:

```r
with_png_device <- function(filename, draw, width = 7, height = 5, dpi = 150) {
  if (!is.function(draw)) stop("'draw' must be a function.", call. = FALSE)
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
  device_id <- grDevices::dev.cur()
  completed <- FALSE
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(which = device_id)
    }
    if (!completed && file.exists(filename)) unlink(filename)
  }, add = TRUE)

  draw()
  grDevices::dev.off(which = device_id)
  if (!file.exists(filename) || !isTRUE(file.info(filename)$size > 0)) {
    stop(sprintf("PNG output was not written or is empty: '%s'.", filename), call. = FALSE)
  }
  completed <- TRUE
  invisible(filename)
}
```

Replace both raw `png()`/`dev.off()` pairs with `with_png_device(...,
draw = function() { ... })`. Do not use an `on.exit(dev.off())` that waits until the
entire module returns; close the exact device immediately after each successful
plot.

In the `pheatmap()` call add:

```r
cluster_rows = nrow(mat_transformed) >= 2L,
cluster_cols = ncol(mat_transformed) >= 2L,
```

Tests must intentionally throw inside `draw`, assert `dev.cur()` is restored and no
partial PNG remains, then prove a normal PNG is non-empty. A cohort composition test
with `top_n_taxa = 1L` must complete and emit the heatmap.

### Patch 5 — Correct NMDS convergence and isolate beta-resampling failures

**Files:** `analysis/05_ordination.R`, `analysis/03_beta_diversity.R`,
`tests/testthat/test-cohort.R`

Add near the top of module 05:

```r
nmds_solution_repeated <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  length(value) == 1L && !is.na(value) && is.finite(value) && value > 0
}
```

After a valid NMDS result, compute:

```r
nmds_repetitions <- as.integer(nmds_res$converged %||% 0L)
nmds_converged <- nmds_solution_repeated(nmds_repetitions)
```

Use `nmds_converged` for both `nmds_scores.tsv` and `nmds_diagnostics.tsv`, and add
`BestSolutionRepetitions = nmds_repetitions` to the diagnostics table. Do not equate
the underlying `monoMDS` stopping code with `metaMDS` best-solution repetition.

Test the helper with `0L` (false), `1L` and `2L` (true), plus `NULL`, `NA`, and a
non-numeric value (false). For a successful synthetic NMDS, require the emitted
logical field to equal `BestSolutionRepetitions > 0`.

In module 03, replace the body of each rarefaction-stability iteration with one
guarded unit:

```r
aligned <- tryCatch({
  rare_counts <- vegan::rrarefy(otu_table, sample = stability_depth)
  rare_rel <- sweep(rare_counts, 1, rowSums(rare_counts), "/")
  rare_dist <- vegan::vegdist(rare_rel, method = "bray")
  rare_fit <- stats::cmdscale(rare_dist, k = 2L, eig = TRUE, add = TRUE)
  rare_points <- as.matrix(rare_fit$points)
  if (ncol(rare_points) < 2L) {
    stop("Rarefied PCoA returned fewer than two axes.", call. = FALSE)
  }
  vegan::procrustes(reference_points, rare_points[, 1:2, drop = FALSE])$Yrot
}, error = function(e) NULL)

if (is.null(aligned)) {
  failed_iterations <- c(failed_iterations, iteration)
  next
}
```

Add a unit or controlled mock regression where Procrustes fails for one iteration;
the module must continue, report that iteration in `FailedIterations`, and retain
the successful rows. Preserve the existing fail-closed behavior when **all**
iterations fail.

### Patch 6 — Synchronize release metadata without another hard-coded verifier

**Files:** `VERSION`, `CITATION.cff`, `README.md`, `CHANGELOG.md`,
`analysis/utils/io.R`, `tests/verify_release_run.R`, and preferably a new
`tests/testthat/test-release-metadata.R`

1. Set `VERSION` and `CITATION.cff` to `0.3.1`.
2. Update all five active README `0.3.0` statements to `0.3.1`.
3. Add a new `## [0.3.1] - 2026-09-07` changelog entry above `0.3.0` summarizing
   the fixes. Do not alter historical version headings or the old v0.2.0 patch-plan
   artifact.
4. Remove release numbers from the two runtime contract errors in `io.R`:

```r
"Unsupported wf-16s classifier '%s'. This pipeline supports minimap2 only; Kraken2/Bracken requires a classifier-specific denominator model."
"Unsupported wf-16s database_set '%s'. This pipeline supports the bundled NCBI database sets only."
```

5. In `tests/verify_release_run.R`, derive the expected version from the repository
   instead of replacing one literal with another:

```r
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not locate verify_release_run.R.")
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
expected_pipeline_version <- trimws(readLines(
  file.path(repo_root, "VERSION"), n = 1L, warn = FALSE
))
```

Then use:

```r
stopifnot(identical(manifest$pipeline_version, expected_pipeline_version))
```

6. Add a metadata-consistency test that reads `VERSION`, `CITATION.cff`, the README
   current-version line, and the first changelog release heading. Require all four
   current values to agree. Match only active metadata—not historical changelog or
   archived patch-plan occurrences.

## Changelog text for v0.3.1

```markdown
## [0.3.1] - 2026-09-07

### Fixed

- Enforce whole-number upstream length and abundance thresholds and reject empty
  or whitespace-padded lineage rank fields before analysis.
- Preserve numeric-looking metadata identifiers as strings and make shared-taxa
  group handling safe for arbitrary validated labels.
- Prevent sample-seed integer overflow and reject non-portable SampleIDs.
- Close PNG devices on plotting errors and disable undefined one-object heatmap
  clustering.
- Report vegan `metaMDS` best-solution repetition correctly and isolate failed
  beta-diversity stability iterations.
- Reject duplicate/unknown modules before expensive input parsing and make release
  verification derive its expected version from `VERSION`.
```

## Required verification gates

Run from a clean checkout of the patched commit. Use fresh temporary output roots;
do not validate release behavior over an old `--overwrite` tree.

```bash
git diff --check

Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
python -m compileall -q analysis tests

Rscript tests/testthat.R
python -m unittest -v tests/test_ncbi_taxonomy.py
python tests/check_committed_whitespace.py

Rscript analysis/00_run_pipeline.R \
  --config config.yml \
  --output-dir "v0.3.1 validation output" \
  --validate-only
test ! -e "v0.3.1 validation output"

Rscript analysis/00_run_pipeline.R \
  --config config.yml \
  --output-dir "v0.3.1 integration output"
Rscript tests/verify_release_run.R "v0.3.1 integration output"

Rscript analysis/00_run_pipeline.R \
  --config tests/fixtures/synthetic_minimap2/config.yml \
  --output-dir "v0.3.1 bamstats output"
Rscript tests/verify_release_run.R "v0.3.1 bamstats output"
```

Run the FAPROTAX integration in an environment containing a compatible microeco:

```bash
Rscript analysis/00_run_pipeline.R \
  --config config.yml \
  --modules qc,alpha,beta,composition,ordination,shared,kreport,faprotax \
  --output-dir "v0.3.1 faprotax output"
Rscript tests/verify_release_run.R "v0.3.1 faprotax output"
```

Release-specific assertions:

```bash
test "$(cat VERSION)" = "0.3.1"
rg -n '0[.]3[.]0|v0[.]3[.]0' \
  VERSION CITATION.cff README.md analysis tests
```

The final `rg` command must return no active stale hit. Historical `0.3.0` entries
in `CHANGELOG.md` are expected and must remain.

The GitHub Actions matrix must finish with these four jobs green on the final SHA:

- Ubuntu / R 4.5 / Python 3.12
- Windows / R 4.5 / Python 3.12
- Ubuntu / R 4.5 / Python 3.12 (FAPROTAX)
- Windows / R 4.5 / Python 3.12 (FAPROTAX)

## Acceptance checklist

- [ ] Unknown and duplicate module requests fail before any assignment, abundance,
      metadata, bamstats, or taxonomy parsing and before output creation.
- [ ] Fractional `min_len`, `max_len`, and `abundance_threshold` fail at
      `read_upstream_params()` with the offending field named.
- [ ] Empty, whitespace-only, and boundary-whitespace taxonomy rank cells fail with
      row and rank context.
- [ ] Numeric-looking SampleIDs and groups remain exact strings; shared-taxa outputs
      preserve validated group labels.
- [ ] Boundary-whitespace, case-colliding, trailing-dot, and Windows-reserved
      SampleIDs fail closed on both supported operating systems.
- [ ] A maximum-range seed produces finite deterministic per-sample seeds; ordinary
      seed-42 rarefaction output stays byte-identical.
- [ ] An exception during base-graphics plotting leaves no leaked device or partial
      PNG.
- [ ] A one-row cohort heatmap completes without attempting row clustering.
- [ ] NMDS `Converged` equals `BestSolutionRepetitions > 0`.
- [ ] One failed Procrustes iteration is recorded without aborting successful
      stability iterations; all failed iterations still fail the module.
- [ ] `VERSION`, CFF, active README version, changelog head, release verifier, and
      manifest agree on `0.3.1`.
- [ ] All existing Ambar Ayunda invariants remain exact: 114,056 total reads; 80,556
      classified; 33,500 unclassified; 46 unresolved taxonomy nodes; 26 conflicts.
- [ ] Cache-only execution does not call the network or mutate the source taxonomy
      cache.
- [ ] Krona direct-count totals and FAPROTAX read-accounting equations remain exact.
- [ ] All four CI jobs are green on the tag target SHA.

## Reviewed but deliberately deferred beyond v0.3.1

These are real maintainability limitations, but addressing them in this patch would
either require a schema version or a higher-risk I/O redesign:

1. **Cardinality-dependent manifest JSON types.** With `auto_unbox = TRUE`, fields
   such as `samples`, module `outputs`, and CLI `modules` can be a scalar for one
   item and an array for several. Stabilize these in a future manifest schema version
   rather than silently changing schema v1 in a patch release.
2. **Whole-file assignment buffering.** `read_assignments_file()` loads every line,
   constructs a second collapsed text buffer, and then parses it. This is safe for
   the tracked 14 MiB fixture but scales poorly to very large PromethION assignment
   files. Replace it with a streaming/schema-pass plus direct typed import in a
   dedicated performance patch with memory benchmarks.
3. **Incomplete requested-module records after fail-fast execution.** When a module
   fails without `--keep-going`, later requested modules are absent rather than
   explicitly recorded as `not_run`. Adding that state should be coordinated with a
   manifest schema revision.
4. **Environment locking.** The lack of a committed `renv.lock` remains accurately
   documented. Do not smuggle a machine-specific lock into this correctness patch.

## Handoff verdict

This plan is intentionally additive to the v0.3.0 architecture: it preserves the
validated classified-only denominators, exact read accounting, offline taxonomy
policy, run-local resolver cache, Kraken/Krona arithmetic, FAPROTAX interpretation,
and existing single/cohort gates. Once the required patches and regressions pass all
four jobs on one final SHA, v0.3.1 is suitable to tag.
