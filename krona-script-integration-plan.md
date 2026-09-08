# Krona Script Integration Plan

## Status

**Planning only.** This document captures the proposed implementation; no
pipeline source files are changed by this handoff. The current repository
remains at version `0.2.0` with its default module behavior unchanged.

## Objective

Add an opt-in Krona-compatible taxonomy visualization output to the existing
`kreport` pathway. The integration should:

1. reuse the already validated abundance table, taxonomy tree, unclassified
   read accounting, and run-local taxonomy resolution;
2. emit a standard tab-delimited input accepted by KronaTools `ktImportText`;
3. optionally invoke `ktImportText` to create a self-contained interactive
   `.krona.html` file;
4. preserve the existing six-column `.kreport` output and all default v0.2.0
   release behavior;
5. fail closed on invalid count/path invariants and make renderer availability
   explicit in provenance and manifest warnings.

## Proposed user experience

### Command-line use

```bash
Rscript analysis/00_run_pipeline.R \
  --config config.yml \
  --modules kreport \
  --krona \
  --overwrite
```

The `--krona` switch enables the extension and implicitly requires the
`kreport` module. It should not add a separate module name to the module
registry: Krona and `.kreport` must be produced from the same validated tree.

### YAML use

```yaml
krona:
  enabled: true
  render_html: true
  executable: "ktImportText"
```

`render_html: false` requests the portable Krona-compatible TSV only. This is
the recommended mode on CI or installations where KronaTools is not present.

### Outputs

For each selected sample, under `07_Kreport/krona/`:

```text
<sample>.krona.tsv
<sample>.krona.html       # when ktImportText is available and rendering is enabled
```

The module should additionally write:

```text
07_Kreport/krona/krona_provenance.json
```

The existing `.kreport`, taxonomy diagnostics, and taxonomy provenance files
remain in their current locations.

## Format and counting contract

### Krona input representation

`ktImportText` consumes one tab-delimited line per contribution:

```text
<magnitude>\t<rank-1 label>\t<rank-2 label>\t...\t<leaf label>
```

The first field is a non-negative integer magnitude. Krona creates the root
itself; therefore the exporter must **not** add a synthetic `root` row or any
ancestor rows whose counts already include descendants.

Example for a sample with 20 unclassified reads and two classified terminal
taxa:

```text
20	Unclassified
30	Bacteria	Bacillati	Bacillota	Bacilli	Bacillales	Bacillaceae	Bacillus	Bacillus_sp1
50	Bacteria	Bacillati	Bacillota	Bacilli	Bacillales	Bacillaceae	Bacillus	Bacillus_sp2
```

The exporter should use the terminal/direct contribution in
`nodes_df$reads_taxon`, not `reads_clade`. This avoids double-counting a
lineage when both its ancestors and leaves are present in the tree. If the
input taxonomy contains a valid internal taxon with direct abundance, its
direct count is emitted using the complete path to that internal node; it is
not discarded merely because it has children.

### Read-accounting invariant

For every sample:

```text
sum(Krona magnitudes) == TotalReads
sum(classified terminal/direct contributions) == ClassifiedReads
unclassified contribution == UnclassifiedReads
```

The existing `validate_kreport_tree()` must pass before Krona output is
written. The Krona builder should independently validate finite,
non-negative, integer-valued totals and the exact final magnitude sum so that
it cannot silently produce a misleading chart if called outside the full
runner.

### Label handling

- Preserve the canonical eight-rank path order already used by the pipeline.
- Convert missing/empty rank labels to a visible placeholder such as
  `Unknown`; do not emit an empty Krona wedge.
- Reject labels containing tab, carriage-return, or newline characters because
  they would corrupt the physical tab-delimited record.
- Do not manually generate XML or HTML. KronaTools is responsible for escaping
  labels in the HTML output.
- Keep `Unclassified` as an explicit top-level contribution whenever its count
  is positive. Do not encode the canonical unclassified placeholder ranks as a
  biological lineage.

## Implementation changes

### 1. `analysis/utils/config.R`

Add this block to `get_default_config()`:

```r
krona = list(
  enabled = FALSE,
  render_html = TRUE,
  executable = "ktImportText"
)
```

Extend `validate_config()` with fail-closed checks:

```r
if (!is.logical(cfg$krona$enabled) || length(cfg$krona$enabled) != 1L ||
    is.na(cfg$krona$enabled)) {
  stop("'krona.enabled' must be true or false.", call. = FALSE)
}
if (!is.logical(cfg$krona$render_html) || length(cfg$krona$render_html) != 1L ||
    is.na(cfg$krona$render_html)) {
  stop("'krona.render_html' must be true or false.", call. = FALSE)
}
assert_nonempty_string(cfg$krona$executable, "krona.executable")
```

In `load_config()`:

- recognize `--krona` through both `cli_opts$krona` and
  `cli_opts[["krona"]]`;
- set `cfg$krona$enabled <- TRUE` when the switch is present;
- record the effective setting in `cfg$cli`, for example:

```r
krona = isTRUE(cli_opts$krona) || isTRUE(cli_opts[["krona"]])
```

Do not resolve `krona.executable` as a filesystem path in configuration: it
may be either a command on `PATH` or an absolute/relative executable path.
Resolve it only through the dependency helper described below.

### 2. `analysis/utils/cli.R`

Add an optparse boolean option:

```r
optparse::make_option(
  c("--krona"),
  action = "store_true",
  default = FALSE,
  dest = "krona",
  help = "Enable Krona-compatible TSV output and optional KronaTools HTML rendering"
)
```

The help text should mention that HTML rendering depends on KronaTools, while
the TSV format itself is generated by R and has no additional R package
dependency.

### 3. `analysis/utils/dependencies.R`

Add executable discovery that works on Unix and Windows and accepts either a
command name or an explicit path:

```r
find_krona_executable <- function(executable = "ktImportText") {
  if (is.null(executable) || length(executable) != 1L ||
      is.na(executable) || !nzchar(trimws(executable))) {
    return(NA_character_)
  }

  if (file.exists(executable) || grepl("[/\\\\]", executable)) {
    if (!file.exists(executable)) return(NA_character_)
    return(normalizePath(executable, winslash = "/", mustWork = TRUE))
  }

  found <- Sys.which(executable)
  if (!nzchar(found)) NA_character_ else
    normalizePath(found, winslash = "/", mustWork = TRUE)
}
```

Add a renderer check with explicit semantics:

- if Krona is disabled, do nothing;
- if only TSV output is requested, do not require KronaTools;
- if HTML rendering is enabled and the executable is unavailable, do not
  pretend an HTML chart was produced. The module should emit the TSV, record
  an explicit warning and `html_status: "renderer_missing"` in provenance,
  and list the warning in the run manifest;
- if the project later chooses strict chart delivery, a future
  `require_html: true` switch can convert this warning into a preflight error.
  Do not introduce that stricter behavior implicitly into the v0.2.0-compatible
  default.

An optional version helper may invoke `ktImportText` with no input and parse a
`KronaTools x.y.z` token from stdout/stderr. Version probing is informational:
a missing or nonstandard version string must not invalidate an otherwise
successful TSV export.

### 4. `analysis/utils/kreport.R`

Add the following helpers.

#### `build_krona_lines()`

Responsibilities:

1. Validate `total_reads`, `uncl_reads`, and all direct counts as finite,
   non-negative integer quantities.
2. Calculate `classified_reads <- total_reads - uncl_reads` and reject a
   negative classified total.
3. Emit one `Unclassified` line only when `uncl_reads > 0`.
4. Select `nodes_df$reads_taxon > 0`, ordered in the already deterministic
   DFS order supplied by `build_kreport_tree()`.
5. Split each node path on `;`, normalize empty labels to `Unknown`, validate
   physical delimiter characters, and produce one magnitude-plus-lineage
   record.
6. Confirm that the sum of emitted magnitudes equals `total_reads` exactly.
7. Return data lines only; comments/metadata belong in provenance, avoiding
   ambiguity for downstream parsers and tests.

Use integer-style formatting (`formatC(..., format = "f", digits = 0)`) so
large whole-read counts are not written in scientific notation or with a
decimal suffix.

#### `write_krona_input()`

Create the parent directory, call `build_krona_lines()`, write the lines with
`writeLines()`, and verify that the file exists and is non-empty. Return the
path or the generated lines invisibly for testability.

#### `render_krona_html()`

Invoke the executable using `processx::run()` with an argument vector, never a
shell-constructed command string:

```r
args <- c("-o", output_path, "-n", sample_id, input_path)
```

On non-zero exit status, missing output, or zero-byte output, stop with a
message that includes the sample and executable. Verify only stable output
properties (for example, that the generated file is non-empty); do not make a
brittle assertion about a particular Krona HTML template or JavaScript
version.

### 5. `analysis/07_kreport_pavian.R`

Extend `run_kreport(context)` after the existing tree validation and `.kreport`
write:

```r
krona_cfg <- cfg$krona %||% list(
  enabled = FALSE,
  render_html = FALSE,
  executable = "ktImportText"
)
krona_enabled <- isTRUE(krona_cfg$enabled)
render_html <- krona_enabled && isTRUE(krona_cfg$render_html)
```

If enabled:

- create `07_Kreport/krona/`;
- discover `krona_cfg$executable` once;
- write one `.krona.tsv` for each sample from the same `nodes_sorted`,
  `total_reads`, and `uncl_reads` used for `.kreport`;
- if rendering is requested and the executable exists, write one
  `.krona.html` per sample;
- if rendering is requested and the executable is absent, issue one clear
  warning, continue with TSV output, and record the skipped renderer status;
- append every created artifact to the module's `outputs` vector.

Write one provenance object after all samples. It should include at least:

```r
list(
  format = "KronaTools ktImportText tab-delimited lineage format",
  renderer = if (html_rendered) "ktImportText" else NULL,
  krona_tools_version = if (html_rendered) krona_version else NULL,
  html_status = html_status,
  standalone_html = if (html_rendered) TRUE else NULL,
  count_model = "direct abundance-table taxon counts plus canonical unclassified count",
  denominator = "TotalReads",
  classified_definition = "sum of direct positive-count classified taxonomy rows",
  samples = sample_records
)
```

Each sample record should contain the sample ID, total reads, classified reads,
unclassified reads, emitted magnitude sum, TSV path, and (when present) HTML
path. Use `jsonlite::write_json(..., pretty = TRUE, auto_unbox = TRUE,
null = "null")`.

The helper should remain safe for current unit-test contexts that construct a
partial `cfg$output$dirs` or omit `cfg$krona`; the `%||%` fallback must leave
Krona disabled in those cases.

### 6. `analysis/00_run_pipeline.R`

Add `--krona` to the effective CLI/config flow, but do not register a new
module. Before output mutation:

- if Krona is enabled and `kreport` is not in `requested_modules`, stop with a
  clear error explaining the dependency;
- do not require an external executable merely to produce the TSV;
- if a strict HTML requirement is added later, perform its preflight here so
  validation-only runs fail before any output mutation.

The existing manifest will capture the effective `cli$krona` setting and the
Krona artifacts through `module_results$kreport$outputs`. Module warnings must
include a missing-renderer warning rather than disappearing into console-only
output.

## Tests to add

### `tests/testthat/test-config.R`

Add assertions that:

- defaults contain `cfg$krona$enabled == FALSE`;
- `render_html` is logical and defaults to the documented value;
- `krona.executable` defaults to `ktImportText`;
- `load_config(..., cli_opts = list(krona = TRUE))` sets both
  `cfg$krona$enabled` and `cfg$cli$krona` to `TRUE`;
- invalid `enabled`, `render_html`, or empty executable values fail closed;
- unknown nested Krona keys are rejected by `merge_config()`.

### `tests/testthat/test-kreport.R`

Add pure-R tests that do not require KronaTools:

1. **Terminal-count accounting:** two species under one lineage must produce
   two leaf contributions plus the unclassified contribution; no ancestor
   clade count may appear as a separate record.
2. **Exact sum:** emitted magnitudes equal `TotalReads`, with the expected
   classified/unclassified partition.
3. **Zero unclassified:** no `Unclassified` line is emitted when its count is
   zero, while the classified sum remains exact.
4. **Path safety:** tabs/newlines in a label are rejected; empty labels are
   represented as `Unknown` if normalization is chosen.
5. **Invalid arithmetic:** negative classified totals, non-integer counts, and
   mismatched total sums fail closed.
6. **File output:** `write_krona_input()` creates a valid physical TSV and
   preserves sample paths containing spaces.
7. **Integrated opt-in:** a small `run_kreport()` context with
   `render_html = FALSE` produces both the usual `.kreport` and the expected
   Krona TSV/provenance outputs.
8. **Backward compatibility:** existing contexts with no `cfg$krona` continue
   to produce the exact pre-Krona output set.

If a `ktImportText` executable is available in the test environment, add an
optional integration test guarded with `skip_if_not()` that checks the HTML
file exists and is non-empty. Do not make the core CI suite depend on a global
KronaTools installation.

### Release-verification logic

Do not alter the default v0.2.0 release assertions merely because the feature
exists. Add a conditional verifier branch only when a manifest was generated
with `manifest$cli$krona == TRUE`; it should then require:

- one Krona TSV per selected sample;
- `krona_provenance.json`;
- exact per-sample magnitude sums and read partitions;
- HTML files only when provenance says HTML rendering succeeded.

This keeps the historical release fixture independent of an optional external
tool while making explicitly Krona-enabled runs auditable.

## Documentation changes

### `README.md`

- Update the architecture line for `07_kreport_pavian.R` to mention the
  Krona-compatible export.
- Add a “Krona export (opt-in)” subsection under the Kraken/Pavian section.
- Explain that R always creates the native `.krona.tsv` when enabled, while
  `.krona.html` requires KronaTools `ktImportText`.
- Document both CLI and YAML examples, output paths, the direct-count model,
  explicit unclassified contribution, and the fact that Krona's root is
  generated by the renderer.
- Add `--krona` to the CLI options table.
- State that a Krona chart is a visualization of the classifier-conditioned
  abundance table, not independent taxonomic validation.
- Remove or qualify any release-scope statement that says automatic Krona HTML
  export is absent once this feature is merged. If HTML remains unavailable
  without an external installation, say so precisely.

### `config.example.yml`

Add a commented opt-in block:

```yaml
krona:
  enabled: false
  render_html: true
  executable: "ktImportText"
```

The tracked `config.yml` may remain disabled so the regression fixture does not
gain a new external-tool dependency.

### `CHANGELOG.md`

If this is merged after v0.2.0 is tagged, add an `[Unreleased]` entry rather
than silently folding a new output contract into the historical release. The
entry should distinguish the always-available native TSV from the optional
KronaTools HTML renderer.

## Edge cases and decisions

| Case | Required behavior |
|---|---|
| No KronaTools installed | Emit TSV; warn explicitly; provenance says renderer missing; no false HTML path is recorded. |
| `render_html: false` | Emit TSV only; no executable lookup is needed. |
| `--krona` without `kreport` | Fail before module execution with a configuration/CLI error. |
| Unclassified reads > 0 | Emit a single top-level `Unclassified` contribution and include it in the total. |
| Unclassified reads = 0 | Omit that line; classified contributions must still equal total reads. |
| Internal taxon has direct reads and children | Emit its direct count at its full lineage path; descendants remain separate. |
| Ancestor clade has no direct reads | Never emit the ancestor's `reads_clade` as a magnitude. |
| Empty rank label | Normalize to `Unknown` or reject consistently; never emit an empty field. |
| Tab/newline in a label | Reject because it corrupts the input record. |
| Sample ID/path contains spaces | Pass arguments through `processx` vectors and use `sanitize_filename()` only for filenames. |
| Zero total reads | Already rejected by `build_context()`; helper should also fail closed. |
| TaxID unresolved | Krona still renders the lineage labels; provenance must not imply TaxID completeness. |
| Taxonomy cache conflict | Preserve the existing kreport diagnostics; do not silently resolve or collapse conflicting paths. |
| Large counts | Write integer text, not scientific notation or floating-point decimals. |
| HTML renderer fails for one sample | Fail the module rather than reporting a successful HTML artifact for that sample; retain the already-created TSVs and record the error through the normal module failure path. |

## Acceptance criteria

The integration is ready to merge when all of the following hold:

- default `config.yml` and the existing release integration output set are
  unchanged;
- `Rscript analysis/00_run_pipeline.R --config config.yml --validate-only`
  still performs zero filesystem mutations;
- core CI passes on Ubuntu and Windows without KronaTools installed;
- pure-R tests demonstrate exact Krona magnitude accounting and rejection of
  malformed records;
- an enabled run produces one valid `.krona.tsv` per selected sample and a
  provenance JSON whose counts reconcile to the existing accounting outputs;
- an installation with official KronaTools produces non-empty interactive
  HTML through `ktImportText` using argument-vector invocation;
- missing renderer availability is visible in the module warning, manifest,
  and Krona provenance rather than inferred from absent files;
- the conditional release verifier detects missing or extra Krona artifacts
  for explicitly Krona-enabled runs;
- README and example configuration describe the external renderer boundary and
  the unclassified/direct-count semantics accurately.

## Suggested implementation order

1. Add configuration and CLI schema/validation.
2. Add executable discovery and renderer helpers.
3. Implement and unit-test the pure-R Krona line builder.
4. Integrate TSV/provenance generation into `run_kreport()`.
5. Add optional HTML invocation and failure handling.
6. Add conditional manifest/release-verifier checks.
7. Update README, example config, and changelog.
8. Run parsing, Python compilation, testthat, validate-only, full fixture
   integration, and—on a KronaTools-enabled environment—an HTML smoke test.

## References

- [KronaTools official repository](https://github.com/marbl/Krona)
- [`ktImportText` input/output documentation](https://nf-co.re/modules/krona_ktimporttext)
- Ondov et al. (2011), *Krona: hierarchical data visualization in a browser*,
  [BMC Bioinformatics 12:385](https://doi.org/10.1186/1471-2105-12-385)
