# ONT wf-16S Amplicon Post-Processing: Audited Phase 1 Implementation Plan

**Document:** `ont-wf16s-amplicon-postprocess-xpand1.md`
**Target repository:** `WhyAdr/ont-wf16S-amplicon-postprocess`
**Audited:** September 2, 2026
**Status:** Ready for implementation, subject to the gates below
**Reference workspace:** `ATW_Sesame_Greenhouse-Examples/` (read-only design reference; not part of this repository)

---

## 1. Audit verdict

The modernization is worthwhile, but the previous draft was not safe to execute as written. It mixed incomplete pseudo-code with implementation instructions, overstated a few analytical capabilities, and omitted several input and failure contracts. This revision is the execution specification.

The implementation must correct these issues:

1. `load_config()` defaulted to `../config.yml`, while the documented command runs from the repository root. That resolves to the parent directory, not this repository.
2. The proposed `--output-dir` override changed only `output.base_dir`; every module still used the original configured subdirectory.
3. The composition pseudo-code referenced a `count` column that does not exist in the aggregated wf-16S table until a sample column is selected or the table is pivoted.
4. The runner caught module errors, continued, and printed a successful finish message. The default must be fail-fast and return a non-zero process status.
5. `mode: auto` treated the mere presence of metadata as cohort mode, even for a one-sample table.
6. One `assignments_file` cannot represent a cohort because wf-16S publishes read assignments per sample.
7. Hellinger PCA is not meaningful with one sample. Single-sample PCA must be skipped, not described as a species decomposition.
8. Repeated `vegan::rrarefy()` draws are rarefaction resamples without replacement, not bootstrap replicates and never biological replicates.
9. PERMANOVA needs replication/design gates and a dispersion diagnostic. It must not run merely because two sample columns exist.
10. The real second lineage field is `kingdom`, confirmed by the tracked alignment table, not a generic `clade`. Kraken's standard rank code for this field is `K`, not `D1`.
11. The NCBI resolver lacked `tool`/`email`, timeouts, retries, atomic cache writes, subprocess status checks, and a deliberate offline mode.
12. The draft promised read-quality diagnostics, but the tracked per-read assignments contain length, not per-read quality. Phase 1 must not claim observed quality analysis.
13. The draft promised a standalone Pavian HTML export, but the proposed code only generated a `.kreport`. Phase 1 will generate and validate reports; Pavian upload/export remains a documented manual step.
14. The package installer could exit successfully after failed installations.
15. Adding `figs/`, `tables/`, and `*.kreport` to `.gitignore` would not untrack existing artifacts and would hide potentially useful future fixtures.
16. A direct `git push origin main` is not a validation step and is not authorized by this plan.

The local ATW project is a useful source of architectural ideas, not a drop-in implementation. In particular, do not copy its current-working-directory assumptions or fail-open runner behavior. Reimplement the needed behavior in this repository and preserve attribution if any non-trivial code is copied after checking its license.

Current wf-16S already produces interactive Sankey and sunburst views in its workflow report. The custom `.kreport` module is therefore a Pavian interoperability/export feature, not a claim that this repository uniquely adds taxonomic visualization.

---

## 2. Ground truth that the implementation must preserve

### 2.1 Repository snapshot at audit time

- Branch: `main`
- Commit: `7ec4d4d`
- Remote state: `main == origin/main`
- Tracked modernization plan: this file
- Untracked content: `ATW_Sesame_Greenhouse-Examples/` only
- Existing root scripts:
  - `analyze_16s_improved.R`: 603 lines
  - `analyze_16s.R`: 379 lines
  - `convert_to_kreport.R`: 164 lines
  - `fetch_ncbi_taxonomy.py`: 180 lines

The ATW directory is a nested, untracked reference checkout. Do not stage it, modify it, or make the target pipeline depend on it.

### 2.2 Real wf-16S fixture contract

The tracked Ambar Ayunda files establish these regression invariants:

| Invariant | Expected value |
|---|---:|
| Abundance rows | 1,837 |
| Sample columns after excluding `tax` and `total` | 1 |
| Sample ID | `AmbarAyunda_minimap2_16S` |
| Total reads | 114,056 |
| Classified reads | 80,556 |
| Unclassified reads | 33,500 |
| Assignment rows | 114,056 |
| Raw assignment status `C` | 89,809 |
| Raw assignment status `U` | 24,247 |
| Rows with status `C` and TaxID `0` | 9,253 |
| Effective classified rows (`taxid > 0`) | 80,556 |

The equality `80,556 + 33,500 = 114,056` and the agreement between abundance and assignment classifications are mandatory validation checks. For this minimap2 output, the status letter alone overstates classification; TaxID `0` is the effective unclassified signal.

The abundance table contract is:

- Tab-delimited.
- `tax` contains exactly eight semicolon-delimited fields for this Phase 1 schema:
  `superkingdom`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`.
- One or more sample columns contain non-negative, integer-valued counts, even if serialized as decimals such as `33500.0`.
- `total` is an aggregate column, not a sample. When present, it must equal the row sum across all detected sample columns before any optional sample subset is selected.
- The synthetic unclassified path is
  `Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown`.

The supplied minimap2 assignment contract is five tab-delimited fields with no header:

1. raw status (`C` or `U`),
2. read ID,
3. TaxID (`0` means effectively unclassified),
4. a length field such as `0|1481`,
5. pipe-delimited lineage.

Parse the read length defensively from either a single integer or the last numeric component of a pipe-delimited value. Reject malformed non-empty values with line-numbered diagnostics; do not silently coerce them to `NA`.

In the real fixture, assignment lineages contain seven named ranks and omit the abundance table's `kingdom` field. Exact cross-file lineage matching therefore uses the normalized abundance key `superkingdom;phylum;class;order;family;genus;species`; it must not use species name alone.

The tracked taxonomy cache contains 3,220 entries: 3,130 non-zero TaxIDs and 90 zero/unresolved entries, including the synthetic unclassified lineage. A successful offline run must report this incompleteness rather than claim that every node has an NCBI TaxID.

### 2.3 Numerical regression targets

For the current real input, classified-only alpha diversity should reproduce the established values within a documented floating-point tolerance:

| Metric | Expected value |
|---|---:|
| Observed species richness | 1,836 |
| Chao1 | 2,983.851 |
| Shannon | 4.603 |
| Effective species number | 99.814 |
| Simpson (`1 - sum(p^2)`) | 0.957 |
| Inverse Simpson | 23.089 |
| Pielou evenness | 0.613 |
| Fisher alpha | 334.543 |
| Berger-Parker dominance | 0.133 |

Detected classified taxa by rank should remain 31 phyla, 69 classes, 128 orders, 281 families, 867 genera, and 1,836 species.

These are regression targets, not claims that every label is biologically validated to species level. They reflect the supplied classifier/database output.

---

## 3. Phase 1 scope and non-goals

### In scope

- A strict, config-driven runner callable from the repository root.
- One validated data model shared by all modules.
- Single-sample parity for the supplied minimap2 dataset.
- Cohort-capable composition, alpha diversity, beta diversity, ordination, and shared-taxa modules when a genuine multi-sample abundance table and metadata are supplied.
- Per-sample assignment mapping for optional QC and TaxID support.
- Deterministic outputs, a resolved-config snapshot, input hashes, package/session information, and a machine-readable run manifest.
- Offline `.kreport` generation using the existing cache, plus an explicit opt-in taxonomy refresh path.
- Unit tests with small synthetic fixtures and a network-free real-input regression test.
- Updated README and configuration examples.

### Out of scope

- Re-running wf-16S, basecalling, demultiplexing, or changing upstream filtering.
- Treating rarefaction iterations as biological replicates.
- Differential abundance testing, supervised learning, functional prediction, phylogenetic/UniFrac analysis, or causal inference.
- Automatic standalone Pavian HTML export.
- Per-read quality analysis unless a future tracked input contract supplies per-read quality values.
- Supporting arbitrary historical wf-16S schemas without an explicit adapter and fixture.
- Deleting or untracking legacy scripts, figures, tables, or the existing Pavian HTML.
- Committing or pushing changes unless the user separately authorizes those Git actions.

The current repository has no real cohort fixture. Cohort code can receive contract and synthetic tests in Phase 1, but it must not be reported as biologically validated until it is run on a genuine cohort with appropriate replication and metadata.

---

## 4. Target layout

```text
.
|-- config.yml
|-- config.example.yml
|-- metadata.example.tsv
|-- README.md
|-- ont-wf16s-amplicon-postprocess-xpand1.md
|-- analysis/
|   |-- 00_run_pipeline.R
|   |-- 01_qc_diagnostics.R
|   |-- 02_alpha_diversity.R
|   |-- 03_beta_diversity.R
|   |-- 04_taxa_composition.R
|   |-- 05_ordination.R
|   |-- 06_shared_taxa.R
|   |-- 07_kreport_pavian.R
|   |-- install_packages.R
|   `-- utils/
|       |-- cli.R
|       |-- config.R
|       |-- dependencies.R
|       |-- io.R
|       |-- metrics.R
|       |-- plotting.R
|       |-- kreport.R
|       `-- ncbi_taxonomy.py
|-- tests/
|   |-- testthat.R
|   |-- fixtures/
|   `-- testthat/
`-- output/                         # generated; ignored as one root-anchored directory
    |-- 01_QC/
    |-- 02_Alpha_Diversity/
    |-- 03_Beta_Diversity/
    |-- 04_Taxa_Composition/
    |-- 05_Ordination/
    |-- 06_Shared_Taxa/
    |-- 07_Kreport/
    |-- resolved_config.yml
    |-- run_manifest.json
    `-- session_info.txt
```

Each numbered R file must define a function such as `run_qc(context)` and must not execute analysis merely because it was sourced. Function scope provides isolation; do not use sourced scripts with hidden global-state contracts. The runner builds one validated `context` object and calls selected module functions explicitly.

Legacy root scripts remain untouched during Phase 1 parity work. After the new runner passes all acceptance gates, the README may mark them as legacy. Replacing them with wrappers or deleting them is a separate decision.

---

## 5. Configuration contract

Use a single configured output root and derive module subdirectories in code. Do not configure both `base_dir` and seven duplicated output paths.

```yaml
schema_version: 1
project_name: "AmbarAyunda_16S_Amplicon"
mode: "auto"                  # auto, single, or cohort
seed: 42

input:
  abundance_table: "output_AAy/abundance_table_species.tsv"
  metadata: null
  params_json: "output_AAy/params.json"
  tax_column: "tax"
  aggregate_columns: ["total"]
  include_samples: null       # null means all non-aggregate sample columns
  assignments:
    AmbarAyunda_minimap2_16S: "output_AAy/reads_assignments/AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv"

output:
  base_dir: "output"

qc:
  display_min_length: 1200
  display_max_length: 1800
  target_min_length: null     # prefer params.json min_len when available
  target_max_length: null     # prefer params.json max_len when available

alpha:
  rarefaction_points: 25
  resample_depth: 50000       # capped below each sample's classified depth
  resample_fraction_cap: 0.90
  resample_iterations: 100

composition:
  top_n_taxa: 15
  heatmap_rank: "genus"
  heatmap_transform: "log10_relative"

beta:
  distances: ["bray", "jaccard"]
  permutations: 999
  strata_column: null
  resampling:
    enabled: false
    iterations: 100
    depth_fraction_of_minimum: 0.75

shared_taxa:
  rank: "species"
  minimum_count: 1
  group_prevalence: 0.5

taxonomy:
  cache: "output_AAy/taxonomy_cache.json"
  network_mode: "cache_only"  # online refresh requires --refresh-taxonomy
  unresolved_policy: "warn"  # warn or error
  email_env: "NCBI_EMAIL"
  api_key_env: "NCBI_API_KEY"
```

All relative paths are resolved against the directory containing the selected config file, never against the caller's current working directory. Absolute paths remain absolute.

CLI precedence is:

1. explicit CLI option,
2. YAML value,
3. documented built-in default.

Required runner options:

```text
--config PATH
--output-dir PATH
--modules qc,alpha,beta,composition,ordination,shared,kreport
--validate-only
--keep-going
--refresh-taxonomy
--overwrite
```

Rules:

- Default execution is fail-fast. `--keep-going` is explicit and the final exit status is still non-zero if any requested module fails.
- `--validate-only` performs parsing, dependency checks, schema checks, sample/metadata alignment, and input reconciliation without creating or changing the output tree.
- `--output-dir` replaces the base root before all module paths are derived.
- Without `--overwrite`, refuse to replace an existing known output file. Never recursively delete the output directory.
- `--refresh-taxonomy` is the only CLI action that permits NCBI network access and overrides `cache_only` for that run.

### Mode resolution

- `auto`: one selected sample means single mode; two or more selected samples means cohort mode.
- Metadata alone never converts a one-sample table into a cohort.
- `single`: require exactly one selected sample.
- `cohort`: require at least two selected samples and metadata.
- Metadata must contain unique, non-empty `SampleID` and `Group` values and must match the selected sample set exactly. Report missing, extra, and duplicate IDs separately.

`metadata.example.tsv` should demonstrate the contract but must be clearly synthetic. It must not be used as a real analysis input.

---

## 6. Shared data and validation layer

Implement this layer before any plotting module.

### 6.1 Context object

The validated context should contain, at minimum:

- resolved configuration and derived output paths,
- taxonomy table with the original full lineage and eight named ranks,
- taxa-by-sample numeric count matrix,
- selected sample IDs in canonical order,
- optional metadata aligned to that exact order,
- optional per-sample assignment paths,
- per-sample total, classified, and unclassified counts,
- input SHA-256 values,
- warnings and validation results.

No module should independently redetect sample columns or reinterpret the taxonomy schema.

### 6.2 Abundance validation

Fail with actionable diagnostics when any of the following occurs:

- missing or duplicate column names,
- missing configured tax column,
- zero selected samples,
- non-numeric, negative, non-finite, or non-integer-valued counts,
- duplicate full lineage rows,
- lineage field count other than eight,
- multiple unclassified rows per sample contract or no recognizable unclassified row,
- an aggregate `total` that differs from the row sum across all detected sample columns before subsetting,
- a sample with zero total reads or zero classified reads where a requested module requires classified counts.

Reject sample IDs that are empty, contain control characters or path separators, equal `.` or `..`, or would collide after output-filename sanitization.

Use full lineage prefixes as taxon keys at every rank. Do not collapse every `Unknown` family/genus from unrelated parents into one taxon. Tables should include both a stable `TaxonPath` and a display `Taxon` label; contextualize ambiguous labels with their parent.

### 6.3 Assignment validation

Assignments are optional because wf-16S only publishes them when requested. If the mapping is `null`, abundance-only modules continue and QC reports a clear skip. If a sample mapping is supplied, a missing or malformed file is an error.

For every supplied assignment file:

- validate exactly five fields per non-empty row,
- validate status and integer TaxID,
- parse length defensively,
- require unique read IDs,
- define effective classification as `taxid > 0`, while retaining raw status for the discrepancy diagnostic,
- reconcile assignment row count, classified count, and unclassified count with the abundance table.

Default reconciliation tolerance is zero. If future wf-16S filtering semantics require a tolerance or documented exception, add a versioned adapter and fixture rather than weakening the default silently.

---

## 7. Module specifications

### 7.1 `01_qc_diagnostics.R`

For each sample with assignments, produce:

- classification donut using effective classification,
- read-length histogram by effective class,
- raw-status versus effective-classification table/plot,
- violin/box plot of read length by effective class,
- `read_length_by_status.tsv`,
- `classification_reconciliation.tsv`.

Length ranges are plot annotations, not downstream filters. Prefer `min_len` and `max_len` from the supplied `params.json`; use config fallbacks only when unavailable. Do not label these values as filters applied by this post-processing pipeline.

In cohort mode, write per-sample metrics and a cohort summary; avoid rendering dozens of unreadable panels into one figure. No per-read quality plot belongs in Phase 1.

### 7.2 `02_alpha_diversity.R`

Compute classified-only metrics independently for every biological sample:

- observed richness,
- Chao1,
- Shannon,
- effective species number,
- Simpson diversity,
- inverse Simpson,
- Pielou evenness,
- Fisher alpha,
- Berger-Parker dominance.

Handle edge cases explicitly: zero classified reads, one detected taxon, failed estimator convergence, and depths below the requested resampling depth.

Produce exact analytical rarefaction curves with `vegan::rarefy()`. Repeated `rrarefy()` output must be named `rarefaction_resamples`, not bootstrap output. Use
`min(configured_depth, floor(classified_depth * fraction_cap))` in single mode. In cohort mode, use one common depth based on the smallest classified library so sample comparisons are standardized. Record the actual depth used.

For cohort mode:

- keep raw per-sample metrics as the primary output,
- optionally summarize rarefied sensitivity across iterations,
- never pass the iterations to group tests as independent observations,
- run group tests only when at least two groups have at least three biological samples each,
- use Kruskal-Wallis for more than two groups and Wilcoxon for an explicitly requested two-group contrast,
- apply Benjamini-Hochberg correction across tested metrics,
- write a skip reason instead of manufacturing a p-value when gates are not met.

Record the upstream `abundance_threshold` from `params.json` when available because richness estimators are sensitive to removal of rare taxa.

### 7.3 `03_beta_diversity.R`

This module is cohort-only.

Primary analysis:

- samples as rows, taxa as columns,
- classified relative abundance for Bray-Curtis,
- presence/absence for Jaccard using the configured minimum-count rule,
- deterministic full-data distance matrices and PCoA coordinates,
- two-dimensional plots only when at least three non-identical samples permit two axes,
- `adonis2` PERMANOVA only with at least two groups, at least two samples per group, and residual degrees of freedom,
- warning in the report when any group has fewer than three samples,
- `betadisper` plus permutation test alongside PERMANOVA,
- optional permutation strata only when a configured metadata column is present and valid.

The deterministic full-data PCoA is the primary result. If optional resampling is enabled, rarefy at the configured fraction of the minimum classified library, Procrustes-align iterations, and export it as a separate stability analysis. Do not label consensus coordinates with eigenvalue percentages from a different ordination without an explicit explanation.

Seed every stochastic operation from the run seed and record the seed and realized resampling depth.

### 7.4 `04_taxa_composition.R`

Use one long-format table produced from the shared count matrix. For each configured rank:

- aggregate by lineage prefix, not label alone,
- calculate per-sample relative abundance with an explicit denominator,
- export complete count and relative-abundance tables,
- group plotted low-abundance taxa into `Other` without altering exported complete tables,
- keep `Unclassified` separate from classified composition.

Required figures/tables:

- all-read classification fraction,
- classified-only phylum, family, genus, and species composition,
- top-taxa bars for a single sample,
- stacked cohort composition when multiple samples exist,
- cohort heatmap only when at least two samples exist.

For heatmaps, select taxa by mean classified relative abundance, apply the configured transform with a documented pseudocount, and annotate samples from aligned metadata. Do not imply significance from a descriptive heatmap.

### 7.5 `05_ordination.R`

This module is cohort-only. Skip with a structured reason in single mode.

- Hellinger-transform classified relative abundances and run PCA across samples.
- Use at most `min(2, n_samples - 1, matrix_rank)` displayed axes.
- Export scores, loadings, and variance explained.
- Label only a configured number of highest-contributing taxa.
- Run Bray-Curtis NMDS only when at least three non-identical samples are available.
- Record convergence, number of tries, and stress; warn prominently when stress is high.
- Do not add automatic `envfit` screening in Phase 1.

### 7.6 `06_shared_taxa.R`

This module requires at least two samples and is most useful with at least two groups.

Export:

- sample-level presence/absence matrix,
- group-level prevalence table,
- group-level membership matrix based on the configured prevalence threshold,
- core and unique taxon tables with the threshold stated in every output,
- an UpSet plot when at least two non-empty sets exist.

Do not define a group taxon as present merely because one read occurs in one sample unless the configured threshold explicitly says so. Do not call a taxon "core" without stating the sample/group prevalence rule. Venn plots are not required in Phase 1.

### 7.7 `07_kreport_pavian.R`

Generate one six-column Kraken-style report per selected sample. Validate:

- `U` unclassified row,
- `R` root row,
- standard rank codes `D`, `K`, `P`, `C`, `O`, `F`, `G`, `S`,
- integer counts,
- `unclassified + root_clade = total`,
- every node satisfies `clade = direct + sum(immediate child clades)`,
- percentages are calculated against total reads,
- indentation matches tree depth,
- output is deterministic.

The existing converter's `D1` kingdom rank is a known defect and must not be preserved.

Taxonomy resolution rules:

1. Reuse non-zero entries in the existing cache.
2. Derive leaf TaxIDs from exact normalized assignment lineages, not species-name-only matching. For the supplied minimap2 schema, normalize the eight-rank abundance path by omitting `kingdom` before matching its remaining seven ranks to the assignment lineage.
3. When multiple TaxIDs map to one exact lineage, record the conflict and deterministic majority choice.
4. In `cache_only` mode, never contact NCBI. Write an unresolved-node TSV and obey `unresolved_policy`.
5. In refresh mode, require `NCBI_EMAIL`; accept `NCBI_API_KEY`; include a fixed tool name; use HTTPS timeouts, bounded retries/backoff, and documented request pacing.
6. Merge refreshed entries into the cache and write through a temporary file followed by atomic replacement. Never destroy a valid cache after a partial network failure.
7. Check the Python subprocess exit code and the expected cache/report outputs before declaring success.

The resolver should keep a provenance sidecar containing refresh time, request mode, source, input hashes, and unresolved/conflicting nodes. Phase 1 documents uploading the resulting report to Pavian; it does not claim to generate a standalone Pavian HTML file.

### 7.8 Minimum output contract

Implementations may add useful artifacts, but these stable paths are required:

| Module | Required outputs |
|---|---|
| QC | `01_QC/classification_reconciliation.tsv`, `01_QC/read_length_by_status.tsv`, and per-sample PNGs under `01_QC/<SampleID>/` |
| Alpha | `02_Alpha_Diversity/alpha_diversity.tsv`, `rarefaction_curve.tsv`, `rarefaction_resamples.tsv`, and corresponding PNGs |
| Beta | Per-distance matrix, PCoA scores, and PCoA PNG; `permanova.tsv`; `betadisper.tsv`; or a structured skip record |
| Composition | Complete count and relative-abundance TSVs per analyzed rank plus single-sample bars or cohort stacks/heatmap as applicable |
| Ordination | PCA scores/loadings/variance TSVs and PNG; NMDS scores/diagnostics and PNG when its gate passes |
| Shared taxa | `sample_presence.tsv`, `group_prevalence.tsv`, `group_membership.tsv`, `core_taxa.tsv`, `unique_taxa.tsv`, and UpSet PNG when its gate passes |
| Kreport | `07_Kreport/<SampleID>.kreport`, `taxonomy_resolution.tsv`, and `unresolved_taxids.tsv` |

Tabular outputs must be tab-delimited UTF-8 with headers, stable row ordering, explicit sample/taxon identifiers, and no row names disguised as an unnamed first column.

---

## 8. Runner, dependency, and provenance behavior

### 8.1 Orchestration

`00_run_pipeline.R` must:

1. parse CLI arguments strictly and reject unknown/missing values,
2. locate its own script directory and repository root without assuming the caller's working directory,
3. load and validate config,
4. check dependencies before creating output,
5. build and validate the shared context once,
6. honor `--validate-only`,
7. derive output directories after CLI overrides,
8. run requested modules in documented order,
9. capture start/end time, status, outputs, warnings, and errors per module,
10. fail fast by default or continue only with `--keep-going`,
11. write a failed run manifest where possible,
12. exit non-zero if any requested module fails.

Optional/skipped modules are not failures only when their prerequisite is genuinely absent, such as QC without optional assignments or beta diversity in single mode. Every skip needs a reason in the manifest.

### 8.2 Dependencies

Keep one authoritative package list. Expected Phase 1 R dependencies are:

`yaml`, `optparse`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `scales`, `vegan`, `RColorBrewer`, `jsonlite`, `pheatmap`, `UpSetR`, `ggrepel`, `digest`, and `testthat` for tests.

The Python resolver should remain standard-library-only.

`analysis/install_packages.R` must be check-only by default and install only with an explicit `--install` flag. Missing or failed packages must cause a non-zero exit. Record installed package versions and `sessionInfo()` in each real run. This is dependency auditing, not a claim of a fully locked environment.

### 8.3 Manifest

`run_manifest.json` should include:

- pipeline/schema version,
- start/end timestamps and final status,
- command-line options,
- resolved config path and output root,
- selected sample IDs and mode,
- input absolute paths, sizes, modification times, and SHA-256 values,
- package and interpreter versions,
- run seed,
- module statuses, skip reasons, warnings, and output files,
- taxonomy cache mode and unresolved count.

Do not include NCBI API keys or other secrets.

---

## 9. Implementation sequence and gates

### Phase A: contracts and test scaffold

Implement CLI/config loading, path resolution, dependency checks, input parsers, context construction, and synthetic fixtures.

Gate A:

- all R files parse,
- Python compiles,
- config paths work from both repository root and another current directory,
- `--output-dir` relocates every derived module directory,
- `--validate-only` produces no filesystem changes,
- malformed fixture tests fail with specific messages,
- the real input produces all Section 2.2 invariants.

### Phase B: single-sample modules

Implement QC, alpha diversity, composition, and deterministic plotting from the shared context.

Gate B:

- real-input regression metrics match Section 2.3 within tolerance,
- fixed seed gives byte-stable tabular resampling output across repeated runs on the same environment,
- all expected tables and non-empty PNG files are produced in a temporary output root,
- no legacy tracked output is modified during tests,
- no network call occurs.

### Phase C: cohort modules

Implement cohort alpha summaries, beta diversity, ordination, composition heatmap, and shared taxa using a synthetic fixture with at least two groups and three samples per group.

Gate C:

- metadata is aligned by `SampleID`, never row position,
- sample permutation in input files does not change keyed numerical results,
- under-replicated fixtures skip inferential tests with explicit reasons,
- PERMANOVA and dispersion outputs are paired,
- single-sample runs skip cohort modules cleanly,
- rarefaction iterations never appear as biological sample rows.

### Phase D: kreport and taxonomy

Implement deterministic report generation and the hardened resolver.

Gate D:

- existing cache supports a completely offline real-input run,
- kingdom rows use `K`,
- report arithmetic/tree invariants pass,
- unresolved and conflicting TaxIDs are exported,
- a simulated resolver/network failure preserves the previous cache and returns non-zero,
- no refresh occurs without explicit opt-in.

### Phase E: integration and documentation

Update README, examples, `.gitignore`, and the run manifest. Run the full validation matrix.

Gate E:

- root invocation works: `Rscript analysis/00_run_pipeline.R --config config.yml`,
- invocation from another directory with an absolute config path also works,
- `--validate-only` is mutation-free,
- default failure behavior returns non-zero,
- `--keep-going` records all failures and still returns non-zero,
- README output claims match the generated manifest,
- `git diff --check` passes,
- only intended target-repository paths are changed.

Do not move to the next phase while its gate is failing.

---

## 10. Required validation commands

Run equivalents appropriate to the installed Windows environment:

```powershell
Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
python -c "import ast,pathlib; ast.parse(pathlib.Path(r'analysis/utils/ncbi_taxonomy.py').read_text(encoding='utf-8'))"
Rscript tests/testthat.R
Rscript analysis/00_run_pipeline.R --config config.yml --validate-only
Rscript analysis/00_run_pipeline.R --config config.yml --output-dir C:/tmp/ont-wf16s-phase1 --overwrite
git -c safe.directory='D:/W/AAy_Amplicon' diff --check
git -c safe.directory='D:/W/AAy_Amplicon' status --short
```

Tests must use a test-created temporary directory rather than a hard-coded shared path. The example `C:/tmp/...` is for a deliberate manual integration run only.

The test suite must include:

- config resolution and CLI precedence,
- output-root derivation,
- abundance schema and total reconciliation,
- assignment parsing and length variants,
- exact real-input count reconciliation,
- rank-prefix handling for contextual `Unknown` taxa,
- alpha regression values,
- mode/statistical gates,
- metadata order invariance,
- deterministic seed behavior,
- kreport structure/arithmetic/rank codes,
- cache-only behavior and atomic refresh failure,
- exit-code behavior.

Do not require internet access in the default test suite. A live NCBI smoke test, if implemented, must be separately named, opt-in, and excluded from normal CI/local validation.

---

## 11. Documentation and Git hygiene

README must document:

- the exact supported wf-16S abundance and minimap2 assignment schemas,
- that assignments are optional and per sample,
- single versus cohort prerequisites,
- classified-only versus all-read denominators,
- rarefaction terminology and limits,
- PERMANOVA/dispersion/replication gates,
- offline default and opt-in NCBI refresh,
- `.kreport` generation plus manual Pavian use,
- legacy script status,
- local validation versus genuine cohort/biological validation.

Append only this generated-output rule to `.gitignore` unless implementation creates another narrowly defined transient path:

```gitignore
# Modular pipeline outputs
/output/
```

Do not add broad `*.kreport`, `figs/`, or `tables/` rules in Phase 1. Do not remove tracked legacy artifacts. Never stage `ATW_Sesame_Greenhouse-Examples/`.

Before any future commit, inspect the exact diff and stage an explicit approved path list. Commit and push are outside this plan and require separate user authorization.

---

## 12. Definition of done for Gemini handoff

Implementation is complete only when all of the following are true:

- every promised file exists and contains executable code rather than placeholders,
- the runner is independent of the caller's current working directory,
- input and metadata contracts fail closed with actionable messages,
- supplied single-sample results reconcile to the audited counts and metrics,
- cohort behavior is covered by synthetic replicated fixtures and is described as not yet biologically validated,
- single-sample PCA/NMDS/PERMANOVA/shared-taxa analyses are skipped rather than fabricated,
- rarefaction resamples are never presented as biological replication,
- kreport uses standard rank codes and passes arithmetic/tree validation,
- default execution and tests make no network request,
- taxonomy refresh is explicit, policy-compliant, retry-bounded, and cache-safe,
- package-install failures and module failures return non-zero,
- run provenance and skip/failure status are recorded,
- README and example config match actual CLI behavior,
- legacy tracked data and unrelated untracked content remain untouched,
- validation commands and `git diff --check` pass,
- no commit or push has been performed without explicit approval.

Useful upstream contracts for implementation review:

- [EPI2ME wf-16S repository and output documentation](https://github.com/epi2me-labs/wf-16S)
- [wf-16S machine-readable output definition](https://github.com/epi2me-labs/wf-16s/blob/master/output_definition.json)
- [Kraken 2 report format](https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format)
- [NCBI E-utilities parameters and policy](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
