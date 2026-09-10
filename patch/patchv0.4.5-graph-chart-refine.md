# `wf16s-postprocess` graph and chart refinement plan (provisional v0.4.5)

**Implementation baseline:** `main` / `origin/main` at `5be37b3`, with `VERSION` `0.4.3`, audited 2026-09-09.

**Requested release target:** `0.4.5`. The repository currently has only a `v0.4.3` tag; complete the release-number gate in Section 10 before editing release metadata.

**Document status:** source-audited Gemini handoff plan. This document does not authorize implementation commits, tags, or pushes.

## 1. Outcome

Deliver three related improvements without changing taxonomy, count, or publication semantics:

1. Produce deterministic, vertically stacked, classified-read-relative composition plots at phylum, family, and genus ranks for both single-sample and cohort runs.
2. Produce top-10 phylum, family, and genus heatmaps in both modes, preserving experimental group/sample order while clustering taxa only when that is informative.
3. Make opt-in Krona-compatible HTML rendering work without an external Perl/KronaTools installation by adding a standard-library Python renderer and pinned, unmodified upstream browser assets.

Every new visual must have an inspectable TSV sidecar, every emitted artifact must be declared by its module, and a requested HTML render must either succeed or fail the module. Do not silently publish a TSV-only result when `render_html: true`.

## 2. Evidence from the current checkout

The implementation must be based on the current contracts, not only on the reference images.

- `analysis/04_taxa_composition.R` already calculates six rank tables using the classified-read denominator and writes `count_<rank>.tsv` plus `rel_abundance_<rank>.tsv`.
- Single-sample mode currently writes four horizontal plots and no heatmap.
- Cohort mode currently writes only `04_phylum_stacked.png` and one configurable heatmap. It does not join `Group` onto the stacked-bar data, and it clusters heatmap columns.
- Cohort mode always has validated `SampleID` and `Group` metadata. `read_metadata_table()` aligns metadata to `context$samples`, so plot code must make any further ordering rule explicit.
- The `kreport` module already requires Python for `ncbi_taxonomy.py`; a standard-library Python HTML builder adds no new system-language dependency to that module.
- Current Krona behavior is opt-in. `--krona` sets `krona.enabled: true`; `render_html` defaults to true; a missing `ktImportText` currently causes a warning and a TSV-only result.
- Module code runs under a private per-module staging directory. `all_outputs` is later used to build ownership, artifact hashes, and the physical file census.
- `config.yml` is a tracked one-sample regression configuration, not a user scratch file. Krona must remain disabled there unless explicitly requested by CLI or YAML.
- The local `ATW_Sesame_Greenhouse-Examples/output_sensitivity/Barplot` and `Heatmap` images are untracked style references only. Do not stage them, read analytical values from them, or make tests depend on them.

The official Krona documentation defines `ktImportText` input as a quantity followed by tab-separated hierarchy labels and documents the Krona 2.0 XML structure. KronaTools 2.8.1 produces standalone charts by default. If its browser assets are redistributed, the upstream copyright/license notice must be retained; altered code must not be presented as official Krona software.

Primary upstream references:

- [Krona text and XML import documentation](https://github-wiki-see.page/m/marbl/Krona/wiki/Importing-text-and-XML-data)
- [Krona 2.0 XML specification](https://github-wiki-see.page/m/marbl/Krona/wiki/Krona-2.0-XML-Specification)
- [KronaTools 2.8.1 release](https://github.com/marbl/Krona/releases/tag/v2.8.1)
- [KronaTools 2.8.1 license](https://raw.githubusercontent.com/marbl/Krona/v2.8.1/KronaTools/LICENSE.txt)

## 3. Scope boundaries

### In scope

- Module 04 composition preparation, static plotting, sidecars, and focused helpers.
- Configuration migration and validation for multi-rank figures.
- A built-in, offline, self-contained Krona-compatible HTML renderer.
- Vendored upstream assets, license/source provenance, and their inclusion in maintained-source hashing.
- Unit, integration, release-verifier, documentation, changelog, and release-metadata updates.
- Cross-platform Windows and Ubuntu CI coverage.

### Out of scope

- Changes to upstream input schemas, taxonomy resolution, TaxID conflict handling, read classification, alpha/beta diversity, ordination, FAPROTAX, or Pavian output.
- Adding unclassified reads to classified composition bars or heatmaps. Unclassified reads remain in `classification_fraction.tsv` and `04_classification_fraction.png`.
- Reproducing BGI/Sesame values, exact fonts, or pixel-identical reference images.
- Runtime downloads, CDNs, an embedded web server, or a new JavaScript build toolchain.
- Modifying legacy root scripts.
- Committing any untracked producer outputs, reference bundles, temporary directories, or this plan unless explicitly included in a later approved delivery.

## 4. Non-negotiable contracts

### 4.1 Scientific accounting

For taxon `T` and sample `s` with positive classified depth:

```text
RelativeAbundance(T, s) = Count(T, s) / ClassifiedReads(s)
```

- Existing `count_<rank>.tsv` and `rel_abundance_<rank>.tsv` schemas and values remain unchanged.
- New composition bars and heatmaps exclude the all-read `Unclassified` row.
- For a valid sample, selected taxa plus `Other` must sum to 1 within floating-point tolerance (`1e-10`).
- A sample with `ClassifiedReads == 0` is not silently converted to a 100% bar. Its new sidecar rows must mark `ValidDenominator = FALSE`, and its plotted bar is zero height.
- Cohort taxon selection and group means exclude zero-classified samples from mean-relative-abundance calculations. Sidecars report both total samples and samples used. If no sample at a rank has a valid positive abundance, write a structured `*_skipped.tsv` instead of a misleading plot.
- Group bars are arithmetic means of per-sample relative abundances, not pooled read-weighted proportions:

```text
GroupMean(T, g) = mean(RelativeAbundance(T, s) for valid s in group g)
```

### 4.2 Taxon identity and labels

- `TaxonPath` remains the analytical key throughout selection, aggregation, joins, and tests.
- Never group by a display name alone.
- Derive the display leaf from the final `TaxonPath` component.
- If leaf names collide, append the immediate parent; if that is still ambiguous, append the shortest unique path suffix. The mapping must be deterministic.
- Keep `Other` as a reserved synthetic key (for example `TaxonPath = "__OTHER__"`) and never allow a biological display label to collide with it.
- Unknown classified rank labels remain classified data and must be contextualized. They are not relabeled as all-read `Unclassified`.

### 4.3 Determinism and ordering

- Rank order follows the configured vector exactly.
- Candidate taxa are ordered by descending mean relative abundance, then `TaxonPath` ascending as the tie-break.
- Cohort group order is stable first appearance in metadata aligned to `context$samples`.
- Samples are grouped by that group order while retaining their original `context$samples` order within each group. Do not alphabetically sort IDs; lexical sorting misorders identifiers such as `S2` and `S10` and can discard experimental intent.
- The same selected-taxon order and color mapping are reused by sample-level bars, group-mean bars, legends, and sidecars for a rank.
- Generated PNG and built-in HTML contents must not contain timestamps, random IDs, absolute staging paths, or network-dependent content.

### 4.4 Transaction, manifest, and provenance behavior

- Modules write only beneath their private staging roots.
- Every file actually written is returned once in `outputs`; no nonexistent or directory paths are returned.
- A module failure leaves no partially published final output. The orchestrator owns rollback/publication.
- Physical files must continue to equal declared owned files (excluding the manifest itself) plus preserved unowned files.
- Vendored JavaScript, images, license, and source metadata must participate in `source_digest_sha256`; do not broaden the source allowlist to arbitrary binary files elsewhere in the checkout.
- `--validate-only` remains mutation-free, including built-in-renderer validation.

## 5. Composition configuration and compatibility

Use additive, validated configuration. Retain the existing `top_n_taxa` setting for the legacy horizontal single-sample plots.

```yaml
composition:
  top_n_taxa: 15

  stacked_bar_ranks: ["phylum", "family", "genus"]
  stacked_bar_max_taxa: 15
  stacked_bar_min_taxa: 5
  stacked_bar_min_mean_relative: 0.005

  heatmap_rank: null  # deprecated scalar compatibility key
  heatmap_ranks: ["phylum", "family", "genus"]
  heatmap_top_n_taxa: 10
  heatmap_include_other: true
  heatmap_transform: "log10_relative"
```

Migration rules must be implemented against the raw YAML before defaults are merged:

1. If only legacy `composition.heatmap_rank` is explicitly supplied and is non-null, warn once and treat it as the complete one-element `heatmap_ranks` vector.
2. If both a non-null legacy key and `heatmap_ranks` are explicitly supplied, fail with a configuration error; do not guess precedence.
3. The tracked `config.yml` and `config.example.yml` move to `heatmap_rank: null` plus the three-rank vector, so normal regression runs exercise the new defaults without a deprecation warning.
4. Keep `schema_version: 1` only if old scalar configurations pass the migration tests unchanged apart from the documented warning. Otherwise bump the config schema and add an explicit migration note.

Validation requirements in `analysis/utils/config.R`:

- Rank vectors are non-empty, unique character vectors drawn from `RANKS_TO_ANALYZE` (or one shared supported-rank constant moved to a utility file).
- `stacked_bar_max_taxa` and `stacked_bar_min_taxa` are positive integers, with minimum less than or equal to maximum.
- `stacked_bar_min_mean_relative` is finite and in `[0, 1]`.
- `heatmap_top_n_taxa` is a positive integer.
- `heatmap_include_other` is one non-missing boolean.
- Keep existing `heatmap_transform` values: `log10_relative` and `none`.
- Unknown keys continue to fail closed.

## 6. Module 04 implementation

### 6.1 Refactor preparation into testable helpers

Keep rank aggregation in `analysis/04_taxa_composition.R`; add only generic theme/palette helpers to `analysis/utils/plotting.R`. Suggested function boundaries:

```r
make_taxon_display_map(rank_table)
select_display_taxa(rel_table, samples, valid_samples, min_mean, min_n, max_n)
collapse_rank_for_display(rel_table, selected_paths, samples, metadata)
summarize_group_means(display_long)
composition_colors(display_paths, display_labels)
build_stacked_taxa_plot(display_long, colors, sample_order, group_order, single)
prepare_heatmap_data(rel_table, samples, metadata, top_n, include_other, transform)
draw_taxa_heatmap(prepared, path, mode)
```

Helpers should return data/plot objects; file writes stay in `run_taxa_composition()`. This permits semantic tests without image comparison.

### 6.2 Exact stacked-taxon selection

For each configured stacked-bar rank:

1. Start from the existing classified-only rank relative-abundance table.
2. Remove taxa with zero count across all samples from the candidate set.
3. Compute each candidate's mean across valid (`ClassifiedReads > 0`) samples.
4. Sort by descending mean, then ascending `TaxonPath`.
5. Select candidates at or above `stacked_bar_min_mean_relative`, capped at `stacked_bar_max_taxa`.
6. If fewer than `min(stacked_bar_min_taxa, number_of_positive_candidates)` were selected, fill from the remaining sorted candidates to that number.
7. Collapse every unselected positive candidate into one `Other` row per sample. Omit `Other` only when its value is zero for every sample.

The long sample sidecar is `04_<rank>_stacked.tsv` with this fixed schema:

```text
Rank, TaxonPath, DisplayTaxon, SampleID, Group,
RelativeAbundance, MeanRelativeAbundance, IsOther,
ValidDenominator, StackOrder, SampleOrder
```

Use tab delimiters and the repository's existing `write.table(..., quote = FALSE, row.names = FALSE)` convention. Assert unique `(SampleID, TaxonPath)` rows, finite values in `[0, 1]`, and exact sample membership before writing.

### 6.3 Sample-level stacked figures

Use the same output names in both modes:

- `04_phylum_stacked.png`
- `04_family_stacked.png`
- `04_genus_stacked.png`

This preserves the existing cohort phylum filename and avoids a second naming scheme for single-sample output. Keep the existing horizontal files unchanged:

- `04a_phylum_composition.png`
- `04b_family_composition.png`
- `04c_genus_composition.png`
- `04d_species_composition.png`

Plot contract:

- `geom_col(width = 0.7)` using already normalized values; do not use `position = "fill"` because it would turn zero-classified samples into artificial full-height bars.
- Y limits `[0, 1]`, breaks at `0, 0.25, 0.5, 0.75, 1`, percent labels, and zero expansion at the baseline.
- Largest selected category at the bottom, then descending selected taxa, with `Other` at the top. Verify this using `ggplot_build()` in tests rather than relying on an untested factor-level assumption.
- In cohort mode, join validated metadata and use `facet_grid(cols = vars(Group), scales = "free_x", space = "free_x")` with the deterministic group/sample factors defined in Section 4.3.
- Rotate cohort sample labels 90 degrees. A single-sample label remains horizontal.
- Add in-slice labels only in single-sample mode and only for values at least `0.03`; do not label zero or invalid-denominator slices.
- New composition figures render at 300 DPI. Suggested size: single `7 x 7` inches; cohort width `max(10, min(24, 6 + 0.32 * n_samples))`, height `8`.
- Keep legends on the right and allow two columns only when needed.

Use an explicit 15-color constant, not a vague 35-color palette when at most 15 named taxa are displayed. Assign colors in selected-taxon order and reuse the mapping for all figures at that rank. Reserve `#F3C58F` for `Other` and `#BDBDBD` only for the separate all-read `Unclassified` status. Candidate named colors:

```r
composition_palette <- c(
  "#FF8C5A", "#E6B66E", "#8DB17E", "#5C9462", "#2F6B34",
  "#34C79B", "#63C4C7", "#3699C5", "#5F65C8", "#9060C2",
  "#BF5BC4", "#ED54B8", "#F184B9", "#C69090", "#996666"
)
```

### 6.4 Cohort group-mean figures

Cohort metadata is mandatory, so every successful cohort composition run emits:

- `04_<rank>_group_mean_stacked.png`
- `04_<rank>_group_mean_stacked.tsv`

The fixed TSV schema is:

```text
Rank, Group, TaxonPath, DisplayTaxon, MeanRelativeAbundance,
SamplesTotal, SamplesUsed, SamplesExcludedZeroClassified,
IsOther, StackOrder, GroupOrder
```

Use the identical selected taxa, `Other`, colors, and stack order as the sample-level figure. Group width is `max(7, min(18, 3 + 1.25 * n_groups))`; height is `7` inches at 300 DPI. Do not label these means as pooled abundance.

### 6.5 Multi-rank heatmaps

For every configured heatmap rank:

1. Select the top `heatmap_top_n_taxa` positive named taxa by descending mean across valid samples, breaking ties by `TaxonPath`.
2. If `heatmap_include_other: true`, append one `Other` row only when the residual is positive in at least one sample. "Top 10" therefore means ten named taxa, optionally plus one clearly marked summary row.
3. Exclude all-read `Unclassified` from the matrix.
4. Derive one pseudo-count per rank from the displayed matrix:

```r
positive <- matrix_values[matrix_values > 0]
pseudo_count <- min(positive) / 2
transformed <- log10(matrix_values + pseudo_count)
```

If there are no positive values, emit `04_heatmap_<rank>_skipped.tsv` with `Rank` and `Reason` and do not fabricate a pseudo-count or PNG.

5. For `heatmap_transform: none`, do not add or report a pseudo-count.
6. Use `rev(RColorBrewer::brewer.pal(11, "RdYlBu"))` expanded to 100 colors, so higher abundance is red and lower abundance is blue. This is a display scale, not evidence of enrichment or depletion.
7. Use 101 finite monotonic breaks spanning the transformed display range. For a constant matrix, expand the range by `0.5` on each side.
8. Set `cluster_cols = FALSE` in all modes. Order columns by Section 4.3 and use a `Group` annotation in cohort mode.
9. In cohort mode, precompute complete-linkage Euclidean `hclust` when at least two display rows exist and pass that object to `pheatmap`. In single-sample mode, set `cluster_rows = FALSE` and sort rows by descending raw abundance; a one-column dendrogram adds no comparative information.
10. Use `border_color = "#D9D9D9"`, 90-degree cohort column labels, and a title that states rank, transform, and the classified-read denominator.

Outputs:

- `04_heatmap_phylum.png` and `04_heatmap_phylum.tsv`
- `04_heatmap_family.png` and `04_heatmap_family.tsv`
- `04_heatmap_genus.png` and `04_heatmap_genus.tsv`

The long heatmap TSV schema is:

```text
Rank, TaxonPath, DisplayTaxon, SampleID, Group,
RelativeAbundance, Transform, PseudoCount, TransformedValue,
IsOther, RowOrder, ColumnOrder, ValidDenominator
```

For a clustered cohort heatmap, `RowOrder` records the actual dendrogram display order. The PNG must be produced through `with_png_device()` or an equivalently cleanup-safe device wrapper. Suggested sizes: single `6 x 7` inches; cohort width `max(8, min(24, 5 + 0.28 * n_samples))`, height `max(6, min(12, 3.5 + 0.45 * n_display_rows))`, at 300 DPI.

## 7. Built-in Krona-compatible HTML renderer

### 7.1 Reproducibility and licensing design

Do not paste a minified upstream script into a Python string. Vendor the minimum unmodified browser resource set under a versioned directory, for example:

```text
analysis/vendor/krona-2.8.1/
  LICENSE.txt
  SOURCE.json
  src/krona-2.0.js
  img/hidden.png
  img/loading.gif
  img/logo-med.png
  img/favicon.ico
```

Before finalizing the asset list, inspect the pinned 2.8.1 standalone HTML generator and include only resources it actually embeds. `SOURCE.json` must record the upstream tag, source URL, retrieval date, and SHA-256 for every vendored file. Preserve upstream files byte-for-byte. The local builder and documentation should say "Krona-compatible" and include attribution; do not imply that the modified wrapper is official Krona software.

Add a repository test that compares every vendored file to its declared SHA-256. No test or runtime path downloads these resources.

Update `maintained_source_files()` so every regular file under this exact vendor directory participates in source provenance. Add a test that mutating a copied `.js` or image asset changes the computed source digest. Do not globally admit every image in `analysis/` as maintained source.

### 7.2 Python builder contract

Add `analysis/utils/krona_builder.py`, using only the Python standard library. CLI:

```text
python krona_builder.py \
  --input <sample.krona.tsv> \
  --output <sample.krona.html> \
  --dataset-name <SampleID> \
  --expected-total <TotalReads> \
  --vendor-dir <analysis/vendor/krona-2.8.1>

python krona_builder.py --validate-only --vendor-dir <vendor-dir>
```

Required behavior:

- Parse UTF-8 tab-delimited input without shell interpretation.
- Require a canonical non-negative integer magnitude and at least one non-empty hierarchy label per line.
- Reject control characters, malformed rows, negative/fractional counts, and totals that differ from `--expected-total`.
- Aggregate duplicate hierarchy paths deterministically.
- Track direct magnitude separately from recursively calculated clade magnitude. Each XML node's displayed magnitude is its clade total; the root equals `TotalReads`.
- Serialize valid Krona 2.0 XML using `xml.etree.ElementTree` or equivalent escaping-safe standard-library code. Escape dataset and taxon labels correctly, including non-ASCII text.
- Embed the XML and pinned browser resources into one HTML file. The finished HTML must have no external script, stylesheet, image, font, fetch/XHR, or CDN requirement. Ordinary clickable attribution links are allowed.
- Write a sibling temporary file, flush/close it, and replace the target atomically. Remove the temporary file on failure.
- Emit no timestamps or absolute local paths. Identical input, sample ID, expected total, and vendor assets must produce identical bytes.
- `--validate-only` verifies the source manifest, asset hashes, and builder prerequisites without creating files.
- Return non-zero with a concise stderr message on every contract failure.

### 7.3 Renderer selection

Replace the ambiguous `fallback_to_builtin` boolean with an explicit renderer policy:

```yaml
krona:
  enabled: false
  render_html: true
  html_renderer: "builtin"  # builtin, kronatools, or auto
  executable: "ktImportText"
```

Semantics:

- `builtin` (new default): use `krona_builder.py`; no Perl/KronaTools installation is required.
- `kronatools`: require and use the configured executable. Missing or failed external rendering is an error.
- `auto`: use `ktImportText` when available; otherwise use the built-in renderer. This explicitly selected mode may yield different HTML bytes on hosts with different tools, and provenance must say which provider was used.
- `render_html: false`: write the validated TSV and provenance only, preserving the supported TSV-only workflow.
- `enabled: false`: emit no Krona directory or provenance, as today.

Keep `enabled: false` in defaults, `config.example.yml`, and tracked `config.yml`. `--krona` continues to enable the extension and, with the default renderer policy, now produces HTML on standard CI hosts.

### 7.4 R integration and fail-closed behavior

In `analysis/utils/kreport.R`:

- Keep `write_krona_input()` and its direct-count accounting unchanged.
- Rename the current external wrapper to `render_kronatools_html()`.
- Add `render_builtin_krona_html()` that invokes Python through `processx::run()` with a true argument vector, reusing the already resolved `python_cmd` from `run_kreport()`.
- Validate status, target existence, non-zero size, and absence of a leftover temporary output.
- Add one resolver function returning the chosen provider and executable/builder details from `html_renderer`.

In `analysis/07_kreport_pavian.R`:

- Resolve the provider once before the sample loop.
- Render each sample immediately after its validated TSV is written.
- If HTML was requested and any render fails, return/raise a module failure. The private module staging directory will then be discarded, so a partial mixed-provider result is never published.
- Add every successful HTML file to `all_outputs`.
- Retain exact total/classified/unclassified/emitted-magnitude checks.

Update `krona_provenance.json` without losing stable array shapes. Required top-level fields:

```text
format, renderer_policy, renderer, renderer_version,
vendor_source, vendor_sha256_manifest, html_status,
standalone_html, count_model, denominator, classified_definition,
requested_executable, resolved_executable, render_html, samples
```

Rules:

- `html_status` is `rendered` or `not_requested` in a successfully published new run. Retain parser compatibility with historical `renderer_missing` provenance if old outputs are inspected, but do not emit it for a successful `render_html: true` run.
- `renderer` is `builtin_krona_compatible` or `ktImportText` when rendered, otherwise null.
- `renderer_version` is the local builder version plus vendored Krona asset version for built-in rendering, or the probed KronaTools version for external rendering.
- `standalone_html` is true only after a successful self-contained HTML write.
- A sample record has `html_path` exactly when top-level status is `rendered`.

### 7.5 Preflight

Update `analysis/utils/preflight.R`:

- `builtin`: AST-parse the builder, call its mutation-free `--validate-only`, and verify the pinned vendor asset manifest.
- `kronatools`: fail `E_KRONA_PREFLIGHT` if the configured executable is missing, a directory, or not runnable.
- `auto`: validate the external executable if found; otherwise validate the built-in provider. Do not warn that HTML will be skipped when the built-in provider is valid.
- Skip renderer checks when `render_html: false`.
- Preserve the existing `kreport` taxonomy resolver preflight and no-write guarantees.

## 8. File inventory

### Modify

- `analysis/04_taxa_composition.R`
- `analysis/07_kreport_pavian.R`
- `analysis/utils/config.R`
- `analysis/utils/kreport.R`
- `analysis/utils/plotting.R`
- `analysis/utils/preflight.R`
- `config.example.yml`
- `config.yml` (composition defaults only; do not enable Krona)
- `tests/testthat/test-config.R`
- `tests/testthat/test-cohort.R` (retain the existing one-row regression)
- `tests/testthat/test-kreport.R`
- `tests/testthat/test-plotting.R`
- `tests/testthat/test-provenance.R`
- `tests/verify_release_run.R`
- `.github/workflows/ci.yml`
- `README.md`
- `CHANGELOG.md`
- `VERSION` (only after the release-number gate)
- `CITATION.cff` (version and actual release date, only after the release-number gate)

### Add

- `analysis/utils/krona_builder.py`
- `analysis/vendor/krona-2.8.1/` minimal pinned asset/license/source bundle
- `tests/test_krona_builder.py`
- `tests/testthat/test-composition.R` for focused data/plot contracts

Do not add a separate single-sample test file solely for naming symmetry; the focused composition test file should cover both modes.

## 9. Test and verification matrix

### 9.1 Configuration tests

- New defaults and exact rank order.
- Invalid/duplicate/empty rank vectors.
- Integer/range/cross-field validation.
- Legacy scalar `heatmap_rank` migration and one warning.
- Conflict when legacy and new heatmap keys are both explicit.
- Krona renderer enum validation.
- Default Krona remains disabled; `--krona` enables it without changing unrelated module defaults.

### 9.2 Composition semantic tests

Use synthetic data with ties, duplicate leaf names under different parents, rare taxa, an `Other` remainder, metadata order that differs from alphabetical order, and at least one zero-classified sample.

Assert:

- Selection threshold/fill/cap and `TaxonPath` tie-break.
- Unique contextual display labels without identity collapse.
- Per-sample conservation and no values outside `[0, 1]`.
- Zero-classified sample remains a zero-height/invalid-denominator record.
- Group means use only valid samples and report sample counts.
- Stable group/sample factor order.
- `Other` is last/top and has the reserved color.
- Plot objects contain the cohort facet and configured percent scale.
- `ggplot_build()` confirms bottom-to-top stack order.
- Sidecar schemas, keys, and order columns are exact.
- All expected single/cohort PNGs and TSVs exist, are non-empty, and are returned in `outputs`.
- Existing horizontal single-sample plots and six base rank tables still exist.

### 9.3 Heatmap tests

- Three default ranks in both modes.
- Exactly top 10 named taxa, or fewer when fewer positive taxa exist; optional `Other` is separate.
- Data-derived pseudo-count and transformed values are reproducible.
- `none` transform has no pseudo-count.
- Cohort columns follow group plus original sample order and are not clustered.
- Cohort rows use the recorded complete-linkage dendrogram order.
- Single-sample rows are abundance-sorted and not clustered.
- One-row and constant-matrix cases render without invalid breaks.
- All-zero classified data emits structured skip TSVs rather than fabricated heatmaps.

### 9.4 Krona builder tests

Python tests must cover:

- Single and nested paths, duplicate-path aggregation, direct-plus-child clade totals, and explicit `Unclassified`.
- Expected-total conservation.
- XML/HTML escaping for `<`, `>`, `&`, quotes, and non-ASCII labels.
- Rejection of malformed, blank, fractional, negative, and control-character input.
- Deterministic byte-identical output across two builds.
- No external resource references.
- Pinned asset SHA-256 validation and mismatch failure.
- Mutation-free `--validate-only`.
- Atomic cleanup after an injected failure.

R tests must cover:

- `builtin`, `kronatools`, and `auto` provider resolution.
- Paths and sample IDs containing spaces.
- Missing explicit external renderer fails preflight.
- Missing external renderer in `auto` selects built-in without a TSV-only warning.
- Non-zero builder exit, missing output, and empty output fail the module.
- Provenance fields and JSON array shapes for one and multiple samples.
- HTML path presence agrees with `html_status`.

### 9.5 Provenance, transaction, and release checks

- Vendored non-R/Python assets are included in maintained-source hashing.
- `--validate-only` creates no output root or renderer temporary files.
- Failed built-in rendering discards the private module stage and preserves a previously published output root during overwrite testing.
- Physical census and per-artifact SHA-256 checks pass with every new PNG, TSV, HTML, license/source asset input, and provenance file handled correctly.
- `tests/verify_release_run.R` derives expected composition artifacts from `resolved_config.yml`, then validates their fixed schemas and numeric conservation.
- When `manifest$cli$krona` and `render_html` are true, the verifier requires exactly one non-empty HTML per sample and rejects historical `renderer_missing` as a successful new-run state.
- Ambar regression anchors remain `TotalReads = 114056`, `ClassifiedReads = 80556`, `UnclassifiedReads = 33500`, 46 unresolved taxonomy nodes, and 26 conflicts.

### 9.6 Commands

Run from the repository root in the locked environment:

```bash
Rscript analysis/install_packages.R --restore
Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
python -m compileall -q analysis tests
Rscript tests/testthat.R
python -m unittest -v tests/test_ncbi_taxonomy.py tests/test_krona_builder.py tests/test_check_committed_whitespace.py
```

Mutation-free validation:

```bash
Rscript analysis/00_run_pipeline.R --config config.yml --output-dir "validation only output" --validate-only --krona
```

Before and after this command, assert that `validation only output` does not exist.

Dirty-worktree smoke integration during implementation:

```bash
Rscript analysis/00_run_pipeline.R --config config.yml --output-dir "v045 worktree output" --krona --allow-dirty
```

Use unit/integration assertions to inspect this development run. Do not claim that `tests/verify_release_run.R` passed against an `--allow-dirty` run; that verifier intentionally requires clean exact-commit provenance.

Run the existing synthetic minimap2 smoke path during implementation and add a six-sample, two-group composition integration fixture or harness that exercises faceting and group means. It must be small, deterministic, and tracked; do not use the large untracked Sesame directory as CI input.

```bash
Rscript analysis/00_run_pipeline.R --config tests/fixtures/synthetic_minimap2/config.yml --output-dir "v045 synthetic worktree output" --allow-dirty
```

After commit authority is granted and the scoped implementation is committed, run canonical release integrations from the clean exact commit, without development escape hatches:

```bash
Rscript analysis/00_run_pipeline.R --config config.yml --output-dir "v045 single release output" --krona
Rscript tests/verify_release_run.R "v045 single release output"
Rscript analysis/00_run_pipeline.R --config tests/fixtures/synthetic_minimap2/config.yml --output-dir "v045 synthetic release output"
Rscript tests/verify_release_run.R "v045 synthetic release output"
```

Final static checks:

```bash
python tests/check_committed_whitespace.py
git diff --check
git status --short
```

Update `.github/workflows/ci.yml` so both Ubuntu and Windows run `tests/test_krona_builder.py` and at least one `--krona` integration using the default built-in provider. Do not make CI depend on an installed `ktImportText`.

### 9.7 Manual visual/browser QA

Automated file-existence tests do not establish visual quality or interactivity. Before sign-off:

1. Render and inspect all new single-sample figures from the Ambar fixture.
2. Render and inspect the small cohort fixture at phylum, family, and genus ranks.
3. Confirm readable labels, correct group order, no clipped legends/titles, `Other` at the top, and no artificial 100% bar for an invalid denominator.
4. Open a built-in HTML file directly from disk in current Chrome or Edge and Firefox. Confirm zoom, focus, back navigation, labels, total magnitude, and offline operation with networking disabled.
5. Treat the untracked Sesame images as visual inspiration only. Record any deliberate deviations, especially exclusion of all-read `Unclassified` from classified composition figures.

## 10. Implementation sequence and release gate

### Phase A: baseline and scope lock

1. Re-run `git status --short --branch` with the checkout-local safe-directory override.
2. Record the starting SHA and explicit approved file list.
3. Run the baseline test suite before editing. Distinguish pre-existing failures from regressions.
4. Leave all unrelated untracked directories and plans untouched.
5. Use `--allow-dirty` only for worktree smoke runs and record that exception; reserve the release verifier for a clean committed tree.

### Phase B: composition data contracts first

1. Add configuration migration/validation tests.
2. Implement pure selection, collapse, label, order, and group-mean helpers.
3. Add semantic tests and TSV sidecars.
4. Add plots only after the data contracts pass.
5. Run single/cohort tests and visually inspect the rendered PNGs.

### Phase C: built-in renderer and provenance

1. Pin and document the minimum upstream Krona 2.8.1 assets and license.
2. Add asset/source-digest tests.
3. Implement and unit-test the Python builder.
4. Add renderer policy, preflight, R invocation, fail-closed behavior, and provenance.
5. Add Windows/Ubuntu CI and offline browser acceptance.

### Phase D: documentation and release metadata

Update README sections for composition outputs, config migration, opt-in Krona behavior, renderer policy, offline guarantees, vendored attribution, and release scope. Update `[Unreleased]` during implementation.

Before changing `VERSION`, `CITATION.cff`, or creating a tag, resolve why the requested target is `0.4.5` when no `v0.4.4` tag exists. Either:

- document that `0.4.4` was intentionally reserved/skipped and proceed with `0.4.5`, or
- retarget this patch consistently to `0.4.4`.

Never silently mix version numbers. At release time, synchronize:

- `VERSION` with LF termination;
- `CITATION.cff` version and actual release date;
- README current-version and release-scope statements;
- a dated changelog section;
- release metadata tests.

### Phase E: delivery

Suggested reviewable commits, if commit authority is later granted:

1. `feat: add deterministic multi-rank composition figures`
2. `feat: add built-in offline Krona-compatible renderer`
3. `docs: document visualization and renderer contracts`
4. `chore: prepare v0.4.5 release` (only after the version gate)

Stage only the approved paths. Do not push, tag, or publish a release without explicit authority. If publication is authorized, verify the remote branch SHA and both operating-system CI results after the push.

## 11. Acceptance checklist

- [ ] Existing classified-read denominators and six base rank tables are unchanged.
- [ ] Existing horizontal single-sample plots remain available.
- [ ] Phylum, family, and genus stacked PNG/TSV pairs are produced in single and cohort modes.
- [ ] Cohort sample plots are faceted in stable group order; group-mean plots state arithmetic-mean semantics.
- [ ] Taxon selection, labels, colors, stack order, group order, and sample order are deterministic and tested.
- [ ] Zero-classified and all-zero edge cases are represented honestly.
- [ ] Phylum, family, and genus heatmap PNG/TSV pairs work for one and multiple samples.
- [ ] Heatmap columns are never clustered; cohort row clustering and single-sample row sorting follow the documented rules.
- [ ] Every visual has a complete, finite, inspectable data sidecar.
- [ ] Krona remains opt-in, but default requested HTML no longer depends on Perl/KronaTools.
- [ ] Built-in HTML is byte-deterministic, self-contained, offline, and interactive in direct-file browser QA.
- [ ] Vendored upstream assets are minimal, pinned, hash-verified, attributed, license-compliant, and included in source provenance.
- [ ] Requested HTML failure fails the kreport module and publishes no partial result.
- [ ] `--validate-only` remains mutation-free.
- [ ] Manifest schema, artifact hashes, ownership, overwrite rollback, and physical census checks pass.
- [ ] R/Python tests and Windows/Ubuntu integrations pass in the locked environment.
- [ ] README, changelog, VERSION, and CITATION metadata agree after the release-number decision.
- [ ] Unrelated untracked reference and producer-output directories remain untouched and unstaged.
