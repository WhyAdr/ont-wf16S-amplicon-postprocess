# Patch plan: `ont-wf16S-amplicon-postprocess` v0.2.0

**Handoff target:** Gemini
**Repository:** `WhyAdr/ont-wf16S-amplicon-postprocess`
**Base version:** v0.1.0 (pending P0-1 composition fix)
**State review:** `ont-wf16s-postprocess-v0.1.0-statereview.md`
**Scope:** Upfront read accounting, classification diagnostics, richness overview,
composition fix, upstream contract enforcement, and expanded integration tests

## 1. Executive summary

v0.2.0 adds three new front-matter diagnostic outputs to every pipeline run:

1. **Read accounting table and pie chart** — a clear, exact breakdown of every
   read into three mutually exclusive categories (`status=U`, `C + TaxID 0`,
   `TaxID > 0`), with a pie chart and annotated legends.
2. **Investigative columns** — median lengths and minimap2 threshold-failure
   partitions (identity-only, coverage-only, both) for each dataset.
3. **Species richness overview** — positive taxa, singletons, taxa with ≤10
   reads, and reads in those taxa.

All existing modules retain their current hardcoded unclassified-read handling:
classified-only for diversity, composition, and ordination; included in the
kreport and the classification-fraction chart.

This patch also fixes the P0-1 single-sample composition bug
(state review §7 P0-1) and adds the upstream contract enforcement
(state review §7 P0-2, P1-1).

---

## 2. Feature specifications

### F1 — Read accounting table and pie chart

**Addresses:** state review §4.2 (exact read accounting), user feature request #1

#### Output: `00_read_accounting.tsv`

Written to the QC output directory. One row per sample. Columns:

| Column | Description |
|---|---|
| `SampleID` | Sample identifier |
| `AbundanceTotal` | Total reads from the abundance table |
| `RawC` | Assignment rows with `status=C` |
| `RawU` | Assignment rows with `status=U` |
| `C_TaxID0` | Reads with `status=C` and `TaxID=0` (aligned but failed QC filter) |
| `TaxID_GT0` | Reads with `TaxID > 0` (effectively classified) |
| `EffectiveClassifiedPct` | `TaxID_GT0 / AbundanceTotal` as percentage |
| `C0_ShareOfUnclassified` | `C_TaxID0 / (AbundanceTotal - TaxID_GT0)` — what fraction of unclassified reads actually aligned |

When assignments are absent, `RawC` through `C0_ShareOfUnclassified` are `NA`;
only `AbundanceTotal` and counts from `context$sample_stats` are emitted.

#### Output: `00_read_accounting_pie.png` (per sample)

A three-segment pie chart with distinct colours and a legend:

| Segment | Colour | Legend text |
|---|---|---|
| `status=U, TaxID=0` | `#bdbdbd` (grey) | "Never aligned: the read did not align/classify against the reference" |
| `status=C, TaxID=0` | `#e6ab02` (amber) | "QC-filtered: the read aligned, but was relabelled unclassified after min_ref_coverage and/or min_percent_identity filters" |
| `TaxID > 0` | `#1b9e77` (teal) | "Effectively classified: eligible for the classified abundance denominator" |

Each segment label includes the count and percentage. The centre annotation shows
the total read count.

#### Implementation

Extend `run_qc()` in `analysis/01_qc_diagnostics.R`:

- After the existing reconciliation loop, build `accounting_df` from the same
  per-sample variables (`n_status_C`, `n_status_U`, `n_qc_reclass`, `n_eff_class`).
- Write `00_read_accounting.tsv` to `qc_dir`.
- Generate the per-sample pie chart inside the existing sample loop, using
  `ggplot2::coord_polar()` as in the existing donut chart pattern but with three
  segments.
- Place the pie PNG at `{sample_out_dir}/00a_read_accounting_pie.png`.

The existing donut chart (`01a_classification_donut.png`) is retained unchanged
for backward compatibility; the new pie provides a finer decomposition.

---

### F2 — Investigative columns (lengths and threshold failures)

**Addresses:** state review §4.2 (read lengths, threshold table), user feature
request #2

#### Output: `00_read_investigation.tsv`

One row per sample. Columns:

| Column | Description |
|---|---|
| `SampleID` | Sample identifier |
| `MedianClassifiedLength` | Median `read_length` for `TaxID > 0` reads |
| `MedianC0Length` | Median `read_length` for `status=C, TaxID=0` reads |
| `MedianRawULength` | Median `read_length` for `status=U` reads |
| `IdentityOnlyFailed` | Count of C0 reads that failed only the identity threshold |
| `RefCoverageOnlyFailed` | Count of C0 reads that failed only the reference coverage threshold |
| `BothFailed` | Count of C0 reads that failed both thresholds |

##### Threshold-failure partitioning

This requires joining assignment read IDs to `bamstats.readstats.tsv.gz` from
the wf-16s output. The join is **optional** and auto-discovered from the wf-16s
output root. The post-processor will scan the output directory tree for a file
matching `*bamstats*readstats*.tsv.gz` and use the first match.

Auto-discovery is triggered when `input.wf16s_output_root` is set in the config:

```yaml
input:
  wf16s_output_root: "wf-16s_Helga16SrRNA/output"  # optional; enables bamstats auto-discovery
```

When no output root is configured or no matching file is found, the three
failure-partition columns are emitted as `NA`. When found:

1. Read the gzipped TSV.
2. Inner-join on read ID with the C0 subset from assignments.
3. Apply the configured `min_percent_identity` and `min_ref_coverage` from
   `params.json` (already available via `context$params`).
4. Partition each C0 read into identity-only-failed, coverage-only-failed, or
   both-failed.

> [!IMPORTANT]
> The `bamstats.readstats.tsv.gz` is produced by `wf-16s` but may not be present
> in the copied output directory. The post-processor must not require it.
> Auto-discovery should log a clear message when the file is found or absent.

#### Implementation

- Add `input.wf16s_output_root` to the config schema in `analysis/utils/config.R`
  (optional, default `null`). Add a `discover_bamstats()` helper that globs for
  `*bamstats*readstats*.tsv.gz` under this root.
- In `run_qc()`, after the reconciliation loop, compute median lengths from
  existing `reads` data.
- If `bamstats_readstats` is provided, read and join; otherwise emit `NA` for
  the three threshold columns.
- Write `00_read_investigation.tsv` to `qc_dir`.

---

### F3 — Species richness overview table

**Addresses:** state review §6.1 (richness driven by low-count tail), user
feature request #3

#### Output: `02_richness_overview.tsv`

One row per sample. Columns:

| Column | Description |
|---|---|
| `SampleID` | Sample identifier |
| `PositiveTaxa` | Number of classified taxa with count > 0 |
| `SingletonTaxa` | Number of taxa with count = 1 |
| `SingletonPct` | `SingletonTaxa / PositiveTaxa` as percentage |
| `TaxaLeq10` | Number of taxa with count ≤ 10 |
| `ReadsInTaxaLeq10` | Sum of reads in taxa with count ≤ 10 |
| `ReadsInTaxaLeq10Pct` | `ReadsInTaxaLeq10 / ClassifiedReads` as percentage |

This uses the classified-only count matrix (`count_matrix[-unclass_idx, ]`),
consistent with the existing unclassified handling.

#### Implementation

Place this in the alpha diversity module (`analysis/02_alpha_diversity.R`),
computed alongside the existing richness indices:

```r
class_counts_per_sample <- count_matrix[-unclass_idx, , drop = FALSE]
for (s in samples) {
  counts_s <- class_counts_per_sample[, s]
  positive   <- counts_s[counts_s > 0]
  n_positive <- length(positive)
  n_single   <- sum(positive == 1)
  n_leq10    <- sum(positive <= 10)
  reads_leq10 <- sum(positive[positive <= 10])
  # ... build row
}
```

Write to `{alpha_dir}/02_richness_overview.tsv`.

---

### F4 — Preserve current unclassified-read handling

**Addresses:** user feature request #4

No changes to the hardcoded unclassified handling in any module. Document the
current behavior explicitly in a new section of `README.md`:

| Module | Unclassified handling |
|---|---|
| `01_qc_diagnostics.R` | Included in reconciliation and all QC plots |
| `02_alpha_diversity.R` | **Excluded** — `count_matrix[-unclass_idx, ]` |
| `03_beta_diversity.R` | **Excluded** — `count_matrix[-unclass_idx, ]` |
| `04_taxa_composition.R` | **Excluded** from rank tables; **included** in classification fraction chart |
| `05_ordination.R` | **Excluded** — `count_matrix[-unclass_idx, ]` |
| `06_shared_taxa.R` | **Excluded** — `count_matrix[-unclass_idx, ]` |
| `07_kreport_pavian.R` | **Included** as the first line per kreport spec |

---

## 3. Bug fixes included

### BF1 — P0-1: Single-sample composition extraction (release blocker)

**Source:** state review §7, P0-1

Replace the dimension-dropping matrix slice in `04_taxa_composition.R` L34:

```diff
-  unclass_counts <- count_matrix[unclass_idx, ]
+  unclass_counts <- context$sample_stats$UnclassifiedReads
+  names(unclass_counts) <- context$sample_stats$SampleID
```

This uses the already-validated `context$sample_stats` as the source of truth,
eliminating the dimension-dropping scalar/name problem entirely.

Required regression assertions (added to `tests/verify_release_run.R`):

- `UnclassifiedReads == 33500` for Ambar
- `ClassifiedReads + UnclassifiedReads == TotalReads`
- No `NA` in `classification_fraction.tsv`
- Two plot-input categories exist and sum to the total
- A clean full run emits no unexpected plotting warning

### BF2 — P0-2: Reject Kraken2 abundance-only input

**Source:** state review §7, P0-2

In `build_context()`, after reading `params.json`, validate:

```r
if (!is.null(params$classifier) && params$classifier != "minimap2") {
  stop(sprintf(
    "Unsupported classifier '%s'. v0.2.0 supports minimap2 only. " %+%
    "Kraken2/Bracken requires a separate accounting model (see state review §4.3).",
    params$classifier
  ), call. = FALSE)
}
```

### BF3 — P1-1: Promote upstream contract fields into the manifest

**Source:** state review §7, P1-1

Extract and record from `params.json`:

- `classifier`, `database_set`, `taxonomic_rank`
- `min_len`, `max_len`, `min_read_qual`
- `min_percent_identity`, `min_ref_coverage`
- `abundance_threshold`
- `wf.agent` (workflow version/revision when available)

Emit these as `upstream_contract` in `run_manifest.json`.

---

## 4. File change summary

### Modified files

| File | Changes |
|---|---|
| `VERSION` | `0.1.0` → `0.2.0` |
| `analysis/utils/config.R` | Add `input.wf16s_output_root` (optional), validate classifier |
| `analysis/utils/io.R` | Add upstream contract extraction and bamstats discovery to `build_context()` |
| `analysis/01_qc_diagnostics.R` | Renumber existing plots (`00a_*` → `01a_*` already; new `00a_*`, `00b_*` for accounting); add read accounting table, pie chart, investigation table |
| `analysis/02_alpha_diversity.R` | Add richness overview table |
| `analysis/04_taxa_composition.R` | Fix P0-1 unclassified extraction |
| `analysis/00_run_pipeline.R` | Add upstream contract to manifest |
| `tests/verify_release_run.R` | Add composition value assertions, warning capture |
| `README.md` | Document unclassified handling per module |
| `CHANGELOG.md` | v0.2.0 entry |

### New files

None — all new outputs are generated by existing modules.

---

## 5. Verification plan

### Automated tests

```bash
# Unit tests
Rscript -e "testthat::test_dir('tests/testthat')"

# Full Ambar regression run
Rscript analysis/00_run_pipeline.R --config config.yml --overwrite

# Full Helga NCBI run
Rscript analysis/00_run_pipeline.R --config config_helga.yml --overwrite
```

### Manual verification

1. **Read accounting table**: cross-check `00_read_accounting.tsv` values against
   state review §4.2 table for Ambar (114,056 / 89,809 / 24,247 / 9,253 / 80,556)
   and Helga (23,099 / 23,007 / 92 / 4,852 / 18,155).
2. **Pie chart**: visually confirm three segments with correct colours, counts,
   and legend text for both Ambar and Helga.
3. **Investigation table**: verify median lengths against state review §4.2 length
   table (Helga NCBI: 1,510 / 1,471.5 / 1,468.5).
4. **Richness overview**: verify against state review §6.1 table (Ambar: 1,836
   positive, 735 singletons, 1,399 taxa ≤10, 4.29%).
5. **Composition fix**: confirm `classification_fraction.tsv` contains no `NA`,
   `UnclassifiedReads == 33500` for Ambar, and the classification fraction plot
   shows two bars summing to 100%.
6. **Kraken2 rejection**: confirm `config_july_kraken2.yml` (if created) fails
   with an explicit unsupported-classifier message.
7. **Manifest contract**: confirm `run_manifest.json` contains `upstream_contract`
   with classifier, database_set, thresholds.

---

## 6. Design decisions (resolved)

1. **bamstats join path**: auto-discover from `input.wf16s_output_root` by
   globbing for `*bamstats*readstats*.tsv.gz`. No direct path needed.
2. **Output numbering**: renumber existing QC plots to maintain consistent
   `00_` → `01_` ordering. The new read accounting outputs get `00_` prefix.
3. **Richness overview placement**: moved to `02_alpha_diversity.R` alongside
   the existing richness indices, as `02_richness_overview.tsv`.
