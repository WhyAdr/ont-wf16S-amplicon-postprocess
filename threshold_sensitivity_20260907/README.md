# BAER, BANAE, BLEA, and BGRN threshold sensitivity

This is a read-level sensitivity reconstruction from the saved native wf-16s
outputs. It does not alter or replace those outputs and does not rerun minimap2.
The coverage threshold is **reference coverage**, matching `min_ref_coverage`.

## Baseline 90/90 read accounting

| Dataset | Total | Classified | Unclassified | Raw unmapped | Mapped TaxID 0 | Identity only | Ref. coverage only | Both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BAER | 74,442 | 49,206 (66.10%) | 25,236 | 61 | 25,175 | 11,774 | 8,650 | 4,751 |
| BANAE | 52,474 | 30,492 (58.11%) | 21,982 | 74 | 21,908 | 10,940 | 6,250 | 4,718 |
| BLEA | 60,694 | 28,996 (47.77%) | 31,698 | 50 | 31,648 | 14,742 | 9,881 | 7,025 |
| BGRN | 124,483 | 107,605 (86.44%) | 16,878 | 78 | 16,800 | 5,052 | 10,231 | 1,517 |

Every baseline-classified read passed both recorded thresholds. Every mapped
TaxID-0 read failed identity, reference coverage, or both. All existing
classified TaxIDs were reproduced exactly by joining the bamstats best reference
to the alignment table.

## Main cause

- **BAER**: 99.76% of unclassified reads were mapped but threshold-rejected; identity was the more frequent implicated threshold.
- **BANAE**: 99.66% of unclassified reads were mapped but threshold-rejected; identity was the more frequent implicated threshold.
- **BLEA**: 99.84% of unclassified reads were mapped but threshold-rejected; identity was the more frequent implicated threshold.
- **BGRN**: 99.54% of unclassified reads were mapped but threshold-rejected; reference coverage was the more frequent implicated threshold.

Read length alone does not explain the rejected reads; median classified and
mapped-TaxID-0 lengths are reported in `baseline_failure_breakdown.tsv`.

## Selected simulations

The scenario name and both numeric columns make the cutoff order explicit:
identity first, reference coverage second.

| Dataset | Scenario | Min identity | Min ref. coverage | Newly rescued | Classified total |
|---|---|---:|---:|---:|---:|
| BAER | identity_85 | 85 | 90 | 10,286 | 79.92% |
| BAER | refcov_85 | 90 | 85 | 4,329 | 71.92% |
| BAER | balanced_89 | 89 | 89 | 3,562 | 70.88% |
| BAER | balanced_88 | 88 | 88 | 6,709 | 75.11% |
| BAER | balanced_85 | 85 | 85 | 16,452 | 88.20% |
| BANAE | identity_85 | 85 | 90 | 9,563 | 76.33% |
| BANAE | refcov_85 | 90 | 85 | 2,800 | 63.44% |
| BANAE | balanced_89 | 89 | 89 | 2,609 | 63.08% |
| BANAE | balanced_88 | 88 | 88 | 5,742 | 69.05% |
| BANAE | balanced_85 | 85 | 85 | 14,105 | 84.99% |
| BLEA | identity_85 | 85 | 90 | 14,267 | 71.28% |
| BLEA | refcov_85 | 90 | 85 | 4,198 | 54.69% |
| BLEA | balanced_89 | 89 | 89 | 5,032 | 56.06% |
| BLEA | balanced_88 | 88 | 88 | 9,678 | 63.72% |
| BLEA | balanced_85 | 85 | 85 | 21,300 | 82.87% |
| BGRN | identity_85 | 85 | 90 | 4,748 | 90.26% |
| BGRN | refcov_85 | 90 | 85 | 5,433 | 90.81% |
| BGRN | balanced_89 | 89 | 89 | 2,833 | 88.72% |
| BGRN | balanced_88 | 88 | 88 | 5,402 | 90.78% |
| BGRN | balanced_85 | 85 | 85 | 10,772 | 95.09% |

## Outputs

- `baseline_failure_breakdown.tsv`: baseline cause counts and QC medians.
- `threshold_grid_80_to_90.tsv`: complete 121-scenario grid per dataset.
- `selected_threshold_scenarios.tsv`: compact set of one-axis and balanced scenarios.
- `unclassified_read_diagnostics.tsv.gz`: every baseline-unclassified read with
  alignment metrics, failure mode, and inferred best-reference taxonomy.
- `selected_scenario_rescued_taxa.tsv`: inferred taxa introduced by selected scenarios.
- `provenance.json`: parameters, input paths, and SHA-256 identities.

The inferred TaxIDs are suitable for cutoff selection and sensitivity analysis.
Use a native wf-16s rerun for canonical workflow reports and final abundance
artifacts after selecting a threshold.
