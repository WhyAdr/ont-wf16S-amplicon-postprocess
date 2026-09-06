# Changelog

All notable changes to this project are documented in this file.

## [0.2.0] - Unreleased

### Added

- Opt-in FAPROTAX 1.2.12 taxon-based functional inference through microeco,
  with explicit non-exclusive function counts, read-accounting coverage, and
  provenance outputs (`--modules faprotax`).

- Exact per-sample read accounting with a three-way diagnostic donut.
- Optional one-to-one bamstats joins for identity-only, reference-coverage-only,
  and dual-threshold minimap2 failures.
- Classified-only richness summaries for singleton and low-count tails.
- Structured upstream producer-contract and bamstats provenance.

### Fixed

- Preserve the single-sample unclassified count in composition outputs instead
  of producing `NA` and a misleading 100%-classified plot.
- Reject unsupported Kraken2, SILVA, non-species, and custom-reference contracts
  before classifier-specific assignment parsing.

### Documentation

- State each module's treatment of unclassified reads and the evidentiary limit
  of low-count richness summaries.

## [0.1.0] - 2026-09-05

### Added

- Strict, config-driven parsing and reconciliation for tracked `wf-16s`
  abundance and minimap2 assignment contracts.
- Single-sample QC, alpha diversity, taxonomic composition, and standard
  six-column Kraken report modules.
- Cohort beta diversity, ordination, and shared-taxa modules with explicit
  sample-size gates and synthetic regression coverage.
- Offline-by-default taxonomy resolution, run-local cache enrichment,
  unresolved/conflict diagnostics, and per-node resolution-source provenance.
- Machine-readable run manifests with semantic version, Git revision, command,
  interpreter/package versions, input hashes, module outcomes, and taxonomy
  provenance.
- Ubuntu and Windows CI with unit, regression, mutation-free validation, and
  full tracked-fixture integration gates.

### Known limitations

- Cohort behavior is tested with synthetic data and is not biological
  validation of cohort statistics or species-level classification accuracy.
- The tracked offline reference run intentionally retains 46 unresolved
  taxonomy nodes and 26 lineage-to-TaxID conflicts.
- This release records environment provenance but does not yet provide a
  committed `renv.lock`; it is not bitwise environment-reproducible.
- Pavian HTML is not generated automatically; `.kreport` files are provided for
  upload to Pavian.
