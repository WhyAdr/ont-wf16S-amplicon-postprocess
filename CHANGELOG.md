# Changelog

All notable changes to this project are documented in this file.

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
