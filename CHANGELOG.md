# Changelog

All notable changes to this project are documented in this file.

## [0.4.1] - 2026-09-07

### Fixed

- Make lock, source, and input provenance truthful; reject dirty or unsynchronized
  publication runs unless an explicit development escape hatch is supplied.
- Add side-effect-free module preflight, immutable-input checks, strict output-root
  ownership, and transactional staged publication for overwrite and failure paths.
- Enforce manifest v2 array/state invariants, exact integer/resource bounds,
  identifier and assignment contracts, pinned upstream database resources, and
  canonical string TaxIDs.
- Record one primary beta-diversity distance and shared constrained permutations;
  skip degenerate rarefaction stability with diagnostics.
- Make taxonomy refresh candidate-based and preserve the source cache on rejected
  or failed resolutions; accept uppercase gzip assignment paths consistently.

### Changed

- Add `--allow-unlocked`, `--allow-dirty`, and `--online-preflight` development
  controls; normal publication runs require synchronized provenance.
- Remove stale tracked runtime output from version control and strengthen CI and
  whitespace checks across the complete commit range.

## [0.4.0] - 2026-09-07

### Added

- Manifest schema v2 with stable JSON arrays for cardinality-dependent fields,
  explicit `not_run` records after fail-fast module errors, and locked-environment
  SHA-256 provenance.
- A committed R 4.5.3 `renv.lock`, restore workflow, and CI restoration checks
  for the full transitive analysis/test/FAPROTAX package graph.
- A manual streamed-assignment benchmark harness. On the release Windows host,
  250,000 synthetic rows completed in 13.33 seconds with a 30.53 MiB final data
  frame and 343.27 MiB peak working set.

### Changed

- Replace whole-file assignment text buffering with bounded two-pass plain or
  gzip parsing while retaining typed output, duplicate-ID checks, reconciliation,
  and physical-line schema diagnostics.
- Remove the unused `ggrepel` runtime dependency from the locked environment.

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

## [0.3.0] - 2026-09-07

### Added

- Opt-in Krona-compatible taxonomy export through `--krona` or the `krona`
  configuration block. R always writes portable `.krona.tsv` files from the
  validated kreport tree; optional `.krona.html` rendering uses KronaTools
  `ktImportText` when available.
- Krona provenance with exact read accounting, direct-count semantics, and
  explicit renderer availability/status.

## [0.2.0] - 2026-09-07

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
