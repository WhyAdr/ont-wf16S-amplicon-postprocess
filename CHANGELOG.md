# Changelog

All notable changes to this project are documented in this file.

## [0.4.3] - 2026-09-08

### Fixed

- **Manifest Schema v2 Revision 2**: Bump schema revision to `2L` with mandatory `artifacts` array (tracking size, sha256, producer module) and `preserved_unowned_outputs` array.
- **Physical Output Census**: Enforce exact equality between physical staged files and declared owned plus preserved files before publication (`physical == (owned \ manifest) ∪ preserved`).
- **Transactional Replace & Crash Recovery**: Replace unsafe rename-over-deleted semantics with atomic temporary replacements and sibling transaction journals with deterministic crash recovery.
- **Output Concurrency Locking**: Acquire cross-platform output-root locks via `filelock` to eliminate races during parallel runs against identical output directories.
- **Taxonomy Cache State Machine & Concurrency**: Implement explicit cache state machine (`unchanged`, `candidate_committed`, `restored`) with tempfile backup, SHA-256 pre-verification, cross-platform file locking (`fcntl` / `msvcrt`), and fail-closed race abortion (`E_TAXONOMY_CACHE_CHANGED`).
- **Deferred Taxonomy Commit & Recovery**: Keep refresh candidates run-local until the kreport module succeeds, journal the source-cache commit, restore it when a later module or publication fails, and forward-recover only from a strictly validated completed output.
- **Conflict Provenance**: Retain the compatibility modal/minimum TaxID tie-break for assignment conflicts while explicitly labelling affected rows `Conflicted` with `assignment_conflict` provenance and a run warning.
- **Private Per-Module Staging**: Replace whole-stage directory copying with isolated per-module staging roots, eliminating quadratic disk I/O and enforcing module output boundaries.
- **Environment & Preflight Hardening**: Make external-CWD loading self-sufficient via `env_loader.R`, enforce project-library containment for all non-base packages, include startup files (`.Rprofile`, `renv/activate.R`, `renv/settings.json`) in source provenance, fingerprint the supplied configuration, and record available lock/install identity metadata for each package.
- **Environment Activation Idempotence**: Avoid redundant `renv::load()` calls when the canonical project library is already active, keeping subprocess validation deterministic on Windows while retaining bootstrap for external environments.
- **Failure Retry & Recovery Validation**: Permit explicit overwrite of a failed current-contract run, reject unowned staging collisions before mutation, and require a complete schema-v2 manifest, artifact hashes, and physical census before trusting a crash-recovery candidate.
- **Release CI Assertions**: Map each GitHub event to an explicit diff base and unit-test root, multi-commit push, new-branch, and pull-request range selection.
- **Cross-Platform Release Restore**: Install the Ubuntu GLPK runtime required by the locked `igraph` binary and force the byte-strict `VERSION` file to LF on Windows checkouts.
- **Exact Release Identity**: Require release-verification outputs to come from the current exact Git commit, a clean source tree, a synchronized lock, and no development escape hatches.
- **AST Syntax Checking**: Replace bytecode compilation in preflight with read-only AST syntax parsing.

## [0.4.2] - 2026-09-07

### Fixed

- Make release verification enforce the complete eight-module registry and
  explicit `not_run` records for unrequested modules.
- Register runner cleanup inside a transactional function and roll back module
  writes on failure, including safe migration of legacy output ownership.
- Activate and record the project library with full lockfile dependency closure,
  stabilize source digests, and enforce immutable resolver inputs.
- Align assignment parsing and NCBI refresh context checks across runtimes.
- Require a minimum beta-resampling success fraction and propagate statistical
  warnings into module and run manifests.
- Correct push-range whitespace validation for main-branch CI runs.

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
