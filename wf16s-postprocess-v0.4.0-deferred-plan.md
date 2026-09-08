# ont-wf16S-amplicon-postprocess v0.4.0 deferred-items implementation plan

## Decision and handoff state

The repository was audited at:

- Repository: WhyAdr/ont-wf16S-amplicon-postprocess
- Current HEAD: 5870c8aba3b516153fce0ac7caee7bbb762a179d
- Current branch: main; working tree clean
- Declared version: 0.3.1
- Baseline plan: wf16s-postprocess-v0.3.1-patch.md
- Current GitHub Actions run: [run 18](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/actions/runs/34102267363)

The v0.3.1 implementation is substantially complete and the current GitHub run
passed all four required jobs. It is not yet a fully closed implementation of the
last plan, however: one valid optional-input edge case remains unhandled and two
required failure-path regressions are not actually exercised by the committed
tests. Apply the short pre-tag gate below before treating v0.3.1 as final.

After that gate, implement the four intentionally deferred items as a coordinated
v0.4.0 minor release. The manifest-schema revision and explicit not_run records
are deliberately grouped because they change the same machine-readable contract.

## Scope boundary

Preserve the v0.3.1 behavior that is already validated:

- classified-only denominators and exact read-accounting equations;
- the minimap2/species/NCBI producer contract;
- offline-by-default taxonomy resolution and run-local cache behavior;
- Kraken/Pavian arithmetic and the direct-count Krona model;
- FAPROTAX 1.2.12 via the existing microeco integration;
- current single/cohort statistical gates and output names.

Do not use this plan to redesign taxonomy resolution, add Kraken2/SILVA support,
change FAPROTAX semantics, change Krona magnitudes, or introduce a biological
validation claim for the synthetic cohort fixture.

## Audit result against the v0.3.1 plan

| Plan area | Current implementation | Audit result |
|---|---|---|
| Patch 1: early module/dependency validation | parse_requested_modules() exists; duplicate/unknown modules and the Krona dependency are checked before build_context() | Implemented; one required test is miswired, described below |
| Patch 2: producer, lineage, metadata, and identifier contracts | Whole-number thresholds, rank-cell validation, character metadata identity columns, portable SampleID checks, YAML-root validation, and padded discovery names are present | Implemented; explicit 1.0 acceptance and trailing-whitespace test coverage should be completed |
| Patch 3: safe sample seeds | Base seed is range-checked and derived seeds use double-precision modular arithmetic | Implemented and covered |
| Patch 4: graphics and one-row heatmaps | Shared PNG wrapper is used by modules 04 and 06; row/column clustering is gated | Implemented and covered |
| Patch 5: NMDS and beta-resampling failures | NMDS repetition count is retained; rarefaction, PCoA, and Procrustes are inside one iteration guard | Code implemented; injected one-failure and all-failure tests are missing |
| Patch 6: release metadata and verifier | Active metadata is 0.3.1; verifier reads VERSION; no active stale 0.3.0 hit remains in VERSION, CITATION.cff, README.md, analysis/, or tests/ | Implemented |

### Verification evidence

Observed locally:

- python -m unittest -v tests/test_ncbi_taxonomy.py: 6/6 passed;
- python -m compileall -q analysis tests: passed;
- python tests/check_committed_whitespace.py: passed;
- git diff --check: passed;
- R was not installed in the audit environment, so local Rscript/testthat and
  end-to-end R runs were not independently rerun.

Observed remotely on the current HEAD:

- Ubuntu / R 4.5 / Python 3.12: success;
- Windows / R 4.5 / Python 3.12: success;
- Ubuntu / R 4.5 / Python 3.12 (FAPROTAX): success;
- Windows / R 4.5 / Python 3.12 (FAPROTAX): success.

No v0.3.1 tag is present in the checked-out repository, so the tag itself has
not yet been verified.

## Pre-tag gate for v0.3.1

### P1 — Preserve numeric-looking bamstats identifiers

Finding: analysis/utils/io.R::partition_minimap2_failures() reads the full
bamstats table with an unconstrained read.delim(). R can type-convert values in
the name and sample_name columns. Consequently, a valid assignment/read pair
such as read ID 0001 and sample ID 01 can become 1 during the bamstats read and
fail matching, even though discovery correctly reads sample_name as character.

This edge is not covered by the current tests. It is a correctness defect in the
optional bamstats path, not merely a cosmetic type issue.

Files:

- analysis/utils/io.R
- tests/testthat/test-io.R

Required implementation:

1. Add one internal read_bamstats_table() helper and use it in both
   discover_bamstats() and partition_minimap2_failures().
2. Read the header first and reject duplicate headers or missing required columns:
   name, sample_name, iden, and ref_coverage.
3. Import all four required columns as character. Convert only iden and
   ref_coverage explicitly with as.numeric() at the point where their finite
   numeric contract is checked.
4. Keep name and sample_name byte-for-byte as supplied; do not trim them for
   matching. Continue rejecting missing, empty, or boundary-whitespace
   sample_name values.
5. Preserve the existing one-to-one read-name and threshold-partition semantics.

The helper may read only the required columns by assigning NULL to unrelated
columns in colClasses, but it must preserve the required columns by name and must
not rely on their physical order.

Regression: create a gzipped bamstats fixture with sample_name = "01" and
name = "0001", an assignment row with read_id = "0001", and thresholds that
classify it as identity_only. Require discovery to map the file to "01" and
partitioning to return one matched read rather than a sample-name or missing-read
error. Include a second ordinary alphanumeric case to prove the refactor preserves
the current behavior.

### P2 — Make the malformed-assignment fail-fast test real

tests/testthat/test-release-process.R writes a file named assignments.tsv in the
fast-fail test but does not place that file in cfg$input$assignments. The test
therefore proves only that an unknown module fails; it does not prove that the
unknown-module diagnostic wins over a malformed assignment.

Pass the malformed path into write_process_config() as an assignment mapping,
then require:

- non-zero process status;
- Unknown module in stderr;
- no assignment-schema error in stderr;
- no output directory created.

### P3 — Add the required beta-resampling failure regression

The current run_beta() guard is correctly broad, but no committed test injects a
Procrustes failure. Add a narrow test seam rather than relying on fragile namespace
monkey-patching:

1. Extract the per-iteration alignment operation into a private helper or allow
   run_beta() to receive an internal defaulted procrustes_fn = vegan::procrustes
   argument. The normal runner must continue to use vegan unchanged.
2. Use a six-sample synthetic cohort, enable beta resampling, and request three
   iterations.
3. Inject a function that fails only on iteration 1 and delegates to vegan for the
   remaining iterations.
4. Require module status completed, SuccessfulIterations == 2,
   FailedIterations == "1", and stability rows for the two successful iterations.
5. Add an all-iterations-fail case and require the existing fail-closed error
   All beta-diversity rarefaction stability iterations failed.

### P4 — Complete small contract-test gaps

Add the following low-cost assertions while the v0.3.1 patch is still open:

- explicitly set a producer threshold to numeric 1.0 and require acceptance;
- test both leading and trailing SampleID whitespace;
- retain the existing tests for fractional min_len, max_len, and
  abundance_threshold.

### Pre-tag acceptance

After P1–P4:

1. Run the full test and verification matrix on a fresh commit.
2. Repeat the active-version stale-hit check.
3. Run both release integrations with fresh output roots.
4. Confirm all four CI jobs are green on the exact tag-target SHA.
5. Only then create the v0.3.1 tag.

## v0.4.0 Phase 1 — Versioned manifest schema with stable arrays

### Problem being deferred

The runner and Krona provenance currently serialize with auto_unbox = TRUE.
Fields such as samples, cli.modules, module outputs, and package/version
collections can therefore be a scalar for one item and an array for multiple items.
This is a real consumer-facing schema instability. Do not silently change schema v1
in a patch release.

### Files

- new analysis/utils/manifest.R;
- analysis/00_run_pipeline.R;
- analysis/07_kreport_pavian.R;
- tests/verify_release_run.R;
- new tests/testthat/test-manifest-schema.R;
- README.md;
- CHANGELOG.md.

### Contract

1. Introduce a distinct manifest schema version. For v0.4.0, set the manifest's
   schema_version to 2 and add config_schema_version: 1; the resolved YAML
   configuration continues to use its own schema_version: 1.
2. The following fields must always be JSON arrays, including when they contain
   zero or one item:

   - manifest.samples;
   - manifest.command;
   - manifest.cli.modules;
   - manifest.warnings;
   - manifest.package_versions;
   - manifest.inputs.assignments and manifest.inputs.bamstats when present;
   - every modules.<name>.outputs and modules.<name>.warnings;
   - krona_provenance.samples.

3. Optional absent objects remain explicit JSON null; they must not become {}.
4. Scalar strings, numbers, and booleans remain JSON scalars.
5. Empty collections are [], not null, omitted keys, or an empty string.
6. Preserve field names and the existing analytical meanings unless a field is
   explicitly versioned below.

### Serialization implementation

Add small constructors in manifest.R, for example:

~~~r
json_array <- function(x) {
  unname(lapply(as.list(x), jsonlite::unbox))
}
~~~

Build the manifest with explicit array containers and serialize with
auto_unbox = TRUE. Do not depend on the cardinality of an atomic vector to imply
an array. Apply the same construction to Krona provenance. Add a single writer that
validates the in-memory manifest before writing it, so future fields cannot silently
reintroduce singleton collapse.

The validator should check the v2 required keys, scalar/array/null shape, and that
module records contain valid statuses and arrays. It should fail with the field path
in the error message.

### Reader and verifier changes

Update tests/verify_release_run.R to require manifest schema 2 and
config_schema_version == 1. Replace implicit unlist() compatibility with an
explicit array reader that rejects a scalar where an array is required. Keep a
small compatibility reader/test for archived v0.3.1 manifests if external users
may still inspect them; do not rewrite old output artifacts in place.

### Tests

Generate one-sample and multi-sample manifests, parse them with
jsonlite::fromJSON(..., simplifyVector = FALSE), and assert identical JSON shapes
for every cardinality-dependent field. Include empty warnings and empty output
arrays. Test one-sample and multi-sample Krona provenance separately. Add a negative
fixture with a scalar samples field and require a schema error naming the field.

### Documentation

Document the v2 manifest contract and the v1-to-v2 boundary in the README. State
that v0.4.0 changes machine-readable manifest shape but does not change analytical
denominators or output semantics.

## v0.4.0 Phase 2 — Streaming assignment parsing

### Problem being deferred

read_assignments_file() currently reads every physical line, creates a second
collapsed text buffer with paste(), and then performs a typed import. This is
acceptable for the tracked fixture but creates avoidable peak memory pressure for
large PromethION assignment files.

### Files

- analysis/utils/io.R, or a new analysis/utils/assignments.R sourced by the
  runner and tests;
- tests/testthat/test-io.R;
- new tests/benchmark_assignment_parser.R;
- README.md.

### Required parser design

Replace the whole-file line vector plus paste() path with a chunked connection
reader. A two-pass implementation is preferred:

1. Schema/count pass: open plain or gzipped input through a connection, read
   bounded chunks with readLines(con, n = chunk_size), maintain physical line
   numbers, skip blank lines exactly as today, validate five fields, status, TaxID,
   length, and duplicate read IDs, and accumulate only counters/validation state.
2. Typed import pass: reopen the connection, count rows already known to be
   valid, preallocate the typed output columns, and fill them chunk by chunk.
3. Reconcile total/effective-classified/effective-unclassified counts against the
   abundance expectations and return the same columns and types as v0.3.1.

Keep these invariants byte-for-byte or message-for-message where currently tested:

- blank-line handling;
- physical line numbers in schema errors;
- exact five-field validation before type normalization;
- read_id as character and unique;
- taxid as a non-negative integer;
- pipe-delimited and plain length fields;
- effective_classified = taxid > 0;
- all existing reconciliation errors.

Do not stream away the final typed data frame: QC still needs the read-level vectors.
The objective is to remove the duplicate raw-text and collapsed-buffer copies, not to
change downstream module interfaces.

### Tests and benchmark

- Compare the streamed parser with the v0.3.1 reference behavior on the tracked
  assignment fixture: row count, column names, types, classification counts,
  length summaries, and a deterministic digest of the returned typed columns.
- Retain malformed-field, duplicate-ID, blank-line, plain-text, and gzipped tests.
- Add a bounded synthetic large-file smoke test that exercises chunk boundaries.
- tests/benchmark_assignment_parser.R must report row count, elapsed time, final
  object size, and peak-memory evidence using the host-appropriate mechanism. Run
  the benchmark manually on Linux and Windows; do not make a large benchmark file
  part of the repository.
- Record the tracked-fixture baseline and the large-file scaling result in the
  v0.4.0 release notes. The benchmark is a performance gate, not a license to
  weaken schema validation.

## v0.4.0 Phase 3 — Explicit not_run module records

### Problem being deferred

When a module fails and --keep-going is not enabled, later requested modules are
absent from the manifest. Consumers cannot distinguish not requested from was
requested but never executed because the pipeline stopped.

### Files

- analysis/00_run_pipeline.R;
- analysis/utils/manifest.R;
- tests/testthat/test-release-process.R;
- tests/verify_release_run.R;
- README.md.

### Required behavior

1. After request validation, initialize one manifest record for every requested
   module in request order with status not_run, empty outputs, empty warnings,
   and null execution timestamps/duration.
2. Replace that record with the real result when the module executes. Valid runtime
   statuses are completed, skipped, failed, and not_run.
3. When a module fails with --keep-going absent, retain later requested modules as
   not_run with a reason such as:
   Pipeline stopped after failure in module 'kreport'.
4. Do not invoke a not_run module and do not create its output directory merely
   because its manifest record exists.
5. Preserve current --keep-going execution semantics. A failed module remains
   failed; later modules execute normally.
6. A run with any failed module remains top-level run_status = failed and exits
   non-zero after writing the manifest.

Define the record contract explicitly:

| Field | Completed/skipped/failed | Not run |
|---|---|---|
| status | one of the three executed statuses | not_run |
| outputs | JSON array | empty JSON array |
| warnings | JSON array | empty JSON array |
| error | null or error string | null |
| reason | null or skip reason | stop-after-failure reason |
| start_time, end_time, duration_seconds | recorded as applicable | JSON null |

### Tests and verifier

Use a process-level fixture where kreport fails under unresolved_policy: error
and composition and alpha are requested after it. Require all three module records
to exist, with kreport = failed and the two later modules = not_run. Require that
their output directories do not exist. Test the existing keep-going case separately
and require subsequent modules to remain executable/completed.

The v0.4.0 release verifier must reject not_run in a successful release integration
run, while general manifest validation must accept it for failed fail-fast runs.

## v0.4.0 Phase 4 — Lock the R analysis environment

### Problem being deferred

The project records package versions and sessionInfo() but has no committed
renv.lock. Installing latest CRAN packages can change numerical behavior, plotting,
or the microeco-embedded FAPROTAX database over time.

### Files

- renv.lock;
- standard renv bootstrap files (.Rprofile, renv/activate.R, and generated
  renv/settings.json if produced);
- .gitignore for renv cache/library directories only;
- .github/workflows/ci.yml;
- analysis/install_packages.R;
- analysis/00_run_pipeline.R;
- tests/testthat/test-release-metadata.R or a new environment test;
- README.md;
- CHANGELOG.md.

### Required lock scope

1. Generate the lock from the final green Ubuntu R 4.5 reference environment and
   commit the complete transitive R dependency graph, including testthat and
   microeco.
2. Pin the exact R patch version used to create the lock in CI; do not leave CI at a
   floating 4.5 if the lock is intended to be reconstructive.
3. Keep Python at the documented 3.12 line. The taxonomy resolver currently uses
   Python standard-library modules only; do not invent a Python dependency lock
   unless a non-standard dependency is added.
4. Record the renv.lock SHA-256 and locked-environment status in the v2 run
   manifest. The lock covers this post-processing environment, not the upstream
   wf-16s workflow, ONT basecaller, NCBI database archive, or optional KronaTools
   executable.
5. Keep the existing runtime check that the installed microeco database is exactly
   FAPROTAX 1.2.12. A lock is complementary to that semantic check.

### CI and user workflow

- Restore from renv.lock in both Ubuntu and Windows core jobs and both FAPROTAX
  jobs before running tests.
- Add a CI assertion that renv::status() reports no drift after restore.
- Avoid replacing the lock with an unconstrained any:: dependency installation.
- Keep analysis/install_packages.R useful for discovery, but add an explicit
  locked/restore path and document renv::restore(prompt = FALSE) as the release
  and publication setup.
- Verify that the lock restores on both supported operating systems and that the
  FAPROTAX jobs still select the locked compatible microeco release.

### Reproducibility tests

- Require the lock file to exist and parse as valid JSON.
- Require the manifest lock hash to equal the repository lock hash.
- Run the tracked single-sample integration twice from fresh output roots and
  compare analytical TSV digests where deterministic; allow timestamps, absolute
  paths, and environment text to differ.
- Document which outputs are provenance-reproducible versus bitwise-stable.

## Optional v0.4.0 release-hygiene item

The repository still tracks output_Helga/run_manifest.json, whose historical
manifest advertises pipeline version 0.1.0 and contains machine-specific absolute
paths. This is outside the active v0.3.1 stale-version scan and is not a current CI
failure, but it can confuse users inspecting repository artifacts.

Choose one explicit treatment in a separate hygiene commit:

- regenerate the complete tracked Helga artifact set with the current release and
  sanitize or document paths; or
- declare it a historical artifact in the README and exclude it from release
  verification; or
- remove the generated artifact set only after confirming it is not a required
  fixture or user reference.

Do not silently rewrite or delete it as part of the schema/performance work.

## v0.4.0 implementation order

1. Complete the v0.3.1 pre-tag gate and tag v0.3.1.
2. Add the manifest v2 model and shape validator.
3. Add explicit not_run records on top of the v2 model.
4. Replace the assignment parser and run correctness/performance benchmarks.
5. Generate and verify renv.lock; update CI to restore it.
6. Update README and changelog only after behavior and CI are green.
7. Run the full four-job matrix from a clean checkout and validate both successful
   and intentionally failed manifest examples.
8. Tag v0.4.0 only after all acceptance items below are true.

## v0.4.0 acceptance checklist

- [ ] The numeric-looking bamstats regression passes before v0.3.1 tagging.
- [ ] The malformed-assignment/unknown-module test actually configures the malformed
      file and proves fail-fast ordering.
- [ ] One injected beta-resampling failure is recorded without aborting successful
      iterations; all failed iterations still fail closed.
- [ ] v2 manifests encode all cardinality-dependent fields as arrays for zero,
      one, and many items.
- [ ] v2 manifest validation rejects scalar/array shape drift with a field path.
- [ ] Successful and skipped module records have stable fields and arrays.
- [ ] Fail-fast runs include later requested modules as not_run; their output
      directories are not created.
- [ ] --keep-going still executes later modules and preserves a failed top-level
      status.
- [ ] Streaming assignment parsing preserves every v0.3.1 schema, type, count,
      and physical-line error invariant.
- [ ] The large-file benchmark demonstrates removal of the duplicate raw-text
      buffering path and records Linux/Windows memory evidence.
- [ ] renv.lock restores on Ubuntu and Windows with the exact supported R patch
      version and the compatible FAPROTAX 1.2.12 microeco data.
- [ ] Run manifests record and validate the lockfile hash.
- [ ] No denominator, taxonomy, FAPROTAX, or Krona-count semantics changed.
- [ ] All four CI jobs are green on the final v0.4.0 SHA.
