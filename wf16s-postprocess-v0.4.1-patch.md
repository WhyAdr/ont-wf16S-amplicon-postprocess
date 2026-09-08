# `wf16s-postprocess` v0.4.1 hardening and fail-closing patch

**Implementation target:** `b6c028c02c364c00e64e2ad6d337c34deda48530` (`VERSION` 0.4.0, `origin/main`)
**Audit date:** 2026-09-07
**Release decision:** **No-go for a v0.4.1 tag until the P1 gates and their regression tests below are complete.**

## 1. Executive summary

v0.4.0 does implement most of the deferred v0.3.1 work: the previously missing regression cases are present, module names are checked early, assignment input is read in bounded chunks rather than buffered as one text object, manifest schema v2 records `not_run`, and an R 4.5.3 `renv.lock` is committed and restored in CI.

The new audit nevertheless found four release-blocking contract failures:

1. A manifest calls the environment `locked` merely because `renv.lock` exists; it does not establish that the running R and package library match that lock.
2. Inputs are read, re-read, and hashed at different times without an immutability check. A run can therefore combine multiple revisions of one input while reporting only one hash.
3. `--overwrite` is not a clean or transactional rerun. Old products can survive a subset rerun and failed modules can leave partial products that the manifest omits.
4. `--validate-only` is not predictive for requested modules. In particular, it can pass before kreport fails for Python/cache/resolver reasons or QC fails during full bamstats reconciliation.

The audit also found strictness defects in manifest v2, integer validation, identifier namespaces, beta-diversity design/provenance, taxonomy refresh, and upstream database verification. These are suitable for v0.4.1 because each closes an invalid-input or misleading-provenance path without intentionally changing the analytical denominators.

Do **not** change the following scientific semantics in this patch:

- abundance denominators and unclassified-read handling;
- lineage rank order and existing valid taxonomy normalization;
- FAPROTAX 1.2.12 non-exclusive function-count semantics;
- six-column Kraken report direct/clade count semantics;
- Krona TSV count semantics;
- existing valid-input results except where beta permutation design was previously inconsistent or provenance was incomplete.

## 2. Audit scope and evidence

Reviewed end-to-end: the maintained runner, all eight modules, all maintained utilities, Python resolver, configuration, test suite, CI workflow, `renv.lock`, release verifier, README/changelog/citation metadata, and tracked fixtures. Root-level legacy scripts were checked for isolation from the maintained entry point, but remain retrospective-only as documented and are not release runtime.

Local checks completed:

| Check | Result |
|---|---|
| Python resolver unit suite | 6/6 passed |
| Python compilation (`analysis`, `tests`) | Passed |
| JSON/YAML parse sweep | Passed |
| gzip fixture integrity | Passed |
| committed-whitespace script | Passed under its current one-commit scope |
| `git diff --check` / repository integrity | Passed |
| Exact audited revision | `b6c028c02c364c00e64e2ad6d337c34deda48530` |
| Tag pointing at audited revision | None found after fetching tags |
| R/testthat/integration execution on this host | Not run: R is unavailable in the audit environment |
| Current remote GitHub Actions result | Not independently verified from this environment |

Two direct adversarial probes reproduced additional defects:

- a valid gzip assignment named `assignments.tsv.GZ` reaches R's case-insensitive gzip path but the Python resolver opens it as plain text and raises `UnicodeDecodeError`;
- in taxonomy `refresh` mode, a mocked partial NCBI resolution followed by remaining unresolved nodes returned exit 1 under `unresolved_policy=error`, while the source cache had already changed from `{}` to `{"Bacteria": 123}`.

## 3. Validation of the v0.4.0 deferred items

| Deferred item | Status at v0.4.0 | v0.4.1 consequence |
|---|---|---|
| Close v0.3.1 missing tests | Implemented | Preserve them; add the adversarial matrix below. |
| Manifest stable arrays | **Partial** | Most fields use arrays, but zero assignment/bamstats collections serialize as `null`, contrary to the README promise. |
| Explicit `not_run` state | Implemented | Complete the state-machine invariants and result-envelope checks. |
| Streamed assignment parser | Implemented, with caveat | Text buffering is bounded, but two passes have no same-file fingerprint gate and the final typed data frame is still O(rows). Document this accurately. |
| Committed reproducible environment | **Partial** | Lock exists and CI restores it; runtime provenance does not prove synchronization, and the installer can target the wrong project. |
| v0.4.0 release metadata | Implemented | No release tag was visible; do not tag v0.4.1 until all gates pass. |

## 4. Findings

### P1 — release blockers

#### P1-01 — `environment.locked` is a false assertion

`analysis/00_run_pipeline.R:318-365` sets `environment.locked` from `file.exists(renv.lock)`. A global/latest package library or `analysis/install_packages.R --install` run is therefore reported as locked even if it differs from the lock. The README strengthens the false claim by saying the recorded lock was “used for the run.”

**Required fix**

- Before reading large inputs or creating output directories, compare the active R version and the required package versions against the repository lock.
- Make synchronized execution the default. Add one explicit escape hatch, `--allow-unlocked`, for development only.
- Record `lock_status` as `synchronized`, `mismatch`, `missing`, or `not_checked`; record expected/actual R and package discrepancies as stable arrays.
- Set `locked: true` only for `synchronized`. If the opt-out is used, set it false, emit a prominent warning, and preserve the discrepancy records.
- For a normal publication run, fail before output mutation on missing or mismatched lock state.

#### P1-02 — input time-of-check/time-of-use can corrupt provenance

The context reads inputs, modules can read them again, and final manifest metadata/hashes are assembled later. The assignment parser itself performs two passes without proving that pass two saw the same bytes as pass one. The Python resolver independently re-reads abundance, assignments, and cache.

**Required fix**

- Resolve and inventory every input before parsing: canonical path, file identity where available, byte size, UTC mtime, and SHA-256.
- Prefer a run-local read-only snapshot for all normal inputs. Otherwise recheck the full fingerprint before and after every multi-pass/re-reader boundary and again before manifest publication.
- Abort on any change; never publish a successful manifest for a mixed-revision run.
- Treat taxonomy refresh as the only intentional mutable-source operation. Record before/candidate/committed hashes and update the source only after every acceptance gate succeeds.
- Reject duplicate canonical assignment paths even when configured under different samples.

#### P1-03 — output replacement is neither clean nor transactional

The README already concedes that `--overwrite` replaces same-named files but does not remove stale artifacts. A subset rerun into an old output root can therefore present old, unrequested products beside new products. A failed module can also leave partial files, while its caught result records `outputs = []`. The run manifest itself is written non-atomically.

**Required fix**

- Validate the output root before mutation: it must not be `/`, a drive root, a home directory, the repository root, an input path/ancestor, or nested within the upstream input root.
- A non-empty output directory is reusable only if it contains a valid prior manifest proving pipeline ownership. Otherwise reject it even with `--overwrite`.
- Stage the whole run, or at minimum each module, under a same-filesystem temporary directory. Publish only a completed/skipped module as one rename operation; discard a failed module's staging area.
- On `--overwrite`, remove/replace the complete set of pipeline-owned products from the prior manifest, including outputs from modules not requested this time. Never delete unowned paths.
- Write manifest, resolved config, session information, TSVs, JSON, and provenance with temp-file-plus-rename semantics. The manifest is the final commit marker.
- A failed or interrupted rerun must leave either the previous complete run or an explicitly failed new run, never a hybrid.

#### P1-04 — `--validate-only` is not predictive

The current preflight does not fully validate the kreport Python executable/script/cache/resolution path and does not run the same bamstats reconciliation used by QC. Thus validation can succeed before an unchanged real invocation fails.

**Required fix**

- Build module-specific preflight hooks and run them after cheap CLI/config checks but before output mutation.
- Kreport preflight must find and version Python, compile/load the resolver, validate the taxonomy cache schema, open every assignment exactly as the resolver will, and perform a no-write resolution plan that applies the configured unresolved policy. Refresh mode may validate credentials/configuration without making network calls unless an explicit online preflight is requested.
- QC preflight must run complete bamstats parsing and `partition_minimap2_failures()` reconciliation, including required `C` reads and finite fields.
- Krona and FAPROTAX preflight must use the same dependency/version/shape checks as execution.
- Add paired tests: if validation succeeds on an immutable fixture, the corresponding real run must get past preflight; every fixture designed to fail during module setup must also fail validation with the same error class.

### P2 — required hardening

#### P2-01 — manifest v2 contradicts its documented array contract and lacks invariants

`analysis/00_run_pipeline.R:259-282` emits `NULL` for zero assignments and bamstats, which becomes JSON `null`; the README says input collections are always arrays even at cardinality zero. `analysis/utils/manifest.R:90-93` explicitly permits the contradiction. The scalar checker also accepts vectors of length greater than one, and element types are not validated.

The validator does not enforce the module state machine: a skipped/failed/not-run reason or error can be absent, a completed record can carry an error, duration can be negative/nonfinite, run status can disagree with module states, module keys can disagree with `cli.modules`, and declared outputs are not checked for uniqueness, containment, or existence.

**Required fix**

- Emit `assignments: []` and `bamstats: []`; never `null` for collection fields.
- Keep `schema_version: 2` as a conformance bug fix to the already documented v2 contract. Add `schema_revision: 1` if consumers need to distinguish corrected manifests.
- Commit a machine-readable JSON Schema or implement an equivalently exhaustive validator.
- Enforce scalar length/type, array element type, unique module/package/sample names, enum values, ISO-8601 UTC timestamps, finite nonnegative durations, and cross-field state/run invariants.
- Validate the in-memory object, write atomically, parse the emitted JSON back without simplification, validate again, then publish it.
- Verify every declared output exists, is a regular file within the owned output root, and occurs once.

#### P2-02 — “integer” gates still accept fractional values and unbounded work

`assert_scalar_number(..., integer=TRUE)` and upstream threshold validation use a floating tolerance. A value such as `1400.000000001` can pass and later reach integer-only formatting. Abundance validation uses an even wider tolerance, allowing fractional counts. `schema_version` is coerced with `as.integer()`, so values such as `1.5` or string-like inputs can be accepted. Several work multipliers have no operational upper bound, making a typo such as `1e12` a valid route to allocation failure or an impractically long run.

**Required fix**

- For all externally supplied integer fields, require numeric, finite, scalar, exactly `x == floor(x)`, and within a declared safe range. Do not use epsilon equality.
- Require `schema_version` to be the exact expected numeric/integer scalar without coercion.
- Apply explicit operational caps to rarefaction points/iterations, beta permutations/resampling iterations, depths, top-N values, and any allocation or loop multiplier. Put constants and rationales in one place.
- Either cap count values at the maximum supported by downstream integer conversions or remove those conversions and use exact-safe doubles up to the documented bound.
- Test `n + 1e-9`, negative zero where relevant, integer maximum boundaries, just-over-bound values, `Inf`, `NaN`, strings, booleans, and length-two vectors.
- Make `--modules` grammar strict: reject leading/trailing/repeated separators such as `qc,,alpha`, not silently normalize them.

#### P2-03 — valid identifiers can collide with generated column schemas

Sample IDs become wide-table column names; values such as `Taxon`, `TaxonPath`, or `SampleID` collide with explicit columns in composition/shared/beta products. Metadata columns such as `PC1`, `PCoA1`, or alpha metric names can cause `.x`/`.y` join suffixes and later failures or silent schema drift. Shared-taxa group labels become column names and can collide with `Taxon`, `TaxonPath`, or `Threshold`. Metadata headers are not centrally checked for emptiness, boundary whitespace, control characters, or reserved names.

**Required fix**

- Introduce one canonical identifier validator and per-output reserved-name registries.
- Prefer long-form outputs or an explicit reversible encoding/mapping over silently repairing names.
- Before every join, assert key cardinality and reject non-key name collisions; set suffixes explicitly only where they are intentionally consumed.
- Validate metadata headers, sample IDs, grouping values, and aggregate columns for boundary whitespace/control characters and TSV/JSON safety.
- Add adversarial tests for every reserved name and a positive test showing the mapping round-trips arbitrary allowed labels.

#### P2-04 — beta-diversity design and output provenance are inconsistent

`analysis/03_beta_diversity.R` selects Bray as primary when present and otherwise silently uses the first configured distance, yet `permanova.tsv` and `betadisper.tsv` do not identify that distance. Rarefaction stability always uses Bray. With `strata_column`, `adonis2` receives restricted strata but `permutest.betadisper` receives an unrestricted permutation count, so paired inference uses different randomization designs. Merely having two strata does not prove that any label-changing, non-identity permutation is possible.

**Required fix**

- Add required config `beta.primary_distance`; require it to occur exactly once in `beta.distances`.
- Record distance, transform/binary setting, seed, requested/effective permutation count, strata column, block sizes, and design status in every inferential output.
- Construct one `permute::how(blocks=...)` design or one validated permutation matrix and use it for both PERMANOVA and betadisper permutation tests.
- Prove there is at least one admissible non-identity, label-changing permutation. Otherwise write a deterministic skipped diagnostic rather than a misleading p-value.
- Make rarefaction stability use the declared primary distance, or explicitly add and record a separate `stability_distance` setting.
- Capture statistical warnings into module diagnostics/manifest rather than blanket-suppressing them.

#### P2-05 — degenerate beta stability can turn a valid dataset into a module crash

Three or more identical sample profiles produce a zero distance matrix. Full PCoA has a skip path, but the stability reference/Procrustes path lacks an equivalent nonzero-rank gate; all iterations can fail and abort the module. A non-`NULL` but malformed/nonfinite/wrong-dimension `Yrot` is also not validated inside the iteration error boundary.

**Required fix**

- Before stability work, require a finite, nonzero primary-distance matrix and enough positive axes for the configured comparison.
- If not, emit `rarefaction_stability_skipped.tsv` with a machine-readable reason and keep the module's other valid products.
- Validate every aligned result inside `tryCatch`: numeric, finite, exact dimensions, and exact sample row names/order.
- Record per-iteration failure reason plus successful/attempted counts. Apply the configured minimum-success rule without leaking partial plots/tables.

#### P2-06 — assignment and compression contracts remain asymmetric

The R parser accepts `.gz` case-insensitively; Python uses `endswith(".gz")`. Classified rows with positive TaxID and blank lineage are accepted. The length field validates only the final pipe component, so malformed prefixes/multiple pipes can pass. Read IDs and lineage components lack one uniform boundary-whitespace/control-character rule. Duplicate canonical assignment paths can be configured for multiple samples.

**Required fix**

- Detect gzip in Python with `str(path).lower().endswith(".gz")`, or preferably by validated compression mode/magic bytes shared with R.
- Require a nonblank normalized lineage for positive TaxID. Define one explicit allowed representation for unclassified TaxID 0.
- Validate the full length token against the intended grammar, not just its suffix; document the exact accepted forms.
- Reject read IDs and lineage components with leading/trailing whitespace, tabs/newlines/control characters, or empty ranks.
- Resolve configured paths canonically and require one physical assignment file per sample.
- Keep R and Python conformance fixtures identical, including `.gz`/`.GZ`, Unicode, CRLF, empty fields, malformed pipes, and long lines.

#### P2-07 — upstream database identity is name-only and spoofable

The upstream contract accepts a supported `database_set` name but does not require and pin the selected `database_sets[[name]]` mapping. A params file can therefore claim a supported label while omitting or replacing reference/ref2taxid/taxonomy resources.

**Required fix**

- Require the selected database-set object and all identity-bearing fields used by the supported wf-16s contract.
- Compare normalized identifiers/URIs and, where upstream supplies them, checksums against a versioned allowlist for the two supported NCBI sets.
- Reject top-level/nested disagreement, missing selected mapping, custom reference overrides, and unexpected classifier/rank combinations.
- Record the complete selected upstream resource identity and checksums in the manifest.
- Put accepted contracts in a data file with tests, not dispersed string checks.

#### P2-08 — TaxIDs and large counts are not safe across R/Python/JSON

Python permits arbitrary nonnegative integers, while `analysis/07_kreport_pavian.R` coerces TaxIDs with `as.integer()`. Values above 2^31-1 become `NA`; values above JSON's exactly representable double range can be silently rounded by R. Assignment TaxIDs currently enter R numerically. Some count paths likewise convert to R integer.

**Required fix**

- Treat TaxIDs as canonical decimal strings end-to-end in R, Python cache JSON, TSVs, and provenance. Use numeric conversion only after a documented safe-bound assertion where an API truly requires it.
- Reject signs, decimals, exponent notation, leading/trailing whitespace, and noncanonical leading zeros (except `0`).
- Choose and document the accepted TaxID domain; test `2147483647`, `2147483648`, `9007199254740991`, and the first rejected boundary.
- Define the maximum supported abundance/count value and ensure every downstream calculation/formatter preserves it exactly enough for the declared semantics.

#### P2-09 — taxonomy refresh commits source state before run acceptance

`analysis/utils/ncbi_taxonomy.py:259-306` writes the source cache when query transport had no failures, before applying `unresolved_policy=error`. The reproduced probe showed a failed run mutating the source cache. Exact-name queries also accept a unique TaxID without checking rank or ancestor context, and placeholder names can be queried.

**Required fix**

- Build a candidate cache without touching the source.
- Apply schema, TaxID-domain, ambiguity, ancestor/rank-context, and unresolved-policy gates to the candidate.
- Write run-local diagnostics/provenance first; commit the source cache atomically only if the module is accepted.
- On any failed module, guarantee byte-for-byte preservation of the source cache.
- Do not query placeholders such as `Unknown`/unclassified markers. For NCBI results, fetch/verify scientific name, rank, and lineage context before accepting a candidate TaxID.
- Record rejected candidates and reasons without installing them.

#### P2-10 — source revision provenance ignores dirty or archive execution

The manifest records `git rev-parse HEAD` but not working-tree modifications. Modified code can therefore report the clean commit ID. A source archive without `.git` loses provenance altogether.

**Required fix**

- Record `git_commit`, `git_dirty`, and a deterministic digest of maintained runtime/config schema files.
- Default release/publication mode to fail on a dirty tree; permit `--allow-dirty` only with `environment/source` warnings and the source digest.
- For source archives, record the deterministic digest even when `git_commit` is null.

### P3 — robustness, CI, and repository hygiene

#### P3-01 — dependency bootstrap can operate on the wrong project

`analysis/install_packages.R --restore` locates the repository lock but calls `renv::restore()` and `renv::status()` without explicit project/lockfile arguments. From another working directory this can target the wrong project. The “runner works outside the repository” test is launched from an R session already activated for this project and therefore does not prove a clean external invocation.

Fix the installer to activate/pass the repository root explicitly; make `--install` either a clearly named unlocked-development mode or remove it from the release path. Add a subprocess test from another directory with a sanitized R environment/library.

#### P3-02 — module result envelopes are trusted too late

The runner treats any result not equal to `failed` or `skipped` as completed. A misspelled/missing status or invalid output list is caught, if at all, only during final manifest work.

Validate each module return immediately: named object, allowed status, state-required fields, character arrays, unique contained outputs, and output existence. Convert an internal contract violation to a failed module record and obey fail-fast/keep-going normally.

#### P3-03 — warnings are intentionally erased

Several statistical calls are wrapped in `suppressWarnings()`, defeating the manifest's warning provenance. Capture warnings with handlers, classify known benign warnings explicitly, and write the rest to per-module diagnostics and the manifest.

#### P3-04 — FAPROTAX result shape is under-validated

The module checks that prepared IDs are present but not that raw binary-table row names are exact and unique. Reject missing, extra, duplicate, reordered-without-explicit-reindex, nonbinary, nonfinite, or wrong-dimension results. Require function names to be nonempty, unique, and TSV-safe. Preserve FAPROTAX 1.2.12 and existing non-exclusive semantics.

#### P3-05 — whitespace CI examines only the last commit

`tests/check_committed_whitespace.py` uses `HEAD^..HEAD`, so a multi-commit pull request can introduce whitespace before its final commit and pass. Use the pull-request merge base through `HEAD`, or scan the maintained source set. Exclude intentional generated/vendor/legacy artifacts explicitly. Also run `git diff --check <merge-base>...HEAD`.

#### P3-06 — tracked generated output is stale and leaks machine-local provenance

`output_Helga/` contains a schema-v1/pipeline-v0.1.0 manifest plus local absolute paths and old package/session data. It is not a current v0.4.0 fixture. Remove generated runtime output from version control, ignore it, and retain only deliberately minimized portable fixtures. Do not delete the tracked real Ambar regression inputs.

#### P3-07 — small release/provenance strictness gaps

- Require `VERSION` to exist and contain exactly one newline-terminated SemVer value; do not ignore extra lines or produce a length-zero condition.
- Store file mtimes as UTC ISO-8601 with `Z`, not locale/timezone-ambiguous strings.
- Expand the metadata test to find every active current-version declaration, not only the first README line.
- For v0.4.1 update `VERSION`, `CITATION.cff`, the five active README declarations, and prepend (do not rewrite) `CHANGELOG.md`.

## 5. Implementation plan

Implement in this order so later phases can rely on earlier invariants.

### Phase 1 — early, side-effect-free preflight

**Files:** `analysis/00_run_pipeline.R`, `analysis/utils/cli.R`, `analysis/utils/config.R`, new `analysis/utils/preflight.R`, `analysis/install_packages.R`.

1. Parse only CLI, version, and configuration first; validate strict module grammar before input parsing.
2. Validate exact config types/ranges and output-root safety.
3. Verify source cleanliness/digest and locked environment.
4. Resolve, deduplicate, fingerprint, and preferably snapshot inputs.
5. Run requested-module preflight hooks.
6. If `--validate-only`, print a structured summary and exit zero without creating or modifying any path.

Exit nonzero at the first preflight failure and include a stable error code such as `E_LOCK_MISMATCH`, `E_INPUT_CHANGED`, `E_OUTPUT_UNSAFE`, or `E_KREPORT_PREFLIGHT` for regression tests.

### Phase 2 — strict external contracts

**Files:** `analysis/utils/config.R`, `analysis/utils/io.R`, new contract data/schema files, tests/fixtures.

- Replace tolerant integer checks with exact bounded checks.
- Add operational caps and strict config mapping/type validation.
- Centralize identifier/header/reserved-name validation.
- Tighten assignment grammar and make R/Python compression behavior identical.
- Validate canonical path uniqueness.
- Pin the full selected upstream database contract.
- Preserve TaxIDs as strings and set explicit count bounds.

All errors must identify file, physical line/field when applicable, rejected value, and expected contract without dumping whole sensitive records.

### Phase 3 — transactional modules and owned output publication

**Files:** runner, every module writer, new `analysis/utils/atomic_io.R` and `analysis/utils/module_result.R`.

- Give each module a staging root and require it to return paths relative to that root.
- Validate the result envelope before publishing.
- Rename staged outputs only after module success/intentional skip.
- Maintain an ownership inventory in the manifest.
- Define exact overwrite behavior for prior owned outputs and protect all unowned files.
- Publish the manifest last and atomically.

Use failure-injection tests after each output write to prove no partial/stale products become part of a committed run.

### Phase 4 — manifest v2 conformance

**Files:** `analysis/utils/manifest.R`, committed schema, runner, `tests/testthat/test-manifest.R`, release verifier.

- Make every documented collection an array at zero/one/many cardinality.
- Enforce the full state machine and cross-field invariants.
- Record truthful environment/source/input states.
- Round-trip validate emitted JSON.
- Make the release verifier validate the schema and declared-output hashes/existence, not just selected values.

### Phase 5 — beta and statistical hardening

**Files:** `analysis/03_beta_diversity.R`, relevant configuration/defaults/docs/tests; review warning handling in modules 02/05/06 as well.

- Declare and record primary/stability distance.
- Share one validated permutation design between PERMANOVA and betadisper.
- Gate unexchangeable strata and degenerate distances.
- Validate Procrustes iteration results inside the iteration boundary.
- Preserve warning/provenance data.

The vegan implementation supports passing a `permute::how` control or explicit permutation matrix to the relevant permutation machinery; use a single materialized design so the two tests cannot diverge accidentally.

### Phase 6 — taxonomy/FAPROTAX transactions and provenance

**Files:** `analysis/07_kreport_pavian.R`, `analysis/08_faprotax.R`, `analysis/utils/ncbi_taxonomy.py`, Python/R tests.

- Add resolver validation-only mode.
- Fix case-insensitive gzip handling.
- Make refresh candidate-based and commit-after-acceptance.
- Verify fetched rank/lineage context and exclude placeholders.
- Use string TaxIDs end-to-end.
- Strengthen FAPROTAX shape/name validation.

### Phase 7 — CI, release metadata, and cleanup

**Files:** `.github/workflows/ci.yml`, test helpers, whitespace checker, README, changelog, citation, `VERSION`, `.gitignore`.

- Test clean external invocation and lock mismatch failure on Ubuntu and Windows.
- Add zero/one/many manifest integrations and default/subset/keep-going/overwrite failure-injection runs.
- Correct PR-range whitespace checks.
- Remove stale generated `output_Helga` artifacts after preserving any deliberately reusable fixture data elsewhere.
- Bump all release metadata to 0.4.1 in one final commit after behavioral work is green.

## 6. Required regression matrix

### Preflight and provenance

- missing lock; mismatched R; one missing package; one wrong package version; synchronized lock;
- dirty source default rejection and `--allow-dirty` provenance;
- input changed between assignment passes, before Python re-read, and before manifest publication;
- duplicated assignment path via different relative/symlink spellings;
- validate-only with missing Python, malformed cache, uppercase `.GZ`, unresolved-policy error, malformed `C` bamstats row, and missing Krona executable when rendering is requested;
- unsafe output roots and input/output overlap.

### Transaction/overwrite

- full run followed by subset overwrite leaves no stale unrequested pipeline products;
- failure after each module's first and last staged write publishes no partial products;
- failed overwrite preserves the previous committed run;
- non-empty unowned directory is rejected and unowned sentinel files are never removed;
- malformed module result status/output path is captured as a module failure.

### Manifest

- zero/one/many samples, assignments, bamstats, warnings, outputs, packages, and Krona samples remain arrays after JSON round trip;
- negative matrix for each invalid module-state combination;
- run-status/module-status mismatch;
- negative/nonfinite duration; invalid timestamp; length-two scalar; wrong array element type;
- output outside root, duplicate output, missing output, and wrong hash;
- not-run records for fail-fast and nonrequested modules exactly match the documented semantics.

### Input/config

- near integers (`1 + 1e-9`), just-over-cap workloads, wrong scalar types, length-two values, malformed nested mappings;
- `qc,,alpha`, `,qc`, `qc,`, duplicate and unknown modules;
- reserved sample/group/metadata names and join collisions;
- blank positive-TaxID lineage, empty rank, padded/control-character read IDs, malformed full length token;
- `.gz`, `.GZ`, CRLF, Unicode, long line, bad UTF-8, and gzip corruption in both R and Python;
- supported database names with missing/mismatched nested resources and top-level overrides;
- TaxID/count boundary values.

### Statistics/taxonomy

- jaccard-only configuration records jaccard everywhere;
- stratified design uses the same materialized permutations for PERMANOVA and betadisper;
- singleton/confounded strata yield an explicit skip, not a nominal p-value;
- identical profiles and rank-one distance matrices skip stability cleanly;
- malformed/nonfinite Procrustes return is isolated and diagnosed;
- taxonomy refresh success commits once; transport failure, ambiguity, context mismatch, and unresolved-policy failure preserve source bytes;
- FAPROTAX raw result with missing/extra/duplicate/nonbinary rows or unsafe function names fails before final output.

## 7. Release acceptance gates

v0.4.1 is releasable only when all of the following hold on the exact candidate commit:

1. The repository is clean and the active environment reports `lock_status: synchronized` on R 4.5.3.
2. All maintained R and Python sources parse/compile.
3. All R/testthat and Python unit/regression tests pass on Ubuntu and Windows.
4. The default tracked-fixture integration, synthetic minimap2 integration, and opt-in FAPROTAX integration pass and are verified against manifest/output invariants.
5. Validate-only parity, transaction failure injection, zero/one/many manifest, uppercase gzip, refresh rollback, identifier-collision, exact-integer, and constrained-permutation tests pass.
6. Baseline comparison confirms no unintended changes to read denominators, valid lineage normalization, kreport/Krona counts, or FAPROTAX 1.2.12 results.
7. `git diff --check <merge-base>...HEAD`, the maintained-source whitespace gate, JSON/YAML parsing, fixture gzip checks, and repository integrity checks pass.
8. `VERSION`, `CITATION.cff`, README declarations, changelog, manifest expectation tests, and release verifier all agree on `0.4.1`.
9. CI is green for all four existing matrix executions: main suite on Ubuntu/Windows and FAPROTAX suite on Ubuntu/Windows.
10. Only after the above, create the signed/annotated v0.4.1 tag at the verified commit and rerun release verification from a clean checkout of that tag.

## 8. Suggested v0.4.1 changelog entry

### Fixed

- Make locked-environment provenance reflect the active R/package library rather than lockfile presence.
- Detect input mutation and prevent mixed-revision analyses.
- Make overwrite and module publication ownership-aware and transactional, eliminating stale and partial products.
- Make validate-only execute the same module preflight contracts as a real run.
- Conform zero-cardinality manifest fields to schema-v2 arrays and enforce module/run state invariants.
- Reject near-integer external values, unsafe resource requests, output-name collisions, malformed assignments, spoofed upstream database metadata, and unsafe TaxID/count conversions.
- Use one recorded permutation design for beta tests and handle degenerate stability inputs without aborting the module.
- Preserve the taxonomy source cache on every failed refresh and make gzip handling consistent across R and Python.

### Changed

- Record synchronized lock status, source cleanliness/digest, immutable input fingerprints, explicit beta distance/permutation design, and transactional output ownership in the run manifest.
- Strengthen CI with external-invocation, failure-injection, manifest cardinality/state, taxonomy rollback, and pull-request-range whitespace gates.

## 9. Handoff note for Luna

Treat P1-01 through P1-04 as indivisible release blockers. Implement Phases 1–4 before touching statistical/taxonomy refinements because they establish the safe execution and publication boundary that the later modules should use. Keep each new rejection behind a focused fixture and stable error code. If a new validation rule changes a previously accepted real fixture, stop and document the exact row/value and scientific reason rather than silently normalizing it.
