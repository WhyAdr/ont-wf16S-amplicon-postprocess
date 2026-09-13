# ONT wf16S postprocess v0.4.7 patch — stress-tested revision

Status: implementation-ready hardening and Pavian-integration plan
Baseline audited checkout: `main` at `8d543e84742302bcc17b0275264bfe48de75cae0`
Existing release tag: annotated `v0.4.6`, targeting `8d543e84742302bcc17b0275264bfe48de75cae0`
Original plan reviewed: `wf16s-postprocess-v0.4.7-patch(1).md`
Refinement date: 2026-09-12 UTC / 2026-09-13 WIB

> This revision supersedes the attached draft where they differ. It keeps the draft's strong conclusions, but closes several architectural contradictions found while stress-testing the current implementation.

## Executive decision

The September 10–11 post-v0.4.6 fixes remain correct and should stay:

1. Python and R now use compatible POSIX record-lock semantics for the taxonomy cache.
2. Python has an explicit process-local re-entrancy guard.
3. The release verifier builds classified-count inputs independently of the optional alpha module, so composition-only and mixed zero-classified runs are verifiable.

The branch is nevertheless **not yet a provenance-grade v0.4.7 release candidate**. The original draft correctly identified lock cleanup, provenance, release identity, and Pavian scope as the main remaining work, but the stress test found six additional design corrections that should be treated as part of the v0.4.7 contract:

- **Lock identity must be tied to the canonical cache target, not merely the spelling of `<cache>.lock`.** A final-component symlink alias can otherwise create a second lock pathname even though both callers modify the same cache file.
- **Runtime source provenance and release-contract provenance should be separated.** The run manifest should inventory producer/runtime files with per-file SHA-256 values; the exact Git commit/tag already binds tests, CI, README, CHANGELOG, and CFF. Do not bloat every run manifest with the entire release harness just to prove the release harness existed.
- **The exact-tag release gate needs two phases.** It is impossible to require `HEAD` to be an exact `v0.4.7` tag *before* creating that tag. Run branch/RC gates first, create the annotated tag locally, run the exact-tag gate locally, then push the tag and let tag CI confirm the same identity before publishing the GitHub Release.
- **A Pavian viewer cannot derive conflict status from `.kreport` alone.** `taxonomy_resolution.tsv` is the authoritative per-sample/per-path sidecar for `Resolved`, `Unresolved`, and `Conflicted`. The viewer builder and verifier must consume it alongside `.kreport`.
- **`taxonomy_resolution.tsv` must become unconditional kreport output.** The current module writes it only when at least one classified node exists. An all-unclassified sample/run therefore has no resolution sidecar, which would make the proposed Pavian contract impossible to verify. v0.4.7 should always emit the header, even when it has zero rows.
- **Do not brand the original static renderer as “Pavian.”** Keep `.kreport` as the official Pavian interoperability artifact, and call the in-repo renderer something explicit such as `builtin_kraken_report_explorer`. Its provenance must say `official_pavian_compatibility = "kraken_report_input_contract_only"`.

With those corrections, v0.4.7 can be a clean stability + provenance + visualization release rather than a second post-release repair cycle.

---

## 1. Verified baseline and release identity

Current repository state was rechecked against GitHub:

| Item | Verified state | Consequence |
|---|---|---|
| `main` | `8d543e84742302bcc17b0275264bfe48de75cae0` | Same checkout assumed by the attached draft. |
| HEAD message | `fix: verify composition without alpha module` | The no-alpha verifier repair is the current branch tip. |
| `v0.4.6` | Annotated tag object `8e78b0d...` | It is not a lightweight tag. |
| `v0.4.6^{commit}` | `8d543e84742302bcc17b0275264bfe48de75cae0` | The tag points to the post-release repaired checkout. |
| tagger timestamp | `2026-09-11T01:22:06Z` | Later than the September 10 release metadata. |
| `VERSION` / CFF / CHANGELOG | `0.4.6`, September 10 metadata | Keep v0.4.6 immutable; document chronology instead of retagging. |

The current `v0.4.6` state is therefore a legitimate historical tag identity, albeit one whose tagger timestamp and human release date differ. **Do not move or recreate that tag.**

### v0.4.7 release identity contract

For v0.4.7:

1. `VERSION` must contain exactly `0.4.7\n`.
2. `CITATION.cff` top-level `version` must be `0.4.7`.
3. `CITATION.cff` `date-released` and the first non-`Unreleased` CHANGELOG heading must agree on one `YYYY-MM-DD` date.
4. The release tag must be **annotated** and named exactly `v0.4.7`.
5. `git rev-parse 'refs/tags/v0.4.7^{commit}'` must equal `git rev-parse HEAD` during the local tag gate and tag CI.
6. The tagger timestamp should normalize to the declared release date in **UTC**. This avoids local-timezone ambiguity near midnight.
7. Tag CI must verify the run manifest's `git_commit == HEAD` and all release flags are strict.
8. Only after tag CI passes should the GitHub Release object be published.

Do **not** require the ordinary branch-level `tests/verify_release_run.R` invocation to have an existing tag; that would make pre-tag CI impossible. Exact tag identity belongs in `tests/check_release_metadata.py --require-tag ...` and the tag-push workflow.

---

## 2. Taxonomy-cache locking: retain the fix, close the interruption and alias holes

### 2.1 What is already correct

`analysis/utils/ncbi_taxonomy.py` now uses `fcntl.lockf()` on POSIX, matching the record-lock family used by R's `filelock`. The process-local `_ACTIVE_CACHE_LOCKS` set is also necessary because POSIX record locks are process-associated and are not themselves a reliable re-entrancy guard.

Keep both behaviors.

### 2.2 P1 — interrupted acquisition can poison the process-local guard

The current order is:

1. normalize lock identity;
2. add identity to `_ACTIVE_CACHE_LOCKS`;
3. poll for the OS lock;
4. only remove the identity on ordinary timeout or normal context exit.

`KeyboardInterrupt`, `SystemExit`, an injected exception from `time.sleep`, or a non-`OSError` failure between registration and successful acquisition can bypass both cleanup sites. In a reusable Python process, later attempts then fail with a false “already held by this process.”

#### Required implementation

Use one outer `try/finally` covering **registration, polling, owner metadata, yield, unlock, and close**. Track `registered` separately so a failed re-entrant attempt never removes the identity belonging to the first holder.

Also:

- use `time.monotonic()` for deadlines;
- ensure each `candidate_fd` is closed in a nested `finally` unless it becomes the acquired descriptor;
- keep unlock/close best-effort and idempotent;
- do not unlink the lock file;
- write owner metadata with a complete-write helper rather than assuming one `os.write()` consumes the whole buffer.

### 2.3 P1 — canonical cache target, not merely canonical lock spelling

`os.path.realpath(cache_path + ".lock")` does not necessarily collapse a **final-component symlink of the cache itself**. For example, `alias.json -> /real/cache.json` yields `alias.json.lock`, while `/real/cache.json` yields `/real/cache.json.lock`.

The same risk exists on the R side because `get_taxonomy_lock_path()` currently builds the lock path from `canonicalize_root_path(cache_path)`, which intentionally canonicalizes the parent and basename rather than resolving the final file target.

For v0.4.7, define one cache identity rule in both runtimes:

- the taxonomy cache must already exist before refresh/cache-only resolution;
- resolve the **cache file itself** to its real target;
- append `.lock` and the transaction-journal suffix to that real target identity;
- retain the user-supplied/configured path separately for provenance and messages.

This makes relative paths, parent-directory symlinks, and final-component symlinks converge on the same lock/journal identity.

Hard-link aliases are a different filesystem identity problem and need not be solved in v0.4.7; document that the cache should not be accessed through multiple hard links.

### 2.4 Required lock tests

Add or retain all of the following:

- blocked acquisition interrupted by injected `KeyboardInterrupt`; after competing process release, the same process acquires successfully;
- same test with injected `SystemExit` or generic `RuntimeError` from the poll/sleep path;
- two Python threads contend for the same cache; the second fails quickly and does not poison the first;
- relative path versus absolute path alias;
- symlinked parent-directory alias on POSIX;
- final-component cache symlink versus real cache path on POSIX;
- existing R-holds/Python-waits and Python-holds/R-waits tests;
- Windows coverage retained, with symlink-specific tests skipped when the runner cannot create symlinks.

---

## 3. Source provenance: make runtime identity auditable without conflating it with the release harness

### 3.1 Current weakness

`source_provenance()` computes SHA-256 values internally and folds them into `source_digest_sha256`, but `run_manifest.json` records only the aggregate digest. The returned `source_files` path list is not emitted by `analysis/00_run_pipeline.R`.

An auditor therefore cannot answer “which file hash contributed to this digest?” from the manifest alone, and the release verifier does not independently recompute the source digest.

### 3.2 Refined v0.4.7 contract

Add an ordered `source_files` array to manifest revision 3. Each element is an object:

```json
{
  "path": "analysis/utils/ncbi_taxonomy.py",
  "sha256": "...64 lowercase hex..."
}
```

Rules:

- path is repository-relative POSIX syntax;
- no backslashes, absolute paths, empty components, `.` or `..`;
- reject control characters;
- reject duplicate and case-fold-colliding paths;
- reject symlinked producer files in the stable release checkout;
- entries are sorted deterministically by repository path;
- aggregate digest is SHA-256 over the exact UTF-8 sequence `path<TAB>sha256`, joined by `\n`, with **no final newline**;
- `source_digest_sha256` remains for compact identity, but it must equal the digest recomputed from `source_files`.

### 3.3 Producer-source scope versus release-contract scope

Do **not** solve release provenance by stuffing every test, fixture, workflow, README, and CHANGELOG into every run manifest.

Use two layers:

**A. Manifest producer-source inventory** — files capable of affecting pipeline output or its in-repo runtime renderers, preserving the existing conceptual scope:

- `analysis/**` maintained R/Python/JSON/YAML files;
- `analysis/vendor/krona-2.8.1/**` pinned assets;
- `config.example.yml`;
- `VERSION`;
- `renv.lock`;
- `.Rprofile`, `renv/activate.R`, `renv/settings.json`.

**B. Release-contract identity** — tests, fixtures, CI workflow, README, CHANGELOG, CFF, and release-check scripts are bound by:

- exact clean `git_commit` recorded in the run manifest;
- exact local/tag-CI tag-to-commit mapping;
- independent release metadata checks.

The Git commit already binds the complete tracked tree. Recording every release-harness file again in each run manifest adds complexity without materially strengthening the trust chain.

### 3.4 Independent verification requirement

`tests/verify_release_run.R` must not call `source_provenance()` and declare victory.

It should independently:

1. derive the expected producer-source path set from Git and hard-coded release scope rules;
2. compare that expected set to `manifest$source_files` exactly;
3. hash each current checkout file independently;
4. compare every per-file hash;
5. recompute the canonical aggregate digest;
6. compare it to `manifest$source_digest_sha256`.

This catches omitted, extra, reordered, replaced, or hash-mismatched source records.

---

## 4. Manifest schema v2 revision 3

v0.4.7 changes the manifest contract enough that it should be unmistakably revision 3. **Do not silently make revision-2 keys optional or mutate revision 2 in place.** Existing v0.4.6 runs should continue validating under revision 2.

### 4.1 New/changed revision-3 fields

Require:

- `source_files`: ordered path/SHA-256 objects;
- `cli.allow_large_workload`: logical;
- `cli.pavian`: logical;
- `exports`: explicit resolved export state for Krona and the built-in Pavian integration;
- existing revision-2 artifacts/ownership fields unchanged.

Recommended `exports` shape:

```json
{
  "krona": {
    "enabled": true,
    "render_html": true,
    "provenance_path": "07_Kreport/krona/krona_provenance.json"
  },
  "pavian": {
    "enabled": true,
    "render_html": true,
    "provenance_path": "07_Kreport/pavian/pavian_provenance.json",
    "integration": "official_pavian_upload_plus_builtin_kraken_report_explorer",
    "official_pavian_compatibility": "kraken_report_input_contract_only"
  }
}
```

When disabled, `provenance_path` is `null` and no **owned** output for that export may appear in `modules$kreport$outputs`.

Do not reject a user-preserved *unowned* file merely because its name resembles a Pavian output; ownership and `preserved_unowned_outputs` remain authoritative.

### 4.2 Release flag policy

The release verifier must require:

- `cli.allow_dirty == false`;
- `cli.allow_unlocked == false`;
- `cli.allow_large_workload == false`;
- `cli.refresh_taxonomy == false`;
- `cli.online_preflight == false`;
- `environment.locked == true`;
- `environment.lock_status == "synchronized"`.

`validate_only`, `keep_going`, and `overwrite` are execution-mode fields, not universal release identity fields. Verify them only where the specific release fixture requires a value; do not invent an unnecessary repository-wide prohibition.

### 4.3 Backward compatibility

Add `validate_manifest_v2_revision3()` and update the dispatcher to accept revisions 1, 2, and 3.

Prefer keeping revision 2 frozen. Shared leaf validators may be refactored, but do not implement revision 3 by temporarily rewriting `schema_revision` to `2` and calling the old validator; that creates brittle hidden coupling.

Update `manifest_is_valid_run()` and prior-output ownership migration to accept revision 3.

---

## 5. Reproducibility contract: decide it now

The repository should explicitly claim **semantic reproducibility**, not universal byte-for-byte reproducibility.

v0.4.7 stable-release contract:

- exact pipeline source commit and per-file producer hashes;
- exact `renv.lock` hash and synchronized R environment;
- exact package versions and installed package identity metadata already captured by the manifest;
- Python 3.12 release requirement;
- input hashes;
- taxonomy cache/resolved-cache hashes;
- deterministic seeds and analysis settings;
- deterministic normalized Pavian JSON/HTML generated by the stdlib builder.

Do **not** claim that all PNGs, R text output, fonts, device rendering, or Windows/Linux file bytes are globally identical unless the project later pins containers/fonts/graphics devices and tests those bytes.

For the new built-in viewer, however, byte determinism is realistic because Python can write explicit UTF-8 bytes with `\n` semantics. Require repeated rendering from the same canonical inputs to produce identical `.pavian.json` and `.pavian.html` hashes on each target OS. A cross-OS hash-comparison job is desirable but can be P2 if CI artifact exchange would materially complicate v0.4.7.

---

## 6. Transaction and journal hardening

The attached draft correctly identified weak publication-journal identity. Tighten both taxonomy and publication journals without breaking recovery semantics.

### 6.1 New journal schema

All **new** journals written by v0.4.7 should include:

- `journal_schema_version: 2`;
- required valid `transaction_id`;
- explicit phase from a closed enum;
- canonical target identity;
- for publication journal, `had_prior` boolean;
- stage/backup paths constrained to the expected parent and filename prefix.

Publication phases:

- `prepared`;
- `prior_moved`;
- `stage_published`.

Taxonomy phases remain the maintained equivalents, e.g. `prepared`, `candidate_committed`, `output_published`.

### 6.2 Legacy journal behavior

A legacy journal without schema/transaction identity may exist after a v0.4.6 crash. Do not destructively “upgrade” it by guesswork.

Recommended behavior:

- parse legacy journal conservatively;
- validate target/path/hash facts that can be established;
- if transaction identity is absent where two completed states could be confused, fail closed with `E_*_RECOVERY_REQUIRED` and retain all recovery material;
- document manual recovery rather than deleting an ambiguous journal.

### 6.3 Provenance writes

Replace direct `jsonlite::write_json()` for `krona_provenance.json` with `atomic_write_json()` and use the same path for Pavian provenance. Hash provenance only after atomic publication succeeds.

---

## 7. Pavian integration: exact boundary

### 7.1 What “built-in Pavian” means in v0.4.7

The official Pavian project is an R/Shiny GPL-3 application. Its parser accepts Kraken-style six-column reports, but its UI is not a small static renderer that can be silently vendored into this MIT repository.

Therefore v0.4.7 should expose two clearly separated capabilities:

1. **Official Pavian interoperability** — the canonical `.kreport` files remain the files users can upload/open in upstream Pavian.
2. **Built-in offline explorer** — an original, dependency-free, self-contained static viewer generated from the same `.kreport` plus the pipeline's taxonomy-resolution sidecar.

Use a renderer identity such as:

```text
builtin_kraken_report_explorer
```

not `official_pavian`, not `builtin_pavian`, and preferably not `pavian_compatible` in places where a reader might infer UI/API equivalence.

Provenance should say:

```json
"official_pavian_compatibility": "kraken_report_input_contract_only"
```

The output directory may remain `07_Kreport/pavian/` because it groups the official interoperability and convenience-viewer feature, but the renderer identity must stay unambiguous.

### 7.2 Do not add a ninth pipeline module

Keep Pavian under the existing `kreport` module, exactly as Krona is an export subextension. The fixed module registry remains eight modules.

### 7.3 Configuration

Add the smallest useful configuration surface:

```yaml
pavian:
  enabled: false
  render_html: true
```

Do **not** add a configurable `renderer: builtin` key while only one provider exists. Every configuration key becomes a public contract; avoid a no-op option until an actual second renderer exists.

Add CLI:

```text
--pavian
```

which sets `pavian.enabled = true`.

The dependency check belongs where the runner already checks Krona after resolving `requested_modules`:

```text
Krona export requires the 'kreport' module.
Pavian export requires the 'kreport' module.
```

This is not a `validate_config()` responsibility because `requested_modules` is computed after YAML/CLI merging.

### 7.4 Behavior

- `kreport` requested, Pavian disabled: canonical `.kreport` and normal taxonomy sidecars only.
- Pavian enabled, `render_html=false`: emit deterministic normalized JSON + Pavian provenance, no HTML.
- Pavian enabled, `render_html=true`: emit JSON + self-contained HTML + provenance.
- `--validate-only --pavian`: validate configuration, builder availability, Python availability, and input/resource contracts without creating any output directory.
- built-in rendering must never perform taxonomy lookup or network access.

---

## 8. Canonical `.kreport` and resolution-sidecar contract

### 8.1 `.kreport` remains the official interoperability artifact

The existing writer emits:

1. `U` unclassified line;
2. `R` classified root line;
3. DFS-ordered taxonomic nodes with two spaces of indentation per depth.

Keep that shape.

### 8.2 Harden the writer against ambiguous labels

Before emitting a report, reject node names containing:

- tab;
- CR or LF;
- NUL or other ASCII control characters;
- ambiguous leading indentation that would be interpreted as hierarchy by a Kraken/Pavian parser.

Do not silently trim or rewrite biological labels at rendering time. Fail with an actionable report-format error.

### 8.3 Always emit `taxonomy_resolution.tsv`

The current code writes `taxonomy_resolution.tsv` only when `length(resolution_rows) > 0`.

Change it so the file always exists with exactly these columns:

```text
SampleID	Depth	RankCode	NodeName	TaxonPath	TaxID	Status	ResolutionSource
```

An all-unclassified run therefore produces a header-only file, which is deterministic and verifiable.

### 8.4 Why the viewer needs the resolution sidecar

`.kreport` contains TaxID but does **not** encode the pipeline's `ResolutionSource` or the distinction between:

- resolved node;
- unresolved TaxID 0;
- assignment conflict resolved by the documented modal/minimum-TaxID tie-break.

A `Conflicted` row may carry a nonzero winning TaxID. Therefore the built-in explorer must not infer conflict status from TaxID alone.

Builder inputs should be:

```text
--kreport PATH
--resolution-tsv PATH
--sample-id ID
--json-out PATH
--html-out PATH        # only when requested
--expected-total INTEGER
```

The builder must join `taxonomy_resolution.tsv` rows for the requested sample to the classified report nodes by **exact `TaxonPath`**, and cross-check `Depth`, `RankCode`, `NodeName`, and `TaxID`.

No extra sidecar row and no unaccounted classified report node is allowed.

---

## 9. Independent Kraken-report parsing rules

The builder and correspondence verifier must each implement their own parser logic. They may share the written **specification**, but the verifier must not import parsing/building functions from the producer module.

### 9.1 Special records

Require exactly:

- first non-empty line: `U`, TaxID `0`, name `unclassified`, depth 0;
- second line: `R`, TaxID `1`, name `root`, depth 0;
- `U.clade == U.direct == unclassified_count`;
- `R.direct == 0`;
- `total == U.clade + R.clade`.

A classified root of zero with no classified children is valid when the sample is entirely unclassified.

### 9.2 Classified-node records

For all remaining rows:

- exactly six tab-separated fields;
- finite non-negative integer clade and direct counts, bounded by `2^53 - 1`;
- canonical unsigned-decimal TaxID spelling, also bounded by `2^53 - 1`;
- allowed rank codes `D,K,P,C,O,F,G,S`;
- rank code must match depth exactly for this fixed eight-rank pipeline (`D=1`, `K=2`, ..., `S=8`);
- indentation is exactly two ASCII spaces per depth;
- depth may not jump by more than one level;
- every non-root classified node has exactly one parent in the active stack;
- no duplicate path;
- no duplicate sibling name/path;
- `clade == direct + sum(child.clade)` for every node;
- root classified clade equals the sum of top-level classified clades;
- sum of all classified direct counts equals `R.clade`;
- `classified direct + unclassified == total`.

### 9.3 Percentage column

Do not treat the report percentage as an authoritative abundance source. Upstream Pavian itself recomputes percentages from counts.

Still validate that the field is a finite decimal in `[0, 100]` and is consistent with `clade / total` within the tolerance implied by the writer's two-decimal formatting. This catches gross corruption without creating cross-language rounding fragility.

### 9.4 Ordering

Require the pipeline's deterministic DFS order and sibling ordering contract rather than merely accepting any Kraken report order. That makes the normalized JSON stable and catches accidental writer drift.

The builder and verifier should independently check ordering. Do not derive one expected order by calling the producer's helper.

---

## 10. Normalized payload and static explorer

### 10.1 Output layout

```text
07_Kreport/
  <sample>.kreport
  taxonomy_resolution.tsv
  pavian/
    <sample>.pavian.json
    <sample>.pavian.html
    pavian_provenance.json
```

HTML is absent when `render_html=false`.

### 10.2 Canonical payload

Recommended payload schema:

```json
{
  "schema_version": 1,
  "sample_id": "Sample A",
  "renderer": "builtin_kraken_report_explorer",
  "renderer_version": "0.4.7",
  "official_pavian_compatibility": "kraken_report_input_contract_only",
  "count_model": "direct abundance-table taxon counts plus canonical unclassified count",
  "denominator": "TotalReads",
  "classified_definition": "sum of direct positive-count classified taxonomy rows",
  "totals": {
    "total": 1234,
    "classified": 1200,
    "unclassified": 34
  },
  "nodes": [
    {
      "path": "Bacteria;Bacillati",
      "parent_path": "Bacteria",
      "name": "Bacillati",
      "depth": 2,
      "rank_code": "K",
      "taxid": "1783272",
      "direct": 0,
      "clade": 1000,
      "status": "resolved",
      "resolution_source": "cache"
    }
  ]
}
```

Design choices:

- store TaxID as a **string**, because it is an identifier and this avoids JavaScript numeric semantics;
- store counts as integers bounded to JavaScript-safe integer range;
- keep nodes in canonical DFS order;
- no timestamps in payload;
- no absolute paths in payload;
- provenance, not payload, carries artifact paths/hashes.

### 10.3 Deterministic serialization

Use Python standard library only.

Canonical JSON bytes should use an explicitly documented serializer, for example:

- UTF-8;
- `allow_nan=False`;
- sorted object keys;
- compact separators `(',', ':')`;
- deterministic array order;
- one chosen final-newline policy, tested byte-for-byte.

The exact serializer is part of payload schema v1 and must be unit-tested.

### 10.4 Safe HTML embedding

Do **not** inject raw JSON text into executable JavaScript or rely on ad-hoc replacement of `</script>`.

Safer design:

1. serialize canonical JSON bytes;
2. Base64-encode those exact bytes;
3. embed the Base64 string in a non-executable data element;
4. decode it in inline JavaScript with `atob` + `TextDecoder`;
5. insert all biological labels with DOM `textContent`, never `innerHTML`.

The correspondence verifier can then Base64-decode the embedded payload and require exact byte identity with `<sample>.pavian.json`.

### 10.5 Browser network hardening

Add a restrictive CSP meta tag, e.g. conceptually:

```text
default-src 'none';
script-src 'unsafe-inline';
style-src 'unsafe-inline';
img-src data:;
connect-src 'none';
font-src 'none';
object-src 'none';
base-uri 'none';
form-action 'none'
```

Do not include `fetch`, `XMLHttpRequest`, `WebSocket`, `EventSource`, `sendBeacon`, CDN assets, remote fonts, remote stylesheets, or remote images.

“No network” should be verified by this explicit source/artifact contract; do not pretend a hosted CI runner has a universal cross-platform egress firewall when it does not.

### 10.6 Viewer scope for v0.4.7

Keep the first stable viewer intentionally smaller than upstream Pavian:

- sample title;
- total/classified/unclassified summary;
- collapsible taxonomy tree;
- direct/clade count toggle;
- percentage view;
- rank/depth filter;
- search;
- visible resolved/unresolved/conflicted state;
- table view with TaxID, direct, and clade counts;
- client-side download of the embedded normalized JSON.

Do **not** promise exact upstream Pavian Sankey/UI equivalence.

A custom Sankey/alluvial renderer is attractive but is a separate visualization algorithm with its own layout and test burden. Defer it unless it can be added without compromising the release hardening work. The official `.kreport` can always be opened in upstream Pavian for Pavian's native visualizations.

The existing sibling `.kreport` file is already published; the HTML does not need to embed a second raw copy merely to offer an “original kreport download” button.

---

## 11. Built-in viewer provenance

Write `07_Kreport/pavian/pavian_provenance.json` atomically.

Recommended shape:

```json
{
  "schema_version": 1,
  "path_basis": "run_dir",
  "integration": "official_pavian_upload_plus_builtin_kraken_report_explorer",
  "renderer": "builtin_kraken_report_explorer",
  "renderer_version": "0.4.7",
  "official_pavian_compatibility": "kraken_report_input_contract_only",
  "html_status": "rendered",
  "standalone_html": true,
  "count_model": "direct abundance-table taxon counts plus canonical unclassified count",
  "denominator": "TotalReads",
  "classified_definition": "sum of direct positive-count classified taxonomy rows",
  "samples": []
}
```

Each sample record should include:

- exact `sample_id`;
- `total_reads`, `classified_reads`, `unclassified_reads`;
- number of resolved, unresolved, and conflicted classified nodes;
- run-relative POSIX path + SHA-256 for `.kreport`;
- run-relative POSIX path + SHA-256 for `.pavian.json`;
- optional run-relative POSIX path + SHA-256 for `.pavian.html`, otherwise explicit `null` fields.

Do not put the provenance file's own hash inside itself. The manifest `artifacts` array already provides that external hash.

Sample record order must exactly follow `manifest$samples` / `context$samples`.

---

## 12. Independent correspondence verifier

Add a verifier such as:

```text
analysis/utils/verify_kraken_viewer_correspondence.py
```

Inputs:

```text
--kreport PATH
--resolution-tsv PATH
--json PATH
--html PATH            # optional when render_html=false
--sample-id ID
--expected-total INTEGER
```

The verifier must **not import** `kraken_report_viewer.py` or any producer parser.

It should independently:

1. parse and validate `.kreport`;
2. parse the sample's `taxonomy_resolution.tsv` rows;
3. validate exact path/depth/rank/name/TaxID/status correspondence;
4. load external JSON and validate its schema;
5. compare exact totals and ordered node records;
6. extract and Base64-decode the HTML payload;
7. require decoded bytes to equal the external JSON bytes exactly;
8. scan the HTML for prohibited remote resource/network constructs;
9. confirm exact sample identity;
10. fail on any extra or missing node, sample row, output, or provenance hash.

Mutation tests must cover:

- sample ID;
- total/classified/unclassified counts;
- direct count;
- clade count;
- parent/path;
- depth;
- rank;
- TaxID;
- `Resolved` ↔ `Unresolved` ↔ `Conflicted` status;
- resolution source;
- sibling order;
- duplicate sibling/path;
- U/root lines;
- malformed indentation;
- embedded payload bytes;
- external JSON bytes;
- injected remote `<script src>`, stylesheet, image, `fetch()`, XHR, websocket, or beacon;
- unsafe label strings such as `</script><script>...` to prove the Base64/textContent strategy remains inert.

Add a static test that the verifier module does not import the producer builder module.

---

## 13. Release verifier behavior

`tests/verify_release_run.R` should remain checkout-bound and release-strict, but not tag-dependent.

Add:

- manifest revision 3 support;
- `allow_large_workload == false`;
- `online_preflight == false`;
- independent per-file producer source inventory verification;
- generic safe run-artifact resolver messages rather than Krona-specific wording;
- `exports.krona` / `exports.pavian` consistency;
- Pavian provenance/hash checks;
- per-sample invocation of the independent viewer correspondence verifier;
- all-unclassified/header-only resolution-sidecar support.

When Pavian is disabled:

- no Pavian path may be declared in `modules$kreport$outputs` or `exports$pavian$provenance_path`;
- preserved unowned files remain outside this rule.

When enabled:

- kreport must be requested and completed;
- Pavian provenance must be an owned kreport artifact;
- provenance sample IDs must equal manifest sample IDs exactly;
- per-sample paths must be unique, safe, run-relative, present, and hash-matched;
- correspondence verifier must pass for every sample;
- `render_html=false` requires no owned `.pavian.html` outputs;
- `render_html=true` requires one HTML per sample.

Do not reimplement the filename-sanitization algorithm in three languages merely to prove that a cosmetic basename matches a sample ID. Semantic identity is carried by exact `sample_id` inside provenance/JSON and is independently checked. The producer still retains its existing portable filename collision checks.

---

## 14. Release metadata checker and tag choreography

Add/expand `tests/check_release_metadata.py` as an independent stdlib-only release checker.

### 14.1 Branch/RC mode

Without `--require-tag`, validate:

- strict `VERSION` bytes;
- CFF version;
- CFF date-released;
- first non-Unreleased CHANGELOG version/date;
- README current-version anchor(s) without hard-coding “exactly four mentions.”

### 14.2 Exact-tag mode

With:

```text
--require-tag v0.4.7
```

also require:

- argument equals `v${VERSION}`;
- `refs/tags/v0.4.7` exists;
- tag object is annotated (`git cat-file -t` returns `tag`);
- `refs/tags/v0.4.7^{commit} == HEAD`;
- entire checkout is clean (`git status --porcelain --untracked-files=all` empty, ignoring Git-ignored files);
- annotated tagger timestamp normalized to UTC date equals CFF/CHANGELOG release date.

The current v0.4.6 tag is annotated, so this policy is consistent with existing release style.

### 14.3 Correct release sequence

The final choreography should be:

1. Merge/prepare clean v0.4.7 RC commit with `VERSION`, CFF, CHANGELOG, README already set to the intended release metadata.
2. Run full branch/RC CI and `check_release_metadata.py` **without** `--require-tag`.
3. On the release machine, verify clean HEAD and run full local release integration.
4. Create local annotated tag:

   ```bash
   git tag -a v0.4.7 -m "v0.4.7"
   ```

5. Run:

   ```bash
   python tests/check_release_metadata.py --require-tag v0.4.7
   ```

6. If it fails, delete the **local unpushed tag**, fix the commit/metadata, and repeat. Do not move a pushed tag.
7. Push the commit and tag.
8. Tag-push CI runs the full test matrix plus exact-tag metadata gate.
9. Only after tag CI is green, publish the GitHub Release.

This removes the circular “must already be tagged before creating the tag” condition from the attached draft.

---

## 15. File-level implementation plan

| File | Required v0.4.7 change |
|---|---|
| `analysis/utils/ncbi_taxonomy.py` | Outer lock cleanup scope; monotonic deadline; guaranteed candidate-FD cleanup; complete owner write; real cache-target lock identity. |
| `analysis/utils/atomic_io.R` | Canonical real taxonomy-cache lock/journal identity; journal schema v2; required transaction IDs for new journals; phase/path validation; rev3 prior-run recognition. |
| `tests/test_ncbi_taxonomy.py` | Interrupted acquisition, thread, relative/absolute, parent-symlink, final-cache-symlink, candidate-FD fault tests. |
| `tests/testthat/test-cross-runtime-lock.R` | Preserve both lock directions; add canonical alias coverage where portable. |
| `analysis/utils/preflight.R` | Emit producer `source_files` path/hash records plus aggregate digest; stable deterministic ordering. |
| `analysis/utils/manifest.R` | Add schema v2 revision 3; source-file object validation; `allow_large_workload`; `pavian`; `exports`; preserve revision 1/2 readers. |
| `analysis/00_run_pipeline.R` | Emit revision 3, source inventory, exports state; add Pavian↔kreport dependency check beside Krona; carry ownership/provenance. |
| `analysis/utils/config.R` | Add default/validation/CLI merge for `pavian.enabled` and `pavian.render_html`; no no-op renderer key. |
| `analysis/utils/cli.R` | Add `--pavian`; correct stale `--krona` help to describe builtin renderer. |
| `config.example.yml` | Add disabled-by-default Pavian block and clear official-vs-builtin wording. |
| `analysis/utils/kreport.R` | Add strict output-label validation; retain canonical U/R + DFS semantics. Do **not** export a parser for the independent verifier to reuse. |
| `analysis/07_kreport_pavian.R` | Always write resolution TSV; keep per-sample kreport map; invoke built-in viewer after resolution TSV exists; atomic Krona/Pavian provenance. |
| `analysis/utils/kraken_report_viewer.py` | New stdlib-only strict parser + resolution join + canonical JSON + self-contained CSP-hardened HTML builder. |
| `analysis/utils/verify_kraken_viewer_correspondence.py` | New independent parser/verifier; no import from builder. |
| `tests/test_kraken_report_viewer.py` | Parser, serialization, HTML/XSS/network, all-unclassified, determinism tests. |
| `tests/test_kraken_viewer_correspondence.py` | Full mutation matrix and independent-parser guard. |
| `tests/verify_release_run.R` | Revision 3, source inventory recomputation, release flags, generic artifact resolver, exports/provenance/viewer checks. |
| `tests/testthat/test-provenance.R` | Source inventory, exports state, moved run, no HTML when disabled, header-only resolution, atomic provenance. |
| `tests/testthat/test-kreport.R` | Label controls, U/R semantics, depth/rank invariants, deterministic order. |
| `tests/testthat/test-release-metadata.R` | Stop requiring exactly four README mentions; assert maintained anchors and date/version consistency. |
| `tests/check_release_metadata.py` | Branch metadata validation + optional exact annotated tag/clean checkout/date gate. |
| `tests/fixtures/...` | Add/reuse a single all-unclassified sample config to exercise header-only resolution + root=0 viewer. |
| `.github/workflows/ci.yml` | New unit tests; validate-only Pavian; combined Krona+Pavian integration; mixed zero-classified; all-unclassified smoke; tag-only metadata job. |
| `README.md` | Official Pavian upload boundary, built-in explorer boundary, outputs, offline/CSP behavior, TaxID/conflict semantics, semantic reproducibility statement. |
| `CHANGELOG.md` | v0.4.7 hardening, manifest rev3, source inventory, viewer, journal and release-gate changes; optionally clarify v0.4.6 chronology. |
| `CITATION.cff` / `VERSION` | Atomic 0.4.7 metadata update with actual release date. |

---

## 16. Rough implementation diffs

These blocks are intentionally **rough structural diffs**, not promised line-exact patches. They are specific enough for an implementation agent to follow while leaving room for local helper names and test fixtures.

### 16.1 `analysis/utils/ncbi_taxonomy.py` — cleanup and monotonic lock acquisition

```diff
@@
 @contextlib.contextmanager
 def acquire_cache_lock(cache_path, timeout=10.0, poll_interval=0.05):
-    lock_path = cache_path + ".lock"
+    cache_identity = os.path.realpath(os.path.abspath(cache_path))
+    lock_path = cache_identity + ".lock"
     os.makedirs(os.path.dirname(os.path.abspath(lock_path)), exist_ok=True)
-    lock_identity = os.path.normcase(os.path.realpath(lock_path))
-    with _ACTIVE_CACHE_LOCKS_GUARD:
-        if lock_identity in _ACTIVE_CACHE_LOCKS:
-            raise SystemExit(...)
-        _ACTIVE_CACHE_LOCKS.add(lock_identity)
-    deadline = time.time() + timeout
+    lock_identity = os.path.normcase(os.path.realpath(lock_path))
+    registered = False
     fd = None
-    acquired = False
-
-    while time.time() < deadline:
-        ...
-        time.sleep(poll_interval)
-
-    if not acquired:
-        with _ACTIVE_CACHE_LOCKS_GUARD:
-            _ACTIVE_CACHE_LOCKS.discard(lock_identity)
-        ...
-        raise SystemExit(...)
-
     try:
+        with _ACTIVE_CACHE_LOCKS_GUARD:
+            if lock_identity in _ACTIVE_CACHE_LOCKS:
+                raise SystemExit(...)
+            _ACTIVE_CACHE_LOCKS.add(lock_identity)
+            registered = True
+
+        deadline = time.monotonic() + timeout
+        while fd is None and time.monotonic() < deadline:
+            candidate_fd = None
+            try:
+                candidate_fd = os.open(lock_path, flags, 0o666)
+                try:
+                    _try_lock_fd(candidate_fd)  # lockf or msvcrt
+                except (OSError, IOError):
+                    pass
+                else:
+                    fd = candidate_fd
+                    candidate_fd = None
+                    break
+            except OSError:
+                pass
+            finally:
+                if candidate_fd is not None:
+                    try:
+                        os.close(candidate_fd)
+                    except OSError:
+                        pass
+
+            remaining = deadline - time.monotonic()
+            if remaining > 0:
+                time.sleep(min(poll_interval, remaining))
+
+        if fd is None:
+            owner_info = _read_lock_owner(lock_path)
+            raise SystemExit(
+                f"[taxonomy] ERROR: E_TAXONOMY_CACHE_BUSY: ... {owner_info}"
+            )
+
         os.lseek(fd, 0, os.SEEK_SET)
         os.ftruncate(fd, 0)
-        os.write(fd, owner_payload.encode("utf-8"))
+        _write_all(fd, owner_payload.encode("utf-8"))
         yield lock_path
     finally:
         if fd is not None:
             ... unlock ...
             ... close ...
-            with _ACTIVE_CACHE_LOCKS_GUARD:
-                _ACTIVE_CACHE_LOCKS.discard(lock_identity)
+        if registered:
+            with _ACTIVE_CACHE_LOCKS_GUARD:
+                _ACTIVE_CACHE_LOCKS.discard(lock_identity)
```

Test the **re-entrant failure before `registered=True`** carefully: its outer `finally` must not discard the first holder's set entry.

### 16.2 `analysis/utils/atomic_io.R` — real cache identity and journal schema

```diff
@@
 get_taxonomy_lock_path <- function(cache_path) {
-  paste0(canonicalize_root_path(cache_path), ".lock")
+  cache_identity <- normalizePath(cache_path, winslash = "/", mustWork = TRUE)
+  paste0(cache_identity, ".lock")
 }

 get_taxonomy_journal_path <- function(cache_path) {
-  paste0(normalizePath(cache_path, winslash = "/", mustWork = FALSE),
+  cache_identity <- normalizePath(cache_path, winslash = "/", mustWork = TRUE)
+  paste0(cache_identity,
          ".wf16s_transaction.json")
 }
@@
 write_publication_journal <- function(..., transaction_id = NULL) {
-  if (!is.null(transaction_id) && !valid_transaction_id(transaction_id)) ...
+  if (!valid_transaction_id(transaction_id)) {
+    stop("Publication journal requires a transaction_id.", call. = FALSE)
+  }
   payload <- list(
+    journal_schema_version = 2L,
     final_root = canonicalize_root_path(final_root),
     stage = ...,
     backup = ...,
+    had_prior = !is.null(backup),
     transaction_id = transaction_id,
     phase = phase,
@@
 recover_publication_journal <- function(final_root) {
   journal <- ...
-  if (is.null(journal) || is.null(journal$final_root)) ...
+  if (is.null(journal) || is.null(journal$final_root)) ...
+  if (identical(journal$journal_schema_version, 2L)) {
+    if (!valid_transaction_id(journal$transaction_id)) fail_closed(...)
+    if (!journal$phase %in% c("prepared", "prior_moved", "stage_published")) fail_closed(...)
+    assert_publication_temp_path(journal$stage, final_root, ".staging-")
+    if (isTRUE(journal$had_prior))
+      assert_publication_temp_path(journal$backup, final_root, ".previous-")
+  } else {
+    # Legacy v0.4.6 journal: conservative recovery only.
+    # Ambiguous identity => retain material and require manual recovery.
+  }
```

Apply the same explicit schema/transaction philosophy to taxonomy journals.

### 16.3 `analysis/utils/preflight.R` — source inventory records

```diff
@@
 source_provenance <- function(repo_root) {
   files <- maintained_source_files(repo_root)
   relative <- substring(files, nchar(normalizePath(repo_root, winslash = "/")) + 2L)
-  entries <- paste(relative, vapply(files, compute_file_hash, character(1)), sep = "\t")
-  digest_value <- digest::digest(paste(entries, collapse = "\n"), ...)
+  hashes <- vapply(files, compute_file_hash, character(1))
+  ord <- order(relative, method = "radix")
+  relative <- relative[ord]
+  hashes <- hashes[ord]
+  records <- lapply(seq_along(relative), function(i) list(
+    path = relative[[i]],
+    sha256 = hashes[[i]]
+  ))
+  canonical_lines <- paste(relative, hashes, sep = "\t")
+  digest_value <- digest::digest(
+    paste(canonical_lines, collapse = "\n"),
+    algo = "sha256", serialize = FALSE
+  )
@@
-  list(..., source_digest_sha256 = digest_value,
-       source_files = json_array(relative))
+  list(...,
+       source_digest_sha256 = digest_value,
+       source_files = json_array(records))
 }
```

Before hashing, validate that each tracked producer path is safe and is not a symlink in the stable release checkout.

### 16.4 `analysis/utils/manifest.R` — revision 3

```diff
@@
+validate_source_inventory <- function(records, path = "source_files") {
+  assert_manifest_array(records, path)
+  paths <- character(0)
+  for (i in seq_along(records)) {
+    rec <- records[[i]]
+    assert_manifest_scalar(rec$path, sprintf("%s[%d].path", path, i), "character")
+    assert_manifest_scalar(rec$sha256, sprintf("%s[%d].sha256", path, i), "character")
+    if (!is_safe_repo_relative_posix(rec$path)) manifest_fail(...)
+    if (!grepl("^[0-9a-f]{64}$", rec$sha256)) manifest_fail(...)
+    paths <- c(paths, rec$path)
+  }
+  if (anyDuplicated(paths) || anyDuplicated(tolower(paths))) manifest_fail(...)
+  if (!identical(paths, sort(paths, method = "radix"))) manifest_fail(...)
+}
+
+validate_manifest_v2_revision3 <- function(manifest, physical_root = NULL) {
+  # Keep revision 2 frozen; validate revision-3 root and revision-2-compatible
+  # core fields using shared leaf helpers rather than rewriting schema_revision.
+  ...
+  if (!identical(manifest$schema_revision, 3L)) manifest_fail(...)
+  validate_source_inventory(manifest$source_files)
+  for (field in c("validate_only", "keep_going", "overwrite",
+                  "allow_unlocked", "allow_dirty", "allow_large_workload",
+                  "online_preflight", "refresh_taxonomy", "krona", "pavian")) {
+    assert one logical ...
+  }
+  validate_exports(manifest$exports, manifest)
+  ... existing revision-2 artifact/census/environment checks ...
+}
@@
 validate_manifest_v2 <- function(manifest, physical_root = NULL) {
   rev <- manifest$schema_revision
-  if (identical(rev, 2L)) {
+  if (identical(rev, 3L)) {
+    validate_manifest_v2_revision3(manifest, physical_root)
+  } else if (identical(rev, 2L)) {
     validate_manifest_v2_revision2(manifest, physical_root)
```

Also extend `manifest_is_valid_run()` and prior-output recognition to revision 3.

### 16.5 `analysis/utils/config.R`, `analysis/utils/cli.R`, `config.example.yml`

```diff
@@ get_default_config()
     krona = list(...),
+    pavian = list(
+      enabled = FALSE,
+      render_html = TRUE
+    ),
@@ validate_config()
+  if (!is.list(cfg$pavian)) stop("'pavian' must be a configuration mapping.")
+  assert one logical cfg$pavian$enabled
+  assert one logical cfg$pavian$render_html
@@ load_config()
   cli_krona <- ...
+  cli_pavian <- isTRUE(cli_opts$pavian) || isTRUE(cli_opts[["pavian"]])
@@
   if (cli_krona) cfg$krona$enabled <- TRUE
+  if (cli_pavian) cfg$pavian$enabled <- TRUE
@@ cfg$cli
     krona = isTRUE(cfg$krona$enabled),
+    pavian = isTRUE(cfg$pavian$enabled),
```

```diff
@@ analysis/utils/cli.R
     make_option("--krona",
-      help = "Enable Krona-compatible TSV output; optional HTML rendering requires KronaTools ktImportText"),
+      help = "Enable Krona export; builtin offline HTML is the default renderer"),
+    make_option("--pavian",
+      action = "store_true", default = FALSE, dest = "pavian",
+      help = "Enable official-Pavian-compatible kreport integration plus the builtin offline Kraken-report explorer"),
```

```diff
@@ config.example.yml
 krona:
   ...
+
+pavian:
+  # .kreport remains the official Pavian interoperability artifact.
+  # The optional HTML is an original offline Kraken-report explorer,
+  # not the upstream Pavian Shiny application.
+  enabled: false
+  render_html: true
```

### 16.6 `analysis/00_run_pipeline.R` — dependency and manifest fields

```diff
@@
 if (isTRUE(cfg$krona$enabled) && !("kreport" %in% requested_modules)) {
   fatal(... "Krona export requires the 'kreport' module.")
 }
+if (isTRUE(cfg$pavian$enabled) && !("kreport" %in% requested_modules)) {
+  fatal(... "Pavian export requires the 'kreport' module.")
+}
@@ manifest <- list(
   git_commit = source_info$git_commit,
   git_dirty = source_info$git_dirty,
   source_digest_sha256 = source_info$source_digest_sha256,
-  schema_version = 2L, schema_revision = 2L,
+  source_files = source_info$source_files,
+  schema_version = 2L, schema_revision = 3L,
@@
   cli = utils::modifyList(cfg$cli, list(modules = json_array(requested_modules))),
+  exports = list(
+    krona = list(
+      enabled = isTRUE(cfg$krona$enabled),
+      render_html = isTRUE(cfg$krona$enabled) && isTRUE(cfg$krona$render_html),
+      provenance_path = if (isTRUE(cfg$krona$enabled))
+        "07_Kreport/krona/krona_provenance.json" else NULL
+    ),
+    pavian = list(
+      enabled = isTRUE(cfg$pavian$enabled),
+      render_html = isTRUE(cfg$pavian$enabled) && isTRUE(cfg$pavian$render_html),
+      provenance_path = if (isTRUE(cfg$pavian$enabled))
+        "07_Kreport/pavian/pavian_provenance.json" else NULL,
+      integration = "official_pavian_upload_plus_builtin_kraken_report_explorer",
+      official_pavian_compatibility = "kraken_report_input_contract_only"
+    )
+  ),
```

If an export provenance file is expected but the kreport module did not actually produce it, manifest construction/validation must fail rather than publishing a contradictory `exports` record.

### 16.7 `analysis/07_kreport_pavian.R` — two-pass kreport/Pavian flow

```diff
@@ run_kreport <- function(context) {
   ...
+  pavian_cfg <- cfg$pavian %||% list(enabled = FALSE, render_html = TRUE)
+  pavian_enabled <- isTRUE(pavian_cfg$enabled)
+  pavian_render_html <- pavian_enabled && isTRUE(pavian_cfg$render_html)
+  pavian_dir <- file.path(kreport_dir, "pavian")
+  kreport_by_sample <- list()
@@ for (s in samples) {
     out_file <- file.path(kreport_dir, sprintf("%s.kreport", sanitize_filename(s)))
     writeLines(kreport_lines, out_file)
     all_outputs <- c(all_outputs, out_file)
+    kreport_by_sample[[s]] <- out_file
@@
   }

-  if (length(resolution_rows) > 0) {
-    res_df <- do.call(rbind, resolution_rows)
-    write.table(res_df, res_summary_file, ...)
-    all_outputs <- c(all_outputs, res_summary_file)
-  }
+  res_df <- if (length(resolution_rows)) {
+    do.call(rbind, resolution_rows)
+  } else {
+    data.frame(
+      SampleID = character(), Depth = integer(), RankCode = character(),
+      NodeName = character(), TaxonPath = character(), TaxID = character(),
+      Status = character(), ResolutionSource = character(),
+      stringsAsFactors = FALSE
+    )
+  }
+  write.table(res_df, res_summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
+  all_outputs <- c(all_outputs, res_summary_file)
+
+  if (pavian_enabled) {
+    dir.create(pavian_dir, recursive = TRUE, showWarnings = FALSE)
+    pavian_records <- list()
+    for (s in samples) {
+      stem <- sanitize_filename(s)
+      json_out <- file.path(pavian_dir, sprintf("%s.pavian.json", stem))
+      html_out <- if (pavian_render_html)
+        file.path(pavian_dir, sprintf("%s.pavian.html", stem)) else NULL
+      render_builtin_kraken_report_explorer(
+        python_cmd = python_cmd,
+        builder_path = file.path(cfg$pipeline_root, "analysis", "utils", "kraken_report_viewer.py"),
+        kreport_path = kreport_by_sample[[s]],
+        resolution_tsv = res_summary_file,
+        sample_id = s,
+        expected_total = sum(count_matrix[, s]),
+        json_out = json_out,
+        html_out = html_out
+      )
+      ... collect paths/hashes/count-status summaries ...
+      all_outputs <- c(all_outputs, json_out, html_out %||% character())
+    }
+    atomic_write_json(pavian_provenance, pavian_provenance_file)
+    all_outputs <- c(all_outputs, pavian_provenance_file)
+  }
@@ Krona provenance
-    jsonlite::write_json(krona_provenance, krona_provenance_file, ...)
+    atomic_write_json(krona_provenance, krona_provenance_file)
```

The important sequencing change is: **write the complete resolution sidecar first, then render Pavian JSON/HTML**.

### 16.8 New `analysis/utils/kraken_report_viewer.py`

```diff
+++ analysis/utils/kraken_report_viewer.py
+#!/usr/bin/env python3
+"""Build deterministic offline Kraken-report explorer artifacts.
+
+The .kreport remains the official upstream Pavian interoperability artifact.
+This file implements an original static viewer and does not vendor Pavian UI code.
+"""
+
+from __future__ import annotations
+import argparse, base64, csv, hashlib, html, json, os, re, tempfile
+from pathlib import Path
+
+MAX_SAFE_INTEGER = 9007199254740991
+RANK_BY_DEPTH = {1:"D", 2:"K", 3:"P", 4:"C", 5:"O", 6:"F", 7:"G", 8:"S"}
+
+def parse_kreport(path: Path) -> dict:
+    # strict six-column parsing
+    # U then R special rows
+    # exact indentation/depth/rank stack
+    # direct/clade arithmetic
+    # deterministic DFS/sibling-order validation
+    ...
+
+def load_resolution_rows(path: Path, sample_id: str) -> dict[str, dict]:
+    # exact header; unique SampleID+TaxonPath; safe fields
+    ...
+
+def merge_resolution(report: dict, resolution: dict) -> list[dict]:
+    # exact path/depth/rank/name/TaxID cross-check
+    # preserve Resolved/Unresolved/Conflicted + ResolutionSource
+    ...
+
+def canonical_json_bytes(payload: dict) -> bytes:
+    text = json.dumps(
+        payload, ensure_ascii=False, sort_keys=True,
+        separators=(",", ":"), allow_nan=False
+    )
+    return text.encode("utf-8") + b"\n"
+
+def build_html(payload_bytes: bytes, sample_id: str) -> bytes:
+    payload_b64 = base64.b64encode(payload_bytes).decode("ascii")
+    safe_title = html.escape(sample_id, quote=True)
+    # Inline CSS/JS only; CSP connect-src 'none'.
+    # DOM labels must use textContent, never innerHTML.
+    ...
+
+def atomic_write_bytes(path: Path, data: bytes) -> None:
+    if path.exists():
+        raise ValueError(f"refusing to overwrite existing builder output: {path}")
+    ... tempfile in sibling dir, fsync optional, os.replace ...
+
+def main() -> int:
+    ...
```

Unit tests should treat any parser normalization that silently fixes malformed report input as a bug. The builder is a validator/renderer, not a repair tool.

### 16.9 New independent verifier

```diff
+++ analysis/utils/verify_kraken_viewer_correspondence.py
+#!/usr/bin/env python3
+"""Independently verify kreport -> resolution sidecar -> JSON -> HTML identity."""
+
+# IMPORTANT: no import from kraken_report_viewer.py.
+
+def parse_report_independently(...):
+    ...
+
+def parse_resolution_independently(...):
+    ...
+
+def extract_embedded_payload_bytes(html_text: str) -> bytes:
+    # locate exactly one data element, validate base64, decode strictly
+    ...
+
+def reject_network_capabilities(html_text: str) -> None:
+    # external src/href where resource-loading applies
+    # fetch/XMLHttpRequest/WebSocket/EventSource/sendBeacon
+    # missing/weak CSP
+    ...
+
+def verify(...):
+    report = parse_report_independently(...)
+    resolution = parse_resolution_independently(...)
+    payload_bytes = Path(json_path).read_bytes()
+    payload = json.loads(payload_bytes.decode("utf-8"))
+    ... exact ordered semantic comparison ...
+    if html_path is not None:
+        embedded = extract_embedded_payload_bytes(Path(html_path).read_text("utf-8"))
+        if embedded != payload_bytes:
+            raise ValueError("embedded payload differs byte-for-byte from JSON")
+        reject_network_capabilities(...)
```

### 16.10 `tests/verify_release_run.R` — source/export checks

```diff
@@
-stopifnot(identical(manifest$schema_revision, 2L))
+stopifnot(identical(manifest$schema_revision, 3L))
 validate_manifest_v2(manifest, physical_root = root)
@@
 stopifnot(identical(manifest$cli$allow_dirty, FALSE))
 stopifnot(identical(manifest$cli$allow_unlocked, FALSE))
+stopifnot(identical(manifest$cli$allow_large_workload, FALSE))
+stopifnot(identical(manifest$cli$online_preflight, FALSE))
@@
-resolve_run_artifact <- function(...) {
-  ... "Krona provenance ..." ...
+resolve_run_artifact <- function(...) {
+  ... generic "run provenance artifact path" messages ...
 }
+
+# Independently derive expected producer-source files from git and scope rules.
+expected_source_records <- recompute_release_source_inventory(repo_root)
+assert_source_inventory_matches(manifest$source_files, expected_source_records)
+stopifnot(identical(
+  manifest$source_digest_sha256,
+  canonical_source_digest(expected_source_records)
+))
+
+verify_export_record(manifest$exports$krona, ...)
+verify_export_record(manifest$exports$pavian, ...)
+
+if (isTRUE(manifest$exports$pavian$enabled)) {
+  pavian_prov <- read/validate ...
+  for (record in pavian_prov$samples) {
+    system2(python, c(
+      viewer_verify_script,
+      "--kreport", resolve_run_artifact(...),
+      "--resolution-tsv", file.path(root, "07_Kreport/taxonomy_resolution.tsv"),
+      "--json", resolve_run_artifact(...),
+      if (render_html) c("--html", resolve_run_artifact(...)) else character(),
+      "--sample-id", record$sample_id,
+      "--expected-total", as.character(record$total_reads)
+    ))
+    stop on nonzero
+  }
+}
```

### 16.11 `tests/check_release_metadata.py` — exact annotated tag

```diff
+++ tests/check_release_metadata.py
+import argparse, datetime as dt, pathlib, re, subprocess
+
+def git(*args):
+    return subprocess.check_output(["git", *args], text=True).strip()
+
+version = read_strict_version(...)
+cff_version, cff_date = read_top_level_cff_fields(...)
+changelog_version, changelog_date = read_first_release_heading(...)
+assert version == cff_version == changelog_version
+assert cff_date == changelog_date
+assert_readme_current_version_anchor(version)
+
+if args.require_tag:
+    expected = f"v{version}"
+    assert args.require_tag == expected
+    assert git("cat-file", "-t", f"refs/tags/{expected}") == "tag"
+    assert git("rev-parse", f"refs/tags/{expected}^{{commit}}") == git("rev-parse", "HEAD")
+    assert git("status", "--porcelain", "--untracked-files=all") == ""
+    tagger_iso = git("for-each-ref", "--format=%(taggerdate:iso-strict)", f"refs/tags/{expected}")
+    assert parse(tagger_iso).astimezone(dt.timezone.utc).date().isoformat() == cff_date
```

Do not parse arbitrary YAML generally; for CFF release metadata, a narrow strict top-level field parser is sufficient and keeps the checker independent of the R environment.

### 16.12 `.github/workflows/ci.yml`

```diff
@@ unit tests
   python -m unittest -v tests/test_krona_builder.py
   python -m unittest -v tests/test_krona_correspondence.py
+  python -m unittest -v tests/test_kraken_report_viewer.py
+  python -m unittest -v tests/test_kraken_viewer_correspondence.py
@@ validate-only
+  test ! -e "ci pavian validation output"
+  Rscript analysis/00_run_pipeline.R \
+    --config config.yml \
+    --output-dir "ci pavian validation output" \
+    --validate-only --pavian
+  test ! -e "ci pavian validation output"
@@ full release integration
-  Rscript analysis/00_run_pipeline.R ... --krona
+  Rscript analysis/00_run_pipeline.R ... --krona --pavian
   Rscript tests/verify_release_run.R "ci integration output"
@@ mixed zero-classified cohort
-  ... --krona
+  ... --krona --pavian
   Rscript tests/verify_release_run.R ...
+
+  # Add a single/all-unclassified config selecting S_zero only.
+  Rscript analysis/00_run_pipeline.R \
+    --config tests/fixtures/synthetic_zero_classified_single/config.yml \
+    --modules qc,composition,kreport \
+    --output-dir "ci all unclassified output" \
+    --pavian
+  Rscript tests/verify_release_run.R "ci all unclassified output"
+
+  python tests/check_release_metadata.py
+
+tag-metadata:
+  if: startsWith(github.ref, 'refs/tags/v')
+  needs: [test, faprotax]
+  ... checkout fetch-depth: 0 ...
+  run: python tests/check_release_metadata.py --require-tag "${GITHUB_REF_NAME}"
```

The main OS matrix already gives Windows/Linux coverage. Do not create a combinatorial explosion of separate jobs for every flag pairing when unit mutation tests can cover the fine-grained cases.

---

## 17. CI and acceptance matrix

### 17.1 Existing hardening that must remain green

On Ubuntu and Windows, R 4.5.3, Python 3.12:

- all Python unit tests and `compileall`;
- all R `testthat` tests;
- default `--validate-only` mutation-free;
- Krona `--validate-only` mutation-free;
- full release integration;
- mixed zero-classified cohort without alpha;
- phylogenetic alpha fixture;
- FAPROTAX and non-FAPROTAX paths;
- committed and commit-range whitespace gates;
- moved-run verifier behavior;
- physical file census and ownership gates.

### 17.2 New v0.4.7 gates

- interrupted taxonomy acquisition never leaves `_ACTIVE_CACHE_LOCKS` poisoned;
- final-component symlink aliases converge on one taxonomy lock identity;
- source inventory has exact safe path/hash records and independently recomputed aggregate digest;
- mutated source inventory path/hash/order fails closed;
- release verifier rejects `allow_large_workload=true`;
- release verifier rejects `online_preflight=true`;
- revision-2 historical fixture still validates as revision 2;
- revision-3 manifest cannot be accepted as revision 2;
- `--pavian` without `kreport` fails before output mutation;
- `--validate-only --pavian` creates no output;
- full `--krona --pavian` run has disjoint owned export sets;
- mixed cohort with classified and zero-classified samples renders/verifies;
- all-unclassified sample emits U/R report + header-only resolution TSV + valid empty-node JSON/HTML;
- unresolved TaxID 0 remains visibly unresolved;
- assignment-conflict node remains visibly conflicted even when winning TaxID is nonzero;
- JSON/HTML exact payload byte identity;
- malicious taxonomy label does not become executable HTML/JS;
- generated HTML has restrictive CSP and no remote resource/network APIs;
- repeated render from same input has identical JSON/HTML hashes on each OS;
- tag checker rejects lightweight tag, wrong tag, wrong target, dirty checkout, date mismatch, CFF mismatch, and CHANGELOG mismatch;
- exact annotated local tag passes before push; tag-push CI independently passes after push.

---

## 18. Recommended implementation order

1. **Fix Python lock cleanup and canonical cache identity first.** Add interruption/thread/symlink tests before touching viewer code.
2. **Harden R taxonomy lock/journal identity** to match Python and introduce journal schema v2.
3. **Add source path/hash inventory** and independent release-verifier recomputation.
4. **Add manifest revision 3** with `allow_large_workload`, `pavian`, `source_files`, and `exports`; freeze rev2.
5. **Add release metadata checker** and remove README “exactly four mentions” brittleness.
6. **Make `taxonomy_resolution.tsv` unconditional** and add strict kreport label validation.
7. **Lock down the written kreport specification** with producer tests; do not yet share a parser with the new verifier.
8. **Implement `kraken_report_viewer.py`**: independent strict parser, resolution join, deterministic JSON, CSP-hardened Base64 HTML.
9. **Integrate the viewer in a second pass of the kreport module** after the resolution sidecar is written.
10. **Implement the independent correspondence verifier** and mutation matrix.
11. **Extend release verifier and provenance ownership** for Pavian outputs.
12. **Add CI integrations**: validate-only Pavian, combined Krona/Pavian, mixed zero-classified, and all-unclassified.
13. **Update README/CHANGELOG/CFF/VERSION** and state semantic reproducibility explicitly.
14. **Run pre-tag RC gates** from clean HEAD.
15. **Create local annotated `v0.4.7` tag**, run exact-tag gate locally, then push.
16. **Publish GitHub Release only after tag CI passes.**

This order deliberately keeps release-hardening work ahead of the visualization feature, so a viewer bug cannot obscure a provenance or transaction regression.

---

## 19. Release gate commands

Illustrative local pre-tag gate:

```bash
Rscript analysis/install_packages.R --restore
Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
python -m compileall -q analysis tests
Rscript tests/testthat.R
python -m unittest -v
python tests/check_release_metadata.py
python tests/check_committed_whitespace.py
git diff --check

git status --porcelain --untracked-files=all
```

Run at least one **full/default maintained-module release integration** with both exports, not only the targeted no-alpha fixture:

```bash
Rscript analysis/00_run_pipeline.R \
  --config config.yml \
  --output-dir "release-full-output" \
  --krona --pavian
Rscript tests/verify_release_run.R "release-full-output"
```

Then targeted mixed/no-alpha and all-unclassified integrations:

```bash
Rscript analysis/00_run_pipeline.R \
  --config tests/fixtures/synthetic_zero_classified_cohort/config.yml \
  --modules qc,composition,kreport \
  --output-dir "release-zero-mixed-output" \
  --krona --pavian
Rscript tests/verify_release_run.R "release-zero-mixed-output"

Rscript analysis/00_run_pipeline.R \
  --config tests/fixtures/synthetic_zero_classified_single/config.yml \
  --modules qc,composition,kreport \
  --output-dir "release-all-unclassified-output" \
  --pavian
Rscript tests/verify_release_run.R "release-all-unclassified-output"
```

Only after the pre-tag gate is green:

```bash
git tag -a v0.4.7 -m "v0.4.7"
python tests/check_release_metadata.py --require-tag v0.4.7
```

If that succeeds, push commit/tag and wait for tag CI before creating the GitHub Release.

The release verifier must additionally establish from each release run:

```text
manifest schema_version == 2
manifest schema_revision == 3
manifest pipeline_version == VERSION
manifest git_commit == HEAD
manifest git_dirty == false
manifest environment.locked == true
manifest cli.allow_dirty == false
manifest cli.allow_unlocked == false
manifest cli.allow_large_workload == false
manifest cli.refresh_taxonomy == false
manifest cli.online_preflight == false
manifest source_files path/hash inventory == independently recomputed producer inventory
manifest source_digest_sha256 == canonical digest(source_files)
taxonomy network mode == cache_only for stable fixtures
all owned artifact hashes and physical census reconcile
all enabled export provenance records reconcile
all Pavian viewer correspondence checks pass
```

---

## 20. Documentation wording that should be explicit

README should state, in substance:

- `.kreport` is the artifact intended for official Pavian upload/use.
- `--pavian` additionally creates an original offline Kraken-report explorer; it is **not** the upstream Pavian Shiny application.
- unresolved TaxID 0 and assignment-conflict states are preserved rather than silently “fixed.”
- built-in viewer rendering is offline and performs no taxonomy lookup.
- the stable release promises **semantic reproducibility** under the recorded source/environment/input contracts; it does not promise all plots are byte-identical across operating systems.
- JSON/HTML viewer payloads are deterministic under identical canonical inputs.

CLI `--krona` help must no longer imply that HTML requires external KronaTools when builtin rendering is the default.

The release-metadata test should verify intended README anchors rather than a magic count of version mentions.

---

## Final disposition

The attached plan was directionally strong, especially on the September lock fixes, composition-without-alpha repair, manifest/source provenance gap, and the licensing boundary around upstream Pavian. After stress-testing it against the current code, the implementation should be adjusted in four important ways:

1. **Treat lock identity and cleanup as a complete cross-runtime transaction**, including cache-file symlink aliases and interruption paths.
2. **Separate producer-source provenance from release-harness provenance**, using per-file runtime hashes plus the exact Git commit/tag for the full release tree.
3. **Make Pavian status semantics depend on both `.kreport` and the always-present `taxonomy_resolution.tsv` sidecar**, because conflict state is not representable in the report alone.
4. **Use a two-phase local-tag/push-tag release choreography**, rather than an impossible “exact tag before tag creation” gate.

With these changes, the v0.4.7 patch becomes internally coherent enough to hand directly to an implementation agent. The safest feature scope is still:

- canonical independently verified six-column `.kreport` for official Pavian interoperability;
- original deterministic offline Kraken-report explorer as an optional kreport subextension;
- exact JSON/HTML/provenance ownership and correspondence checks;
- no vendoring or relabeling of the upstream Pavian Shiny UI;
- manifest revision 3 with auditable producer-source inventory;
- strict pre-tag and tag-push release identity gates.

A future optional launcher for a user-installed upstream Pavian package/container can be planned separately; it should remain outside the static artifact renderer and outside the v0.4.7 critical path.
