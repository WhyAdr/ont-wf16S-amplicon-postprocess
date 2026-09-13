# v0.4.8 Taxonomy Sankey Renderer — Refined Production Plan

**Status:** implementation-ready **only after v0.4.7 is green, tagged, and released**
**Repository:** `WhyAdr/ont-wf16S-amplicon-postprocess`
**Audited `main`:** `11790b4ce41f093e3d2379489a2f199ff5cd6f2a`
**Current pipeline version:** `0.4.7`

This document **supersedes** the earlier `wf16s-postprocess-v0.4.8-sankey-renderer-plan.md`.

The overall architecture remains valid: keep `.kreport` as the official upstream Pavian interoperability artifact, keep the v0.4.7 `.pavian.json/.html` explorer, and add an original offline `builtin_taxonomy_sankey` renderer. This revision closes two conservation flaws, removes a verifier trust-loop, makes IDs/order exact, extends provenance/manifest contracts, defines overwrite behavior, and adds transient-view/stress tests.

---

## 0. Mandatory prerequisite: repair and release v0.4.7 first

The Sankey feature must not start from the current red baseline.

At the release-gate refresh:

- remote `main` is `11790b4ce41f093e3d2379489a2f199ff5cd6f2a`;
- `VERSION` is `0.4.7`;
- Ubuntu and Windows primary jobs are green;
- both FAPROTAX jobs are green;
- Pavian/Krona integration stages execute and pass;
- annotated `v0.4.7` points to the audited `main` commit and its tag workflow is green;
- publication of the GitHub Release object remains an explicit release-admin step.

The earlier red baseline was caused by the test-only symlink assertion described
below. The release metadata check was also moved before mutating integration
stages so its clean-checkout assertion is meaningful.

The production lock implementation is already conceptually correct:

```python
cache_identity = os.path.realpath(os.path.abspath(os.fspath(cache_path)))
lock_path = cache_identity + ".lock"
```

The failing test incorrectly resolves `alias + ".lock"` after appending the suffix. A final-component cache symlink cannot be resolved that way.

### Required test fix

```diff
 def test_final_cache_symlink_shares_lock_identity(self):
     alias = self.work / "cache-alias.json"
     try:
         os.symlink(self.cache, alias)
     except (OSError, NotImplementedError) as exc:
         self.skipTest(f"cache symlink unavailable: {exc}")

+    expected_lock = pathlib.Path(
+        os.path.realpath(os.path.abspath(str(alias))) + ".lock"
+    )
+    target_lock = pathlib.Path(
+        os.path.realpath(os.path.abspath(str(self.cache))) + ".lock"
+    )
     self.assertEqual(
-        os.path.normcase(os.path.realpath(str(alias) + ".lock")),
-        os.path.normcase(os.path.realpath(str(self.cache) + ".lock")),
+        os.path.normcase(str(expected_lock)),
+        os.path.normcase(str(target_lock)),
     )

     with taxonomy.acquire_cache_lock(str(alias), timeout=0.2):
+        self.assertTrue(target_lock.is_file())
+        self.assertFalse(pathlib.Path(str(alias) + ".lock").exists())
         with self.assertRaises(SystemExit):
             with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
                 pass
```

The behavioral contention test is the authoritative contract. Do **not** change the production lock-path algorithm merely to satisfy the old assertion.

### v0.4.7 completion gate

Before starting v0.4.8:

```text
full R suite green
full Python suite green
Ubuntu primary CI green
Windows primary CI green
both FAPROTAX jobs green
Pavian/Krona integration stages actually execute and pass
annotated local v0.4.7 tag created
exact-tag metadata check passes
tag pushed immutably
tag CI green
GitHub Release published
```

The final two lines are separate: a green annotated-tag workflow proves the
repository gate, while the GitHub Release object still requires authenticated
release publication. Do not mark v0.4.8 implementation-ready until both are
recorded.

---

## 1. Feature boundary

Add an original renderer:

```text
builtin_taxonomy_sankey
```

Do **not** vendor/copy Pavian or `fbreitwieser/sankeyD3` code. Both are GPL-family projects. The new renderer may emulate the useful behavior—rank columns, abundance-proportional flow, top-N display, highlighting, count/percentage labels—but must be independently implemented.

Runtime remains:

```text
Python 3.12 standard library
vanilla JavaScript
SVG
```

No runtime D3, npm, Node, CDN, Shiny, or new R plotting dependency.

Node is allowed only as a pinned **CI test runtime** for DOM-free client model tests.

---

## 2. Output architecture

Keep:

```text
07_Kreport/
  <sample>.kreport
  taxonomy_resolution.tsv

  pavian/
    <sample>.pavian.json
    <sample>.pavian.html
    pavian_provenance.json
```

Add:

```text
07_Kreport/
  pavian/
    sankey/
      <sample>.sankey.json
      <sample>.sankey.html
```

Semantics:

- `.kreport` = official upstream Pavian input artifact;
- `.pavian.json` = canonical verified normalized source payload;
- `.pavian.html` = built-in hierarchy/table explorer;
- `.sankey.json` = canonical verified taxonomy-flow model;
- `.sankey.html` = original offline Sankey viewer.

Do not call the built-in HTML "Pavian Sankey". Prefer "Taxonomy Sankey" or "Pavian-inspired taxonomy flow".

---

## 3. The corrected conservation model

The earlier plan allowed residual mass to stop at an intermediate rank. That contradicts the stronger invariant that every rank transition and the final column represent all classified reads.

This revision chooses **persistent residual lanes**.

Let:

```text
C = classified reads
selected ranks = R0, R1, ..., Rm-1
```

Every displayed selected-rank column must satisfy:

```text
sum(node.value) == C
```

Every transition must satisfy:

```text
sum(link.value crossing transition) == C
```

and therefore:

```text
rightmost_flow == C
```

### Why residual lanes must persist

Example:

```text
C = 100
Domain = 100

visible Kingdom = 60
hidden Kingdom  = 40
visible Phylum under visible Kingdom = 60
```

Incorrect model:

```text
D -> visible K = 60
D -> Other K   = 40
visible K -> P = 60
```

The Phylum column totals only 60.

Correct model:

```text
D -> visible K      = 60
D -> Other K        = 40

visible K -> P      = 60
Other K -> carry    = 40
```

Phylum column:

```text
60 + 40 = 100
```

The carry is explicitly synthetic/non-biological. It means "this abundance was no longer expanded in the visible taxonomy at Kingdom", not "this is a Phylum assignment".

---

## 4. Explicit classified-entry node

Every view begins from:

```text
synthetic:entry:classified
```

with:

```text
entry.value = C
```

This fixes two cases:

1. top-N pruning at the **first** selected rank;
2. selected ranks beginning below Domain, e.g. `P,G`.

For first selected rank `R0`:

```text
ALL_R0     = all biological nodes at R0
VISIBLE_R0 = retained biological nodes at R0
HIDDEN_R0  = ALL_R0 - VISIBLE_R0

all_mass_R0     = sum(clade(ALL_R0))
hidden_mass_R0  = sum(clade(HIDDEN_R0))
assigned_above  = C - all_mass_R0
```

Required:

```text
assigned_above >= 0
```

Entry links:

```text
ENTRY -> each visible R0 node        value = node.clade
ENTRY -> Other R0                    value = hidden_mass_R0
ENTRY -> Assigned above R0           value = assigned_above
```

Omit zero-valued residual nodes/links.

Exact:

```text
sum(ENTRY outgoing links) == C
```

---

## 5. Biological transition semantics

For retained real source `S` at `Ri`, targeting `Ri+1`:

```text
ALL_TARGETS(S)
    = all biological descendants of S at target rank

VISIBLE_TARGETS(S)
    = retained descendants at target rank

HIDDEN_TARGETS(S)
    = ALL_TARGETS(S) - VISIBLE_TARGETS(S)
```

Then:

```text
visible_mass =
    sum(target.clade for VISIBLE_TARGETS)

hidden_mass =
    sum(target.clade for HIDDEN_TARGETS)

all_target_mass =
    visible_mass + hidden_mass

assigned_above_next =
    S.clade - all_target_mass
```

Required:

```text
assigned_above_next >= 0
```

Emit:

```text
S -> visible target(s)
S -> source-specific Other <target rank>
S -> source-specific Assigned above <target rank>
```

Exact per-source invariant:

```text
S.clade
==
sum(all outgoing link values)
```

For skipped ranks, e.g. `P -> G`, `ALL_TARGETS(S)` means every Genus descendant under the Phylum regardless of omitted Class/Order/Family columns.

---

## 6. Persistent residual-lane semantics

Residual subtypes:

```text
other_hidden
assigned_above
```

When a residual originates at target rank `Rj`, create a residual node there. For every later selected rank, create a **carry node** with exactly the same lane identity and value.

Required for every residual/carry node before the final column:

```text
outdegree == 1
outgoing.value == node.value
target.value == node.value
target.lane_id == node.lane_id
```

A residual lane never re-enters a biological node.

### Why this is safe with top-N

After initial top-N selection, apply **ancestral closure**. Any retained lower-rank biological node forces its ancestor to be visible at every earlier selected rank. Thus a branch hidden into `Other Family` can never later "reappear" as a visible Genus/Species. Persistent residual lanes remain semantically closed.

---

## 7. Top-N + ancestral closure

For each selected rank independently:

```text
sort by:
  1. descending integer clade
  2. canonical TaxonPath UTF-8 byte order
retain first max_taxa_per_rank
```

Mark:

```text
selection_reason = "top_n"
```

Then for every retained lower-rank node, retain its ancestor at every earlier selected rank.

Ancestors added only for closure:

```text
selection_reason = "ancestor_closure"
```

A previous rank may therefore exceed N. That is intentional and preferable to breaking lineage connectivity.

Never use locale-aware ordering in canonical selection.

---

## 8. Exact IDs

### Biological node

```python
def sha256_utf8(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()

def taxon_node_id(rank_code, path):
    return f"taxon:{rank_code}:{sha256_utf8(path)}"
```

Full path remains stored in the node record.

### Entry

```text
synthetic:entry:classified
```

### Residual lane

Canonical lane key:

```text
subtype NUL origin_source_id NUL origin_target_rank
```

```python
def residual_lane_id(subtype, origin_source_id, origin_target_rank):
    raw = "\0".join([subtype, origin_source_id, origin_target_rank])
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()
```

### Residual node by column

```python
def residual_node_id(lane_id, column_rank):
    return f"synthetic:residual:{lane_id}:{column_rank}"
```

Origin and carry nodes share the same `lane_id`.

### Link

```python
def link_id(source_id, target_id, kind):
    raw = "\0".join([source_id, target_id, kind])
    return "link:" + hashlib.sha256(raw.encode("utf-8")).hexdigest()
```

Do **not** include value in link identity. Mutating a value should mutate the same edge semantically rather than creating a new identity.

---

## 9. Synthetic-node audit metadata

`other_hidden` origin node records:

```json
{
  "kind": "residual",
  "subtype": "other_hidden",
  "biological": false,
  "carried": false,
  "lane_id": "...",
  "origin_source_id": "...",
  "origin_target_rank": "K",
  "column_rank": "K",
  "value": 40,
  "member_count": 8,
  "member_paths_sha256": "..."
}
```

Digest membership as:

```python
blob = "\n".join(sorted(hidden_paths, key=lambda x: x.encode("utf-8")))
member_paths_sha256 = sha256_utf8(blob)
```

`assigned_above` uses:

```text
member_count = null
member_paths_sha256 = null
```

Carry node copies lane metadata/value and sets:

```text
kind = residual_carry
carried = true
column_rank = current selected rank
```

Synthetic nodes have **no TaxID**. Never use TaxID `0` for presentation-generated nodes; `0` is reserved for unresolved biological taxonomy.

---

## 10. Exact canonical ordering

Serialization order and visual order are separate.

### Source biological index

`source_order_index` is the zero-based index in canonical `.pavian.json` `nodes`.

### Canonical node serialization

Columns:

```text
ENTRY = -1
R0 = 0
R1 = 1
...
```

Kind order:

```text
entry          = 0
taxon          = 1
residual       = 2
residual_carry = 3
```

Sort nodes by:

```text
column_index
kind_order
canonical_anchor
id
```

where:

```text
taxon canonical_anchor = source_order_index
residual/carry anchor   = origin source_order_index
ENTRY-origin residual   = -1
```

### Visual order

Store explicit integer `visual_order`.

Biological nodes inherit deterministic DFS/source ordering.

Residuals are placed after the visible descendants belonging to their source branch; subtype order is:

```text
other_hidden
assigned_above
```

Entry-origin residuals follow all visible first-rank biological nodes.

Carry nodes preserve the lane's established relative order.

The independent verifier reconstructs exact `visual_order`.

Use this total-order algorithm; prose such as "DFS/source ordering" is not
itself a contract:

```text
1. For each selected biological column, sort children by
   (-integer clade, path UTF-8 bytes, source_order_index, id).
2. Sort roots with the same key and visit each root in pre-order DFS.
   This produces biological indices 0..B-1.
3. For a residual/carry lane, let D be the biological descendants of its
   origin source in this column. Its insertion slot is max(D.index + 1), or,
   when D is empty, the bisect-right position of the origin path in the same
   canonical pre-order path list. Entry-origin lanes use slot B.
4. Sort a column's items by:
      (insertion_slot, item_kind_marker, subtype_order, lane_id, id)
   where biological items have item_kind_marker=1, residual/carry items have
   item_kind_marker=0, and subtype_order is other_hidden=0,
   assigned_above=1. A carry uses its lane's original subtype_order.
5. Assign visual_order as the resulting zero-based sequence. ENTRY alone has
   visual_order=0. No ties remain because lane_id and id are final keys.
```

The canonical source path comparison in steps 1 and 3 is bytewise UTF-8, not
locale-aware. Carry lanes retain their origin slot and lane key at every later
column, so adding a new residual lane cannot reorder an existing lane.

### Link ordering

Link kind order:

```text
biological     = 0
other_hidden   = 1
assigned_above = 2
carry          = 3
```

Sort by:

```text
transition_index
canonical source node order
canonical target node order
kind_order
id
```

Reject duplicate `(source, target, kind)`.

No zero-valued links.

---

## 11. Canonical global invariants

For `C > 0`:

```text
entry.value == C
sum(entry outgoing) == C

for each real non-final source:
    sum(outgoing) == source.clade

for each residual non-final source:
    exactly one outgoing carry
    outgoing.value == source.value

for each real target:
    incoming biological-path value == target.clade

for every selected-rank column:
    sum(node.value) == C

for every transition:
    sum(link.value) == C

final selected column:
    sum(node.value) == C

rightmost_flow == C
```

Separate sample accounting:

```text
total == classified + unclassified
```

Unclassified remains outside taxonomy flow.

---

## 12. Zero-classified behavior

For:

```text
classified = 0
```

canonical JSON contains:

```text
entry.value = 0
default_view.nodes = []
default_view.links = []
all conservation totals = 0
```

Do not emit zero-width biological or residual nodes.

HTML displays:

```text
No classified reads available for taxonomy flow.
```

and still shows Total/Classified/Unclassified summaries.

---

## 13. Configuration

```yaml
pavian:
  enabled: false
  render_html: true

  sankey:
    enabled: true
    render_html: true
    ranks: ["D", "K", "P", "C", "O", "F", "G", "S"]
    max_taxa_per_rank: 10
```

Rules:

```text
ranks:
  unique canonical subsequence of D,K,P,C,O,F,G,S
  at least 2 ranks

max_taxa_per_rank:
  integer 1..100
```

`pavian.enabled=false` means effective Sankey execution state is false even if default child settings remain present in YAML.

Compute and persist the effective states once; do not let individual output
builders re-read raw child flags:

```text
pavian_enabled_effective = isTRUE(pavian.enabled)
pavian_html_effective = pavian_enabled_effective && isTRUE(pavian.render_html)
sankey_enabled_effective =
  pavian_enabled_effective && isTRUE(pavian.sankey.enabled)
sankey_html_effective =
  sankey_enabled_effective &&
  pavian_html_effective &&
  isTRUE(pavian.sankey.render_html)
```

The manifest, Pavian provenance, output census, and release verifier must all
record these effective values. In particular, the renderer gate is exactly
`pavian_enabled_effective && sankey_enabled_effective`; the child setting alone
must never activate Sankey output. JSON may be emitted when
`sankey_enabled_effective` is true and `sankey_html_effective` is false, while
HTML path/status fields must be `null`/`not_requested` in that case.

---

## 14. Trust boundary: verifier settings must not come from output

The previous proposal used `observed["defaults"]` to tell the verifier what ranks and N to expect. That is self-referential.

Trusted expected values must be reconciled from:

```text
resolved_config.yml
run_manifest.json
pavian_provenance.json
```

Then invoke the independent verifier with explicit arguments:

```text
--expected-ranks D,K,P,C,O,F,G,S
--expected-max-taxa-per-rank 10
--expected-renderer-version 0.4.8
```

The verifier rejects any `.sankey.json` whose own defaults differ from those trusted values.

---

## 15. Manifest: introduce schema revision 4

Do not silently expand revision 3. Rev3 remains the v0.4.7 contract.

v0.4.8:

```text
schema_version: 2
schema_revision: 4
```

Extend:

```json
"exports": {
  "pavian": {
    "enabled": true,
    "render_html": true,
    "provenance_path": "07_Kreport/pavian/pavian_provenance.json",
    "integration": "official_pavian_upload_plus_builtin_taxonomy_viewers",
    "official_pavian_compatibility": "kraken_report_input_contract_only",
    "sankey": {
      "enabled": true,
      "render_html": true,
      "ranks": ["D","K","P","C","O","F","G","S"],
      "max_taxa_per_rank": 10
    }
  }
}
```

Invariants:

```text
Pavian false => Sankey enabled/render_html false
Sankey render_html true => Sankey enabled true
Sankey enabled => kreport module requested
ranks canonical
N integer 1..100
```

Prior manifest revisions remain readable unchanged.

---

## 16. Pavian provenance schema v2

Bump `pavian_provenance.json` schema to 2.

Add explicit taxonomy-resolution identity:

```json
"taxonomy_resolution_path": "07_Kreport/taxonomy_resolution.tsv",
"taxonomy_resolution_sha256": "..."
```

Add:

```json
"explorer": {
  "renderer": "builtin_kraken_report_explorer",
  "renderer_version": "0.4.8"
},
"sankey": {
  "enabled": true,
  "render_html": true,
  "renderer": "builtin_taxonomy_sankey",
  "renderer_version": "0.4.8",
  "schema_version": 1,
  "ranks": ["D","K","P","C","O","F","G","S"],
  "max_taxa_per_rank": 10,
  "count_model": "classified clade-read flow with persistent explicit residual lanes"
}
```

Per sample record:

```text
kreport path/hash
pavian JSON path/hash
pavian HTML path/hash
Sankey JSON path/hash
Sankey HTML path/hash
```

Nullability must match effective enabled/render state exactly.

---

## 17. `.sankey.json` schema

Recommended top-level object:

```json
{
  "schema_version": 1,
  "sample_id": "S1",
  "renderer": "builtin_taxonomy_sankey",
  "renderer_version": "0.4.8",
  "source_renderer": "builtin_kraken_report_explorer",
  "source_schema_version": 1,
  "source_payload_sha256": "...",
  "count_model": "classified clade-read flow with persistent explicit residual lanes",
  "denominator": "TotalReads",
  "totals": {
    "total": 100000,
    "classified": 92000,
    "unclassified": 8000
  },
  "defaults": {
    "ranks": ["D","K","P","C","O","F","G","S"],
    "max_taxa_per_rank": 10
  },
  "rank_order": {},
  "source_nodes": [],
  "entry": {
    "id": "synthetic:entry:classified",
    "value": 92000
  },
  "default_view": {
    "nodes": [],
    "links": [],
    "conservation": {
      "classified": 92000,
      "column_totals": [],
      "transition_totals": [],
      "rightmost_flow": 92000
    }
  }
}
```

`source_nodes` remain embedded so the offline browser can construct transient alternate rank/top-N views without network access.

### Canonical byte and shape contract

The producer and independent verifier implement the same serialization contract
without sharing producer functions:

```text
encoding: UTF-8, no BOM
object member order: the schema order below; never alphabetical sorting
array order: explicitly defined by the containing field
separators: comma and colon only; no insignificant spaces
Unicode: emit non-ASCII code points literally; escape only JSON-required
         controls, quotation marks, and backslashes
newline: exactly one LF after the complete JSON document
numbers: decimal, base-10, integer tokens only for all count/order fields
forbidden: NaN, Infinity, -0, exponent notation, fractional counts
```

The canonical writer must reject any integer outside
`0 <= value <= 9007199254740991` (`2^53 - 1`) before writing JSON. The browser
must validate every parsed count with `Number.isSafeInteger` and apply the same
bound; silently rounding unsafe values is forbidden.

The exact top-level key order is:

```text
schema_version, sample_id, renderer, renderer_version, source_renderer,
source_schema_version, source_payload_sha256, count_model, denominator, totals,
defaults, rank_order, source_nodes, entry, default_view
```

The exact nested shapes are:

```text
rank_order:
  array of selected rank codes in displayed-column order; ENTRY is implicit
source_nodes:
  array in canonical source order; each record contains
  id, taxid, rank, name, path, parent_id, clade, direct, status,
  source_order_index
default_view.conservation:
  classified, column_totals, transition_totals, rightmost_flow
column_totals:
  [{"rank": <rank>, "value": <integer>}] in rank_order
transition_totals:
  [{"from": <rank>, "to": <rank>, "value": <integer>}] in transition order
```

Node records use one fixed key order and explicit `null` values for fields that
do not apply to that node kind:

```text
id, kind, subtype, biological, carried, lane_id, taxid, rank, name, path,
parent_id, clade, direct, status, selection_reason, source_order_index,
visual_order, origin_source_id, origin_target_rank, column_rank, value,
member_count, member_paths_sha256
```

Link records always use:

```text
id, source, target, kind, value, transition_index
```

The verifier compares both the parsed object and the exact canonical bytes.
The CI release fixture uploads the canonical `.sankey.json` from Ubuntu and
Windows and fails if their SHA-256 bytes differ. HTML must embed those exact
JSON bytes, not a separately reserialized object.

---

## 18. Independent verifier must compare the entire object

Add:

```text
analysis/utils/verify_taxonomy_sankey_correspondence.py
```

It must not import the producer.

It independently reconstructs the entire expected object from:

```text
source .pavian.json
trusted ranks
trusted max_n
trusted renderer version
```

Then:

```python
if observed != expected:
    fail("Sankey JSON differs from independently reconstructed object")

if observed_bytes != canonical_json_bytes(expected):
    fail("Sankey JSON is not canonical byte serialization")
```

This comparison includes:

```text
top-level metadata
source hash
totals
defaults
rank_order
source_nodes
entry
node IDs
link IDs
synthetic metadata
selection reasons
visual_order
conservation summaries
```

For HTML:

```text
embedded payload bytes == .sankey.json bytes
CSP valid
no forbidden network capability
```

---

## 19. Release verifier trust chain

`tests/verify_release_run.R` must:

```text
load resolved_config.yml
load run_manifest.json
load pavian_provenance.json
cross-check Sankey enabled/render state
cross-check ranks
cross-check max_n
verify taxonomy_resolution.tsv path/hash
verify every per-sample Sankey path/hash
reconcile with owned_outputs/artifacts/physical census
invoke Python verifier with trusted settings
```

Conceptual call:

```r
args <- c(
  sankey_verifier,
  "--payload", pavian_json,
  "--sankey-json", sankey_json,
  "--expected-ranks",
  paste(manifest_ranks, collapse = ","),
  "--expected-max-taxa-per-rank",
  as.character(manifest_max_n),
  "--expected-renderer-version",
  manifest$pipeline_version
)
```

Add `--sankey-html` only when effective render_html is true.

---

## 20. Existing-output policy

`taxonomy_sankey_renderer.py` has **no per-file overwrite mode** in v0.4.8.

If JSON or HTML destination exists or is a symlink:

```text
fail closed
```

The pipeline's existing top-level `--overwrite` remains a **whole-run transactional publication** policy. Module builders operate inside fresh private staging.

This avoids a second, weaker ownership model.

Standalone developer use must choose new paths or deliberately remove old developer artifacts.

Before final publication, re-check that the destination did not appear unexpectedly.

---

## 21. New files

```text
analysis/utils/taxonomy_sankey_renderer.py
analysis/utils/verify_taxonomy_sankey_correspondence.py
analysis/utils/taxonomy_sankey_client.js

tests/test_taxonomy_sankey_renderer.py
tests/test_taxonomy_sankey_correspondence.py
tests/test_taxonomy_sankey_client.mjs
```

The JS file contains DOM-free transient graph-model logic. The Python HTML builder embeds this exact committed JS inline.

---

## 22. Transient browser views

Canonical default view:

```text
Python-built
independently verified
release-grade
```

Interactive changes to:

```text
rank selection
max taxa per rank
```

produce a transient browser view.

UI state must clearly show:

```text
Verified default view
```

or:

```text
Interactive transient view
```

Reset restores `payload.default_view`.

### Search

Search only highlights/dims nodes. It must not remove graph mass.

### Branch collapse

Defer true branch collapse from the first v0.4.8 implementation. A collapse operation requires its own residual transformation to preserve conservation. Do not ship a visually convenient but mathematically ambiguous collapse.

---

## 23. Test the actual transient JS model

Use pinned Node in CI **only for testing**.

`tests/test_taxonomy_sankey_client.mjs` runs the same committed DOM-free JS that is embedded in HTML.

Matrix:

```text
ranks:
  full
  D,P,G,S
  P,G
  K,F,S
  G,S

max_n:
  1
  2
  10
  40
```

For every transient view:

```text
entry outflow == classified
every real non-final source outflow == clade
every residual non-final source has one equal carry
every column total == classified
every transition total == classified
final column total == classified
no zero-valued link
no duplicate node/link ID
all endpoints exist
```

Also require JS reconstruction of the configured default to be semantically identical to Python `default_view` on shared fixtures.

---

## 24. Required regression fixtures

### First-rank top-N

```text
C = 100
Domain A = 60
Domain B = 25
Domain C = 15
max_n = 2
```

Expected:

```text
ENTRY -> A = 60
ENTRY -> B = 25
ENTRY -> Other Domain = 15
```

The 15 persists through every later column.

### First selected rank below Domain

```text
C = 100
80 reads reach Phylum
20 terminate above Phylum
selected ranks = P,G
```

Expected:

```text
ENTRY -> Phylum partition = 80
ENTRY -> Assigned above Phylum = 20
```

At Genus, the 20 is carried.

### Intermediate residual regression

```text
C = 100
Domain = 100
visible Kingdom = 60
hidden Kingdom = 40
visible Phylum = 60
```

Expected Phylum column:

```text
visible P = 60
Other K carry = 40
total = 100
```

This exact Luna-discovered case gets a named test.

### Internal direct count

```text
Genus clade = 100
Species A = 50
Species B = 30
direct at Genus = 20
```

Expected:

```text
G -> A = 50
G -> B = 30
G -> Assigned above Species = 20
```

### Ancestral closure

Construct a lower-rank top-N taxon whose ancestor is not top-N at the previous rank. The ancestor must be added explicitly; the lower node may never emerge from an Other lane.

### Other required fixtures

```text
equal abundance ties
unresolved TaxID=0
conflicted node with nonzero winner
zero classified
Unicode
hostile HTML-like labels
same display name under different parent paths
rank subsets
```

### Configuration and source-validation fixtures

These are release-blocking regression fixtures, not only documentation
examples:

```text
pavian.enabled=false, pavian.sankey.enabled=true:
  no Sankey JSON/HTML, effective sankey.enabled=false,
  effective sankey.render_html=false, provenance paths null
pavian.enabled=true, pavian.sankey.enabled=false:
  no Sankey JSON/HTML and the same explicit disabled-state record
pavian.render_html=false, pavian.sankey.render_html=true:
  Sankey JSON may exist, Sankey HTML must not, and effective HTML=false
missing biological ancestor:
  source validation fails before rendering; no guessed parent or Other lane
parent rank/path mismatch, duplicate source ID, invalid rank, negative or
fractional count, and source payload/hash mismatch:
  each fails closed with a deterministic diagnostic
```

The disabled-state tests must inspect manifest, provenance, physical output
census, and verifier behavior together. The missing-ancestor tests must prove
that the producer and independent verifier reject the same invalid source
without silently repairing it.

---

## 25. Generated stress fixture

Generate, do not commit, a large source payload with approximately:

```text
5,000–10,000 taxonomy nodes
many shared ancestors
many >N branches at several ranks
equal-count ties
internal direct counts
unresolved/conflicted records
Unicode
duplicate display names under distinct paths
```

Assertions:

```text
no recursion failure
deterministic repeated output
all conservation invariants
bounded visible default graph
max_n=40 transient view valid
full source payload canonical
```

Avoid brittle wall-clock pass/fail thresholds in normal CI.

---

## 26. Mutation tests

Fail closed on mutation of any of:

```text
sample_id
renderer identities/version
source schema/hash
count model
denominator
totals
defaults
rank_order
source_nodes
entry
biological ID/path/rank/TaxID/status/direct/clade
selection_reason
source_order_index
visual_order
residual subtype/lane/origin/value/member digest
carried flag
link ID/source/target/kind/value/transition
conservation records
HTML embedded payload
CSP/network capabilities
```

Also reject:

```text
duplicate node IDs
duplicate link IDs
unknown endpoints
negative/fractional counts
unsafe JS integer
unknown synthetic subtype
carry value changes
residual re-entry into biological target
zero-valued links
```

---

## 27. Security/offline contract

CSP:

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

No:

```text
fetch
XMLHttpRequest
WebSocket
EventSource
sendBeacon
external script/style
http://
https://
eval
new Function
```

User-derived text goes through `textContent` / SVG text nodes, never executable HTML.

Hostile labels such as `<script>alert(1)</script>` must remain inert text.

---

## 28. Deterministic SVG layout

No general Sankey solver.

X:

```text
ENTRY, R0, R1, ..., Rm-1
```

Y:

```text
explicit deterministic visual_order
```

Height/ribbon thickness:

```text
proportional to integer node/link value under one common view scale
```

Ribbons:

```text
deterministic cubic Bézier closed SVG paths
```

No iterative relaxation, random jitter, or persisted drag state.

Residual carry lanes are visually subdued and explicitly non-biological.

---

## 29. Required HTML interactions

v0.4.8 minimum:

```text
Sankey canvas
rank labels
node labels
count/percentage toggle
hover tooltip
click-to-pin
ancestor/descendant highlight
search highlight
rank controls
top-N control
Reset to verified default
unresolved/conflicted indicators
download canonical Sankey JSON
download current SVG
responsive sizing
keyboard-accessible node selection
```

Defer:

```text
branch collapse
dragging
multi-sample side panels
remote NCBI links
PNG export
server-backed filtering
```

---

## 30. Rough config diff

```diff
 pavian:
   enabled: false
   render_html: true
+  sankey:
+    enabled: true
+    render_html: true
+    ranks: ["D", "K", "P", "C", "O", "F", "G", "S"]
+    max_taxa_per_rank: 10
```

Validate canonical rank subsequence and N 1..100.

---

## 31. Rough manifest diff

```diff
 manifest <- list(
   ...
-  schema_revision = 3L,
+  schema_revision = 4L,
   ...
   exports = list(
     ...
     pavian = list(
       enabled = pavian_enabled,
       render_html = pavian_html_enabled,
       provenance_path = ...,
-      integration =
-        "official_pavian_upload_plus_builtin_kraken_report_explorer",
+      integration = if (sankey_enabled_effective) {
+        "official_pavian_upload_plus_builtin_taxonomy_viewers"
+      } else {
+        "official_pavian_upload_plus_builtin_kraken_report_explorer"
+      },
       official_pavian_compatibility =
         "kraken_report_input_contract_only",
+      sankey = list(
+        enabled = sankey_enabled_effective,
+        render_html = sankey_html_effective,
+        ranks = json_array(cfg$pavian$sankey$ranks),
+        max_taxa_per_rank =
+          as.integer(cfg$pavian$sankey$max_taxa_per_rank)
+      )
     )
   )
 )
```

Rev3 validator remains unchanged.

---

## 32. Rough Module 07 integration

After `.pavian.json` has been successfully built and independently verified:

```diff
+pavian_enabled_effective <- isTRUE(pavian_cfg$enabled)
+pavian_html_effective <- pavian_enabled_effective && isTRUE(pavian_cfg$render_html)
+sankey_enabled_effective <-
+  pavian_enabled_effective && isTRUE(pavian_cfg$sankey$enabled)
+sankey_html_effective <-
+  sankey_enabled_effective && pavian_html_effective &&
+  isTRUE(pavian_cfg$sankey$render_html)
+
+if (isTRUE(sankey_enabled_effective)) {
+  sankey_dir <- file.path(pavian_dir, "sankey")
+  dir.create(sankey_dir, recursive = TRUE, showWarnings = FALSE)
+
+  builder <- file.path(
+    cfg$pipeline_root,
+    "analysis", "utils",
+    "taxonomy_sankey_renderer.py"
+  )
+
+  sankey_json <- file.path(
+    sankey_dir,
+    sprintf("%s.sankey.json", sample_filename)
+  )
+
+  sankey_html <- if (isTRUE(sankey_html_effective)) {
+    file.path(
+      sankey_dir,
+      sprintf("%s.sankey.html", sample_filename)
+    )
+  } else NULL
+
+  args <- c(
+    builder,
+    "--payload", json_out,
+    "--json-out", sankey_json,
+    "--ranks", paste(pavian_cfg$sankey$ranks, collapse = ","),
+    "--max-taxa-per-rank",
+    as.character(pavian_cfg$sankey$max_taxa_per_rank)
+  )
+
+  if (!is.null(sankey_html)) {
+    args <- c(args, "--html-out", sankey_html)
+  }
+
+  result <- processx::run(
+    python_cmd,
+    args,
+    error_on_status = FALSE
+  )
+
+  if (!identical(result$status, 0L)) {
+    stop(sprintf(
+      "Builtin taxonomy Sankey failed for sample '%s': %s",
+      sample_record$sample_id,
+      trimws(paste(result$stderr, result$stdout))
+    ), call. = FALSE)
+  }
+
+  all_outputs <- c(
+    all_outputs,
+    sankey_json,
+    if (!is.null(sankey_html)) sankey_html else character(0)
+  )
+}
```

---

## 33. Rough renderer core

```python
def build_view(source_payload, selected_ranks, max_n):
    C = int(source_payload["totals"]["classified"])

    if C == 0:
        return empty_view(C, selected_ranks)

    nodes = validate_source_nodes(source_payload)

    retained = select_top_n(
        nodes,
        selected_ranks,
        max_n,
    )
    retained = apply_ancestral_closure(
        nodes,
        retained,
        selected_ranks,
    )

    view_nodes = []
    links = []

    build_entry_transition(
        C,
        selected_ranks[0],
        nodes,
        retained,
        view_nodes,
        links,
    )

    for i, (source_rank, target_rank) in enumerate(
        zip(selected_ranks, selected_ranks[1:]),
        start=1,
    ):
        build_biological_transition(
            source_rank,
            target_rank,
            nodes,
            retained,
            view_nodes,
            links,
            transition_index=i,
        )

        carry_all_existing_residual_lanes(
            target_rank,
            view_nodes,
            links,
            transition_index=i,
        )

    assign_visual_order(
        view_nodes,
        nodes,
        selected_ranks,
    )
    canonical_sort_nodes(view_nodes, selected_ranks)
    canonical_sort_links(links, view_nodes)

    conservation = assert_conservation(
        C,
        selected_ranks,
        view_nodes,
        links,
    )

    return {
        "nodes": view_nodes,
        "links": links,
        "conservation": conservation,
    }
```

`carry_all_existing_residual_lanes()` must never duplicate a node already emitted in the target column.

---

## 34. Rough verifier CLI

```diff
 parser.add_argument("--payload", required=True, type=Path)
 parser.add_argument("--sankey-json", required=True, type=Path)
 parser.add_argument("--sankey-html", type=Path)
+parser.add_argument("--expected-ranks", required=True)
+parser.add_argument(
+    "--expected-max-taxa-per-rank",
+    required=True,
+    type=int
+)
+parser.add_argument(
+    "--expected-renderer-version",
+    required=True
+)
```

Then reconstruct entire object independently and compare object + canonical bytes.

---

## 35. Rough provenance diff

```diff
 pavian_provenance <- list(
-  schema_version = 1L,
+  schema_version = 2L,
   path_basis = "run_dir",
+  taxonomy_resolution_path =
+    krona_artifact_relpath(
+      res_summary_file,
+      cfg$output$base_dir
+    ),
+  taxonomy_resolution_sha256 =
+    compute_file_hash(res_summary_file),
+  explorer = list(
+    renderer = "builtin_kraken_report_explorer",
+    renderer_version = pipeline_version
+  ),
+  sankey = list(
+    enabled = sankey_enabled_effective,
+    render_html = sankey_html_effective,
+    renderer = "builtin_taxonomy_sankey",
+    renderer_version = pipeline_version,
+    schema_version = 1L,
+    ranks = json_array(cfg$pavian$sankey$ranks),
+    max_taxa_per_rank =
+      as.integer(cfg$pavian$sankey$max_taxa_per_rank),
+    count_model =
+      "classified clade-read flow with persistent explicit residual lanes"
+  ),
   samples = json_array(pavian_records)
 )
```

---

## 36. Rough release-verifier flow

```r
manifest_sankey <- manifest$exports$pavian$sankey
resolved_sankey <- resolved_cfg$pavian$sankey
prov_sankey <- pavian_prov$sankey

assert_same_rank_vector(
  manifest_sankey$ranks,
  resolved_sankey$ranks,
  prov_sankey$ranks
)

stopifnot(
  identical(
    manifest_sankey$max_taxa_per_rank,
    resolved_sankey$max_taxa_per_rank
  ),
  identical(
    manifest_sankey$max_taxa_per_rank,
    prov_sankey$max_taxa_per_rank
  )
)

resolution <- resolve_run_artifact(
  root,
  pavian_prov$taxonomy_resolution_path
)

stopifnot(
  identical(
    sha256_file(resolution),
    pavian_prov$taxonomy_resolution_sha256
  )
)
```

Then invoke the Python verifier using the manifest-derived trusted values, not the observed Sankey output.

---

## 37. CI amendments

First fix v0.4.7.

For v0.4.8 add:

```diff
 - name: Run unit and regression tests
   run: |
     ...
+    python -m unittest -v tests/test_taxonomy_sankey_renderer.py
+    python -m unittest -v tests/test_taxonomy_sankey_correspondence.py

+ - uses: actions/setup-node@v4
+   with:
+     node-version: '24.11.1'

+ - name: Test transient Sankey model
+   run: node --test tests/test_taxonomy_sankey_client.mjs
```

Add integration fixture with non-default:

```yaml
pavian:
  enabled: true
  render_html: false
  sankey:
    enabled: true
    render_html: false
    ranks: ["P", "G"]
    max_taxa_per_rank: 2
```

Run repeated-build hash checks in one target environment.

The implementation must replace the illustrative action tag with the full
commit SHA for the reviewed `actions/setup-node` release (and pin every other
action in the workflow the same way). A major-only Node selector such as `24`
is not acceptable. Record the action tag-to-SHA mapping in the CI change.

Add a cross-OS comparison job that downloads the Ubuntu and Windows fixture
artifacts and requires byte-identical SHA-256 values for canonical
`.sankey.json` (and the embedded payload portion of HTML). A semantic-only
cross-OS comparison is insufficient for the canonical serialization contract.

Add producer and client tests for the numeric boundary: `2^53 - 1` is valid,
`2^53` and above are rejected, and negative, fractional, exponent-form, NaN,
or Infinity values fail closed. JavaScript must not rely on unsafe `Number`
rounding; this plan chooses the bounded-safe-integer contract rather than
BigInt JSON.

---

## 38. Determinism

Require same-target repeated builds:

```text
sankey.json SHA256 run1 == run2
sankey.html SHA256 run1 == run2
```

No timestamps, random IDs, locale sorting, or browser-generated canonical state.

Semantic equality and canonical-byte equality are required across OSes. The
release fixture comparison in Section 37 is mandatory, not a best-effort
check.

---

## 38A. Documentation and release-note work

Before implementation is considered complete, update all of:

```text
README.md:
  built-in Taxonomy Sankey description, offline/security boundary,
  effective configuration rules, artifact paths, and explicit distinction
  from the real Pavian service
config reference / example:
  pavian.enabled, pavian.render_html, sankey.enabled,
  sankey.render_html, ranks, and max_taxa_per_rank with effective-state
  examples
CLI/help text:
  Sankey request, disabled behavior, and JSON-vs-HTML output semantics
CHANGELOG.md and release notes:
  v0.4.8 schema/provenance/renderer changes, deterministic-byte contract,
  cross-OS verification, and the no-network guarantee
```

Documentation examples must match the disabled-state and non-default-rank
fixtures. Do not describe the built-in renderer as the real Pavian service or
claim Pavian-compatible Sankey internals.

---

## 39. Production acceptance gate

Do not sign off until all are true:

```text
v0.4.7 baseline fixed and released
manifest revision 4 implemented
Pavian provenance schema v2 implemented
taxonomy_resolution.tsv hash recorded/verified
explicit Sankey effective state recorded
trusted config/provenance/manifest cross-checks
entry-node model implemented
first-rank residual model implemented
persistent residual carry implemented
ancestral closure implemented
path-based IDs exact
synthetic IDs exact
canonical synthetic/link order exact
entire Sankey JSON independently reconstructed
HTML embedded payload byte-equal to JSON
no per-artifact overwrite mode
generated large-branch stress fixture passes
transient JS conservation matrix passes
zero-classified case passes
hostile-label/security tests pass
Ubuntu + Windows CI green
release verifier passes exact run artifacts
disabled-state and missing-ancestor/source-validation fixtures pass
safe-integer boundary tests pass
canonical JSON bytes match across Ubuntu and Windows
README/config/CLI/CHANGELOG release documentation is updated
```

---

## 40. Final implementation recommendation

The corrected architecture is:

```text
abundance table
    ↓
.kreport + taxonomy_resolution.tsv
    ↓
independent v0.4.7 explorer correspondence verification
    ↓
canonical <sample>.pavian.json
    ↓
trusted resolved config / manifest / provenance
    ↓
builtin_taxonomy_sankey
    ↓
canonical <sample>.sankey.json
    ↓
independent full-object Sankey verifier
    ↓
self-contained <sample>.sankey.html
```

The defining visualization invariant is now:

```text
EVERY SELECTED COLUMN
==
EVERY TRANSITION
==
CLASSIFIED READS
```

because the graph explicitly contains:

```text
visible biological flow
+ Other hidden flow
+ Assigned-above-rank flow
+ persistent carry of all prior residual flow
```

No interactive control may make classified mass disappear.

That rule should be treated as the Sankey renderer's equivalent of the pipeline's read-accounting contract.
