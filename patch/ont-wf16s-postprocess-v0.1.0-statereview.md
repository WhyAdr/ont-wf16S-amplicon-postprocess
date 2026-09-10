# `ont-wf16S-amplicon-postprocess` v0.1.0 state review

**Review date:** 2026-09-05

**Code snapshot:** `58341edaa1460caf7d540f3eb1fb818ce61fc9a1` (`main`)

**Additional producer outputs examined:** five local, untracked `wf-16s` v1.6.1 runs

**Review method:** source inspection, exact streaming recounts, direct parser probes,
two temporary full-pipeline smoke runs, official producer-contract review, and an
independent data probe by the Luna subagent

## 1. Executive assessment

The project now has a coherent and well-tested **narrow core**: it can post-process
an NCBI-taxonomy, minimap2-classified, species-rank `wf-16s` abundance table plus
its five-field per-read assignment file. It validates inputs strictly, separates
raw alignment status from effective classification, provides single-sample and
cohort analytics with explicit gates, exports Kraken/Pavian-compatible reports,
and records substantial run provenance. The tracked Ambar Ayunda fixture and two
new NCBI/minimap2 datasets reconcile exactly.

It is not yet a general post-processor for every output that `wf-16s` v1.6.1 can
produce. The newly supplied files demonstrate three distinct upstream contracts:

1. NCBI + minimap2 + species abundance: supported.
2. SILVA + minimap2 + genus abundance: rejected because it has seven ranks, not
   the hard-coded eight-rank species schema.
3. NCBI + Kraken2/Bracken + species abundance: the six-field assignment file is
   rejected, while the abundance table alone can be accepted with a biologically
   incorrect denominator interpretation.

The new data also expose a release-blocking single-sample composition bug that
the existing release integration test did not detect:

- `analysis/04_taxa_composition.R` drops the sample name when extracting the
  one-column unclassified count.
- `classification_fraction.tsv` therefore writes `UnclassifiedReads=NA`.
- `04_classification_fraction.png` drops the missing row and renders the sample
  as 100% classified.
- The classified fraction itself, the QC reconciliation table, alpha-diversity
  denominators, and `.kreport` arithmetic remain correct.

**Current release judgment:** do not create the `v0.1.0` tag until the
single-sample composition defect is fixed and its output content—not just its
existence—is asserted. The implementation is close to a sound narrow v0.1.0,
but this result changes the state from “green release candidate” to “one P0
correctness patch required.”

## 2. What the project provides now

| Area | Current capability | Important boundary |
|---|---|---|
| Configuration and CLI | YAML configuration; root-independent path resolution; module selection; output override; `--validate-only`; fail-fast default; `--keep-going`; `--overwrite`; explicit `--refresh-taxonomy` | No automatic discovery of files inside a `wf-16s` result directory |
| Input validation | Exact abundance schema, numeric/integer/non-negative counts, aggregate-column reconciliation, one canonical unclassified row, safe sample IDs, metadata alignment, unique assignment read IDs, status/TaxID/length validation | Abundance is hard-coded to eight ranks and assignments to the minimap2 five-field layout |
| Analytical modes | Automatic or explicit single/cohort mode; cohort-only modules produce structured skip records in single mode | Added outputs are separate technical reruns, not a biological cohort |
| QC | Per-read classification reconciliation; raw `C` and `U`; `C + TaxID 0` filter failures; effective `TaxID > 0` classification; read-length summaries and four plots | Requires assignments; composition's separate classification plot currently has the single-sample bug |
| Alpha diversity | Observed richness, Chao1, Shannon, effective species, Simpson, inverse Simpson, Pielou, Fisher alpha, Berger-Parker, analytical rarefaction, seeded rarefaction resamples | Labels and logic assume species-rank input; estimates are conditional on upstream filtering and taxonomy |
| Beta diversity | Bray-Curtis and binary Jaccard distances, PCoA, optional rarefaction stability, PERMANOVA and paired betadisper | Cohort only; no biologically validated cohort fixture; the primary model is group-focused |
| Composition | Classified-only count and relative-abundance tables at phylum through species; contextualized ambiguous labels; single/cohort plots; genus heatmap | Species schema is mandatory; current single-sample all-read classification visualization is wrong |
| Ordination | Hellinger-transformed PCA and Bray-Curtis NMDS with explicit skip/failure artifacts | Cohort only; ordination is exploratory, not evidence of causal separation |
| Shared taxa | Presence, group prevalence, core taxa, group-unique taxa, and UpSet export | Cohort only and threshold-dependent |
| Kraken/Pavian | Standard six-column `.kreport`, standard rank codes, clade arithmetic validation, offline cache, assignment-derived resolution, optional NCBI refresh, unresolved/conflict reports, per-node resolution provenance | Resolver semantics are NCBI-centric; a SILVA TaxID is not an NCBI TaxID |
| Provenance | Semantic pipeline version, Git commit, command, resolved config, R/Python/package versions, file hashes, module outcomes, taxonomy provenance | `params.json` is hashed but classifier, database, rank, thresholds, and upstream workflow revision are not promoted into structured manifest fields |
| Testing and CI | 139 R assertions, 6 Python tests, deterministic checks, mutation-free validation, a full Ambar run, and Ubuntu/Windows CI | The full-run verifier checks that `classification_fraction.tsv` exists but does not validate its values or plotting warnings |

The legacy root scripts remain for retrospective comparison. New work should use
`analysis/00_run_pipeline.R` and the modular `analysis/` implementation.

## 3. Added `wf-16s` evidence

### 3.1 Provenance and storage status

All five new runs report `wf-16s v1.6.1` and Nextflow revision `e2b77380f8` in
their logs. They use `min_len=1300`, `max_len=1700`, `min_read_qual=10`, and
`abundance_threshold=1`. The minimap2 runs use 90% minimum reference coverage
and 90% minimum identity. The Kraken2 run uses confidence 0.1, species-rank
Bracken, and `bracken_threshold=10`.

The directories total approximately 281.62 MiB and are all untracked:

| Local directory | Classifier/database | Files | Size | Abundance shape | Assignment shape | v0.1.0 result |
|---|---|---:|---:|---|---|---|
| `wf-16s_Helga16SrRNA/` | minimap2 / `ncbi_16s_18s` | 24 | 49.90 MiB | 757 rows, eight ranks, species | 23,099 rows, five fields | Parser and full smoke run pass |
| `wf-16s_Helga-SILVA/` | minimap2 / `SILVA_138_1` | 24 | 82.84 MiB | 384 rows, seven ranks, genus | 23,099 rows, five fields | Rejected at abundance rank validation |
| `wf-16s_JulyBiofilm-NCBI/` | minimap2 / `ncbi_16s_18s_28s_ITS` | 24 | 44.17 MiB | 679 rows, eight ranks, species | 17,920 rows, five fields | Parser and full smoke run pass |
| `wf-16s_July-Biofilm-SILVA/` | minimap2 / `SILVA_138_1` | 24 | 79.45 MiB | 300 rows, seven ranks, genus | 17,920 rows, five fields | Rejected at abundance rank validation |
| `wf-16s_July-Biofilm-Kraken2/` | Kraken2/Bracken / `ncbi_16s_18s_28s_ITS` | 17 | 25.26 MiB | 6 rows, eight ranks, species | 17,920 rows, six fields | Assignment rejected; abundance-only use is semantically unsafe |

These bundles should remain local integration evidence rather than being added
wholesale to Git. Their logs and parameter files contain workstation paths and
EPI2ME instance identifiers, and their size is unnecessary for ordinary CI.
A tracked, checksum-aware local-fixture manifest is preferable.

### 3.2 Direct compatibility probes

The production parser independently accepted and reconciled both new
NCBI/minimap2 datasets:

- Helga: 23,099 assignment rows = 18,155 effectively classified + 4,944
  effectively unclassified.
- July: 17,920 assignment rows = 9,668 effectively classified + 8,252
  effectively unclassified.

Both datasets also completed all applicable single-sample modules in temporary
full runs. Cohort modules skipped as designed. Using the existing Ambar NCBI
cache as a portability smoke test, Helga produced 344 unresolved taxonomy nodes
and 7 conflicts; July produced 207 unresolved nodes and 9 conflicts. Those
numbers are **not suitable regression anchors**, because the cache is not a
dataset-specific or complete representation of either upstream taxdump.

Both full runs emitted the same composition warning:

```text
Removed 1 row containing missing values or values outside the scale range
(`geom_col()`).
```

Inspection showed that the removed row was `UnclassifiedReads=NA`. The existing
Ambar full-run artifact confirms the defect:

```text
SampleID                       TotalReads  ClassifiedReads  UnclassifiedReads  ClassifiedFraction
AmbarAyunda_minimap2_16S      114056      80556            NA                 0.706284632110542
```

The safest repair is to build this table from `context$sample_stats`, which is
already the validated source of truth, instead of independently re-indexing a
dimension-dropping matrix slice.

## 4. What “unclassified” means in these outputs

### 4.1 Three quantities must not be collapsed

For minimap2 outputs, the per-read file preserves a raw alignment status and a
post-filter TaxID:

1. `status=U, TaxID=0`: the read did not align/classify against the reference.
2. `status=C, TaxID=0`: the read aligned, but the assignment was relabelled
   unclassified after minimum reference-coverage and/or percent-identity filters.
3. `TaxID>0`: effectively classified and eligible for the classified abundance
   denominator.

This interpretation agrees with the pinned `wf-16s` v1.6.1 documentation, which
states that minimap2-mapped reads failing `min_ref_coverage` or
`min_percent_identity` are relabelled unclassified. The post-processor's decision
to use `TaxID > 0`, rather than raw `C`, is therefore correct.

### 4.2 Exact read accounting

| Dataset | Abundance total | Raw `C` | Raw `U` | `C + TaxID 0` | `TaxID > 0` | Effective classified | `C0` share of effective unclassified |
|---|---:|---:|---:|---:|---:|---:|---:|
| Ambar NCBI/minimap2 | 114,056 | 89,809 | 24,247 | 9,253 | 80,556 | 70.628% | 27.621% |
| Helga NCBI/minimap2 | 23,099 | 23,007 | 92 | 4,852 | 18,155 | 78.596% | 98.139% |
| Helga SILVA/minimap2 | 23,099 | 23,009 | 90 | 3,930 | 19,079 | 82.597% | 97.761% |
| July NCBI/minimap2 | 17,920 | 17,785 | 135 | 8,117 | 9,668 | 53.951% | 98.364% |
| July SILVA/minimap2 | 17,920 | 17,787 | 133 | 4,376 | 13,411 | 74.838% | 97.050% |
| July NCBI/Kraken2 | 17,920 assignment rows | 17,574 | 346 | 0 | 17,574 | 98.069% raw Kraken classification | 0% |

The new minimap2 datasets differ from Ambar in an important way: almost all of
their effectively unclassified reads did align and were subsequently set to
TaxID 0. This makes “database absent” only one possible explanation. The more
immediate hypotheses are marginal identity, insufficient reference coverage,
ambiguous best hits, taxonomy/reference mapping behavior, or classifier/database
specificity.

Read length alone does not explain the pattern. The upstream data are already
filtered to 1,300–1,700 bp, and classified versus `C0` medians are close:

| Dataset | Median classified length | Median `C0` length | Median raw-`U` length |
|---|---:|---:|---:|
| Helga NCBI | 1,510 bp | 1,471.5 bp | 1,468.5 bp |
| Helga SILVA | 1,509 bp | 1,451 bp | 1,468.5 bp |
| July NCBI | 1,474 bp | 1,472 bp | 1,411 bp |
| July SILVA | 1,476 bp | 1,459 bp | 1,413 bp |

An exact read-ID join to each `bamstats.readstats.tsv.gz` establishes which
minimap2 threshold was missed:

| Dataset | Identity only failed | Reference coverage only failed | Both failed |
|---|---:|---:|---:|
| Helga SILVA | 265 | 3,496 | 169 |
| Helga NCBI | 972 | 3,210 | 670 |
| July SILVA | 455 | 3,682 | 239 |
| July NCBI | 3,536 | 2,863 | 1,718 |

Every `TaxID > 0` minimap2 read passed both configured 90% thresholds, and every
`C0` read failed at least one. Coverage failure dominates three runs; July NCBI
also has a large identity-failure component. This is much more informative than
calling all TaxID-0 reads merely “unclassified.” An optional QC join should
surface these three failure partitions directly.

### 4.3 Kraken2/Bracken needs a different denominator model

The July Kraken2 assignment file has 17,920 reads, of which 17,574 are classified
at some taxonomic rank and 346 are raw `U`. Only 440 classified assignments end
at a named species; most stop at higher ranks and receive labels such as
`unclassified Alphaproteobacteria species`.

Bracken then reports five species with 6,173 estimated reads. Adding the same 346
unclassified reads produces an abundance-table total of 6,519, leaving 11,401
raw Kraken-classified reads outside the species abundance table. The five
Bracken estimates consist of 412 directly assigned reads plus 5,761 redistributed
reads.

Consequences for the current post-processor:

- With the six-field assignments configured, v0.1.0 fails early and safely.
- Without assignments, the eight-rank abundance file passes validation.
- It would then report 6,519 as total reads and 6,173/6,519 = 94.69% classified,
  rather than distinguishing 17,920 input reads, 17,574 Kraken-classified reads,
  and 6,173 Bracken species-level estimated reads.
- Bracken estimates would also be labelled as direct “read counts.”

Until a classifier-aware accounting model is implemented, `params.json` with
`classifier=kraken2` should be rejected explicitly. Silent abundance-only
acceptance is more dangerous than the current six-field parser error.

## 5. Database and classifier sensitivity in the same reads

The Helga NCBI and SILVA files contain the same 23,099 read IDs; the three July
runs contain the same 17,920 read IDs. They therefore provide unusually useful
technical-sensitivity tests.

### 5.1 Minimap2 database transition counts

| Sample | Classified by both | Unclassified by both | SILVA only | NCBI only | Exact genus label among both-classified |
|---|---:|---:|---:|---:|---:|
| Helga | 17,882 | 3,747 | 1,197 | 273 | 5,922 / 17,882 (33.1%) |
| July | 9,116 | 3,957 | 4,295 | 552 | 5,782 / 9,116 (63.4%) |

Exact TaxID agreement between NCBI and SILVA is zero, which is expected rather
than alarming: `wf-16s` documents that SILVA 138.1 uses its own TaxID namespace.
The genus-label disagreement combines real reference coverage differences,
different taxonomic concepts, and nomenclature differences. It must not be
presented as biological turnover.

Raw minimap2 alignment status is almost invariant across the database pairs:
Helga has 23,007 reads `C` in both and 90 `U` in both; July has 17,785 `C` in
both and 133 `U` in both, with only two raw-status transitions per sample. The
large effective-classification differences therefore arise after alignment,
principally from reference identity/coverage filtering.

### 5.2 Dominant taxonomic signals

The following are observed classifier labels, not independently validated
organism detections:

- Helga NCBI/minimap2 is led by `Methylobacillus` (26.3% of classified reads),
  `Flavobacterium` (19.8%), `Pseudomethylobacillus` (11.6%), and `Methylovorus`
  (8.3%). SILVA instead assigns 47.0% to `UBA6140`, followed by
  `Flavobacterium` (18.9%) and `Rhodopirellula` (8.2%). This supports a
  methylotroph-associated hypothesis in the NCBI labeling, but not a functional
  claim about methylotrophy without orthogonal evidence. Broader family totals
  are more concordant: Methylophilaceae is 9,050 SILVA versus 8,775 NCBI,
  Flavobacteriaceae 3,629 versus 3,625, and Pirellulaceae 1,616 versus 1,460.
- July NCBI/minimap2 is led by `Sphingobium` (9.1%), `Azospirillum` (8.2%),
  `Proteiniphilum` (7.9%), `Neoroseomonas` (7.5%), and `Xanthobacter` (6.4%).
  SILVA includes related labels but also places 18.1% in an unclassified `A4b`
  genus-level group.
- July Kraken2/Bracken collapses its represented species pool to five estimates:
  `Aggregatilinea lenta` (42.2%), `Petrimonas sulfuriphila` (25.4%),
  `Proteiniphilum saccharofermentans` (22.3%), `Bdellovibrio bacteriovorus`
  (8.7%), and `Aminivibrio pyruvatiphilus` (1.5%). That profile is not directly
  comparable to minimap2 classified-read proportions because both the estimator
  and denominator differ.

Several July genus totals are nevertheless similar between the two minimap2
databases: `Proteiniphilum` is 785 SILVA versus 759 NCBI, `Petrimonas` 321 versus
339, `Sphingobium` 841 versus 884, and `Azospirillum` 845 versus 791. Those broad
signals are more defensible than terminal labels. Conversely, SILVA assigns
2,433 reads to its `A4b` lineage, Kraken/Bracken estimates 2,605 as
`Aggregatilinea lenta`, and NCBI/minimap2 assigns only 20 reads to
`Aggregatilinea`. This is method sensitivity, not three independent confirmations.
Similarly, only 27 reads were directly assigned to `Bdellovibrio bacteriovorus`
before Bracken estimated 534; neither the label nor the estimate establishes
predatory activity.

These paired runs are best used as **method sensitivity controls**. They must not
be combined as samples in PCA, NMDS, PERMANOVA, or shared-taxa analyses.

## 6. Biological and bioinformatic points to scrutinize

### 6.1 Richness is strongly driven by the low-count tail

| Output | Positive taxa | Singleton taxa | Taxa with ≤10 reads | Reads in taxa with ≤10 reads |
|---|---:|---:|---:|---:|
| Ambar NCBI species | 1,836 | 735 (40.0%) | 1,399 | 4.29% |
| Helga NCBI species | 756 | 348 (46.0%) | 658 | 8.24% |
| Helga SILVA genus | 383 | 139 (36.3%) | 309 | 4.34% |
| July NCBI species | 678 | 276 (40.7%) | 555 | 14.46% |
| July SILVA genus | 299 | 104 (34.8%) | 213 | 4.38% |
| July Kraken/Bracken species | 5 | 0 | 0 | 0% |

Singletons comprise 35–46% of detected taxa in the minimap2 outputs. This does
not prove that they are artifacts—the rare biosphere can be real—but it means
observed richness and especially Chao1 are highly sensitive to sequencing error,
reference redundancy, near-tied alignments, taxonomy splitting, and the chosen
abundance threshold. Exact Ambar alpha values are excellent software regression
anchors; they are not biological truth standards.

Recommended reporting should include threshold-sensitivity summaries (for
example, taxa retained at 1, 2, 5, and 10 reads) while preserving the unfiltered
table. Do not silently select the threshold that gives the preferred ecology.

### 6.2 Relative abundance is not cell abundance

The pipeline correctly uses a classified-read denominator for community
composition, but the resulting proportions remain compositional. They are also
affected by DNA extraction, primer-template mismatches, PCR competition,
sequencing, reference coverage, classifier behavior, and inter-genome variation
in 16S rRNA operon copy number. A read proportion should not be described as
biomass, absolute abundance, or cell fraction without external quantitation and
appropriate correction.

For cross-database sensitivity, report both:

- fraction of all input reads assigned to each broad taxon; and
- composition conditional on the effectively classified pool.

Classified-only normalization is useful within one contract, but it can hide the
large yield difference between July NCBI (54.0%) and SILVA (74.8%).

### 6.3 Species calls need calibrated language

Full-length 16S generally improves resolution over short variable regions, but
resolution is lineage-dependent. Some species share indistinguishable or nearly
indistinguishable 16S sequences, some genomes contain heterogeneous operon
copies, and long-read errors plus database naming can move the best hit. A
species label should be treated as an upstream classifier assignment, not an
isolate identification, pathogenicity result, strain call, or phenotype.

The NCBI/SILVA paired runs make this limitation visible: the same read can be
classified in one database and filtered in another or assigned to a differently
named genus. For claims that matter biologically, inspect alignment identity and
coverage, alternative close hits, lineage consistency, and—when possible—use a
mock community, isolate sequence, targeted marker, or shotgun evidence.

### 6.4 Cohort inference still needs biological validation

The cohort code has reasonable software gates, but the new directories are
single samples reanalysed with alternative methods. They do not validate:

- biological replicate handling;
- batch, extraction, barcode, or run effects;
- negative-control contamination;
- repeated-measures designs;
- the stability of PERMANOVA and betadisper under realistic sparsity; or
- whether ordination separation persists after technical covariates are modeled.

PERMANOVA should be interpreted with dispersion results, adequate replication,
pre-specified design/strata, and effect size—not only a p-value. Differential
abundance and contaminant modeling remain outside v0.1.0.

## 7. Findings by priority

### P0-1 — Single-sample classification composition is wrong

`analysis/04_taxa_composition.R` uses:

```r
unclass_counts <- count_matrix[unclass_idx, ]
```

For a one-column matrix this becomes an unnamed scalar. Later indexing by sample
name returns `NA`. Reuse `context$sample_stats` or preserve and explicitly name
the dimension. Required regression assertions:

- `UnclassifiedReads == 33500` for Ambar;
- `ClassifiedReads + UnclassifiedReads == TotalReads`;
- no `NA` in `classification_fraction.tsv`;
- two plot-input categories exist and sum to the total; and
- a clean full run emits no unexpected plotting warning.

### P0-2 — Kraken2 abundance-only input can pass with false semantics

The current eight-rank parser cannot distinguish direct minimap2 counts from
Bracken rank-level estimates. When assignments are absent, the July Kraken2
table is accepted even though it represents only 6,519 of 17,920 input reads.
For v0.1.0, validate `params.json` and fail with an explicit unsupported-contract
message when `classifier != minimap2`. Full Kraken2 support requires an adapter,
not a relaxed five-to-six-column parser alone.

### P1-1 — Upstream contract identity is not enforced or surfaced

`params.json` is parsed and hashed but its contract-defining fields are not
validated or promoted into the manifest. Record at least:

- upstream workflow name/version/revision when available;
- classifier;
- database set and taxonomy namespace/release;
- selected taxonomic rank;
- minimum length, quality, identity, and coverage;
- abundance threshold; and
- Kraken confidence and Bracken threshold when applicable.

The supplied SILVA runs even request `taxonomic_rank=S` in `params.json` while
producing a genus table, because SILVA 138.1 does not extend below genus in this
workflow. The actual file schema and producer/database contract must be checked
together.

### P1-2 — SILVA and Kraken2 require explicit, separate adapters

SILVA support needs a rank-aware data model and taxonomy-namespace handling.
Kraken2 support needs a six-field assignment parser and separate accounting for
input reads, Kraken assignment rank, and Bracken-estimated abundance. Neither
should be implemented by padding genus rows with a fake species or by discarding
the Kraken k-mer field.

### P1-3 — Taxonomy portability is incomplete

The new runs do not include a post-processor cache. Reusing the Ambar cache lets
NCBI smoke tests complete but leaves hundreds of unresolved nodes. A publication
run should use a cache generated against the same upstream taxonomy release or
ingest the workflow's taxdump directly. SILVA's own TaxIDs must never be sent to
NCBI Entrez or presented as NCBI IDs.

The read-only smoke resolution broke down as follows:

| NCBI run | Total nodes | Source-cache | Assignment-derived | Unresolved | Conflicting lineages |
|---|---:|---:|---:|---:|---:|
| Helga | 1,436 | 497 | 595 | 344 | 7 |
| July | 1,197 | 520 | 470 | 207 | 9 |

Some conflicts are near ties, not trivial noise. For example,
`Mycolicibacterium mageritense` split 28 reads to TaxID 53462 and 26 to TaxID
1209984; `Agrobacterium tumefaciens` split 20 to TaxID 358 and 22 to TaxID
1183401. Retaining the modal mapping is operationally deterministic, but the
near-tie must remain visible and should not be described as an unambiguous
species identity. The probe also found abundance/assignment label transformations
such as `*_Incertae_sedis` versus `unclassified ... family`; reconciliation must
remain lineage-aware rather than depend on literal leaf-name equality.

The upstream taxdump is discoverable but not shipped. Each `wf-16s` output
contains an `output/params.json` whose `database_sets` block records the exact
taxonomy source URI. Both NCBI runs used
`s3://ont-open-data/workflow-databases/wf-metagenomics-dbs/ncbi_16s_18s/new_taxdump_2025-01-01.zip`;
the SILVA runs used
`s3://ont-open-data/workflow-databases/wf-metagenomics-dbs/SILVA_138_1/taxonomy.tar.gz`.
The `prepare_databases` workflow in `wf-metagenomics/modules/local/databases.nf`
downloads and unpacks this archive into Nextflow's `storeDir`
(`/home/prom/epi2melabs/data/{database_set}/`), where `taxonkit` (v0.20.0) reads
`nodes.dmp`, `names.dmp`, and `merged.dmp` via `--data-dir` to produce lineage
strings. The taxdump directory persists on the PromethION but is not copied into
the per-run output directory.

The recommended resolution priority for the post-processor should therefore be:

1. **Taxdump ingest** — parse the workflow's own `nodes.dmp`/`names.dmp` from
   `storeDir` or re-download the exact S3 URI recorded in `params.json`. This is
   the ground truth the classifier operated against and resolves every TaxID the
   classifier could have emitted, with zero version skew.
2. **Assignment-derived** — extract `(taxon_name → TaxID)` mappings from per-read
   assignment files as a supplement for any path not in the taxdump.
3. **NCBI Entrez** — opt-in last resort for merges, renames, or nodes missing
   from the frozen taxdump. Must warn that the online taxonomy may reflect a
   newer release than the one the classifier used, producing lineage
   disagreements.

Placing Entrez last rather than first avoids a silent version mismatch: if the
classifier used the `2025-01-01` freeze and NCBI has since merged or split a
taxon, an online query will return a different lineage from the one the
classifier operated on, and the `.kreport` tree will disagree with the
classification logic that produced it.

### P1-4 — Full integration checks artifact presence more than meaning

`tests/verify_release_run.R` asserts the Ambar QC reconciliation and taxonomy
counts, but only checks that `classification_fraction.tsv` exists. This allowed
an NA-containing table and misleading plot through Ubuntu and Windows CI. Every
high-value table should have schema, finiteness, conservation, and expected-value
assertions.

### P2-1 — Assignment parsing will become memory-heavy for cohorts

The R parser reads the entire assignment file into a character vector, pastes it
into another full string, parses a data frame, and retains every sample's rows
when QC is enabled. This is acceptable for the current 17,920–114,056-read
fixtures, but it scales poorly toward multi-sample, million-read projects. Add a
performance budget and consider streaming summaries, chunked duplicate checks,
or processing one sample at a time.

### P2-2 — The project lacks a truth-known biological fixture

The Ambar and added environmental outputs test producer compatibility and
internal conservation, not taxonomic accuracy. A small mock community with
known expected members, expected absences, and tolerated abundance ranges is the
right integration layer for biological performance.

## 8. Recommended integration-test expansion

Yes—the additional outputs materially strengthen the test strategy, provided
they are used as versioned local fixtures rather than copied into the repository.

### Layer A — fast committed tests on every push

1. Fix and assert exact single-sample composition contents for Ambar.
2. Add minimal producer-native schema snippets for:
   - minimap2 five-field assignments;
   - Kraken2 six-field assignments;
   - eight-rank species abundance; and
   - seven-rank genus abundance.
3. Add a `params.json` contract detector and test supported and unsupported
   combinations with actionable messages.
4. Assert conservation in every summary table and reject `NA`, `NaN`, and
   infinite values unless explicitly allowed.
5. Capture warnings in the full-run test and fail on unexpected warnings.

### Layer B — opt-in local producer integration matrix

Track a small manifest, not the datasets themselves. Suggested fields:

```text
fixture_id  output_dir  wf16s_version  classifier  database_set  rank  expected_result
```

The test should activate through an environment variable such as
`WF16S_INTEGRATION_MANIFEST`; a configured but missing path or changed checksum
should fail, while an unset manifest should skip cleanly.

| Fixture | Expected v0.1 behavior | Assertions |
|---|---|---|
| Ambar NCBI/minimap2 | Full pass | Existing exact 114,056/80,556/33,500 and 46/26 checks plus composition values/no warnings |
| Helga NCBI/minimap2 | Full pass | 23,099/18,155/4,944, `C0=4,852`, finite tables, all promised single-sample artifacts |
| July NCBI/minimap2 | Full pass | 17,920/9,668/8,252, `C0=8,117`, finite tables, all promised single-sample artifacts |
| Helga SILVA/minimap2 | Explicit unsupported failure | Message identifies genus/seven-rank/SILVA boundary, not a generic parse failure |
| July SILVA/minimap2 | Explicit unsupported failure | Same, with exact schema fingerprint |
| July Kraken2/Bracken | Explicit unsupported failure | Message identifies classifier and explains assignment-versus-Bracken denominators |

Do not assert the 344/7 or 207/9 taxonomy diagnostics unless a compatible,
versioned dataset-specific cache is added to the fixture contract.

### Layer C — supported-contract variation

Add synthetic or safely derived cases for:

- abundance-only NCBI/minimap2 input;
- multiple samples with partial and complete assignment mappings;
- all-classified, mostly-unclassified, and zero-usable-classified failures;
- multiple `C0`/raw-`U` mixtures;
- one-to-one assignment/bamstats joins that partition identity-only,
  coverage-only, and both-threshold failures;
- taxa containing ambiguous repeated names under different parents;
- single-row and singleton-heavy communities;
- changed sample order and file order;
- paths with Unicode, spaces, and metacharacters; and
- large assignment files with an explicit memory/time budget.

### Layer D — biological validation

Use a truth-known mock community processed through the complete wet-lab and
`wf-16s` path. Predefine:

- expected taxa at genus and, only where resolvable, species;
- expected negative taxa/contaminants;
- allowable abundance deviation bands;
- minimum classification yield;
- database/classifier matrix; and
- failure criteria for false-positive low-count tails.

This layer should be described as biological performance validation, unlike the
current software regression fixtures.

## 9. Recommended next patch sequence

### Patch 0.1.1 — protect the existing narrow contract

1. Fix the single-sample composition extraction and strengthen the full-run
   verifier.
2. Promote relevant upstream `params.json` fields into validated context and the
   manifest.
3. Explicitly reject unsupported classifier/database/rank combinations before
   analysis; in particular, block Kraken2 abundance-only misinterpretation.
4. Add the optional local-fixture manifest and exact Helga/July NCBI smoke tests.
5. Re-run Ubuntu and Windows CI, then rerun the two local NCBI integrations with
   warning capture.

### Version 0.2 — broaden producer support deliberately

1. Introduce rank-aware abundance objects whose terminal rank is explicit.
2. Add classifier-specific assignment adapters.
3. Model `input reads`, `classifier-positive reads`, `rank-resolved reads`, and
   `estimated abundance` as different quantities.
4. Introduce an explicit taxonomy namespace and release identifier.
5. Generate rank-appropriate metrics and labels (`Observed genera`, not
   `Observed species`) and disable invalid modules cleanly.
6. Validate SILVA and Kraken2 with the supplied native outputs and a mock
   community before declaring support.

## 10. Practical treatment of the unclassified reads

The supplied runs used `output_unclassified=true`, so the most useful next
analysis is evidence-preserving and stratified:

1. Export separate read-ID lists for raw `U` and minimap2 `C0` reads.
2. For `C0`, join read IDs to alignment statistics or BAM-derived identity,
   aligned length, and reference coverage; quantify which threshold each read
   missed and by how much.
3. For raw `U`, inspect length/quality and screen a reproducible subset against
   an alternative curated database or broader sequence search.
4. Compare NCBI and SILVA transition groups: both classified, both unclassified,
   NCBI-only, and SILVA-only.
5. Preserve the original thresholds and results; treat any relaxed-threshold run
   as a sensitivity analysis, not a replacement chosen after seeing preferred
   taxa.
6. Include extraction blanks, negative controls, and a mock community in future
   experiments so low-count detections and unclassified fractions have an
   experimental baseline.

Unclassified reads can reflect missing references or undescribed diversity, but
they can also reflect errors, off-target amplicons, low-complexity sequence,
contamination, ambiguous mappings, or thresholds. The present tables cannot
distinguish those explanations on their own.

## 11. Bottom line

The project has moved from a collection of scripts to a disciplined,
cross-platform post-processing pipeline with strong validation, provenance, and
internal numerical conservation. Its reliable v0.1.0 product boundary is NCBI +
minimap2 + species-rank `wf-16s` output.

The newly added datasets are valuable because they reveal what the original
fixture could not:

- a real single-sample composition bug;
- valid NCBI/minimap2 portability beyond Ambar;
- SILVA genus-rank incompatibility;
- Kraken2/Bracken denominator incompatibility;
- large database-dependent classification shifts in the same reads; and
- the dominance of singleton taxa in richness estimates.

Fix the P0 composition error, enforce the narrow producer contract, and add the
local integration matrix before tagging v0.1.0. Broader SILVA and Kraken2 support
should follow as explicit adapters with namespace- and denominator-aware tests,
not as permissive parsing.

## References

- EPI2ME Labs, [`wf-16s` v1.6.1](https://github.com/epi2me-labs/wf-16s/tree/v1.6.1):
  supported classifiers, database characteristics, output contracts, minimap2
  filter semantics, and Kraken2/Bracken behavior.
- EPI2ME, [Unexpected results, so now what?](https://epi2me.nanoporetech.com/post-meta-analysis/):
  recovery and follow-up analysis of unclassified workflow reads.
- Kembel et al. (2012), [*Incorporating 16S Gene Copy Number Information Improves
  Estimates of Microbial Diversity and Abundance*](https://doi.org/10.1371/journal.pcbi.1002743).
- Gloor et al. (2017), [*Microbiome Datasets Are Compositional: And This Is Not
  Optional*](https://doi.org/10.3389/fmicb.2017.02224).
- Matsuo et al. (2021), [*Full-length 16S rRNA gene amplicon analysis of human gut
  microbiota using MinION nanopore sequencing confers species-level
  resolution*](https://doi.org/10.1186/s12866-020-02094-5): full-length resolution
  gains and primer-associated bias demonstrated with a mock community.
