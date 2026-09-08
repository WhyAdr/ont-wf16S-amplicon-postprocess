# FAPROTAX 1.2.12 integration plan

**Handoff target:** Luna

**Audited base:** `c7a7881` on `main` (`ont-wf16s-postprocess` 0.2.0)

**Recommended release target:** 0.3.0, after the current 0.2.0 release is cut

**Scope:** add an opt-in `faprotax` module using `microeco` and its embedded
FAPROTAX 1.2.12 database. Do not add FAPROTAX2-db, PICRUSt2, differential
testing, or claims of gene/pathway abundance.

## 1. Locked decisions

- Use the current `microeco::trans_func$cal_func(prok_database = "FAPROTAX")`
  API. Target CRAN `microeco >= 2.3.0` and fail unless the embedded database
  reports exactly `1.2.12`.
- Keep the module opt-in for its first release: invoke it with
  `--modules faprotax` or include it in a comma-separated module list. Do not
  add it to the existing default modules until its real-fixture concordance
  review is accepted.
- Check `microeco` only when `faprotax` is requested. Core runs must not fail
  merely because this optional dependency is absent.
- Feed only positive-count Bacteria and Archaea rows. Exclude the canonical
  unclassified row and classified non-prokaryotic rows from mapping, but retain
  both in explicit read accounting.
- Map pipeline `superkingdom` to microeco `Kingdom`. Do not pass the NCBI
  `kingdom` rank (`Bacillati`, `Pseudomonadati`, etc.) as microeco `Kingdom`.
- Use classified taxon counts directly. A function's abundance is the sum of
  counts for taxa mapped to that function; do not rarefy, copy-number-correct,
  or renormalize across functions.
- Functions overlap. Percentages are independent fractions of explicit read
  denominators and therefore must not be stacked or expected to sum to 100%.
- Describe results as **taxon-based functional assignments/inference**. They do
  not establish gene presence, pathway completeness, expression, or activity.

The API and embedded 1.2.12 version are documented in the
[`microeco` tutorial](https://chiliubio.github.io/microeco_tutorial/explainable-class.html).
The current CRAN release is listed on the
[`microeco` package page](https://cran.r-project.org/package=microeco), and the
original database/software remains available from the
[`FAPROTAX` download page](https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/lib/php/index.php?section=Download).

## 2. Output contract

Write to `08_FAPROTAX/`:

| Artifact | Contract |
|---|---|
| `faprotax_function_abundance.tsv` | Long table: `SampleID`, `Function`, `FunctionReadCount`, `PctEligibleProkaryoticReads`, `PctClassifiedReads`, `PctTotalReads` |
| `faprotax_mapping_coverage.tsv` | Exact, mutually reconciling read and taxon coverage per sample |
| `faprotax_taxon_function_assignments.tsv` | Positive taxon-function mappings with stable `FeatureID` and original `TaxonPath` |
| `faprotax_top_functions.png` | Independent bars/facets for the configured top functions; never a 100% stack |
| `faprotax_provenance.json` | microeco version, FAPROTAX version, method, database source, denominator definitions, and overlap warning |

Required coverage columns and invariants:

```text
TotalReads
UnclassifiedReads
ClassifiedReads
EligibleProkaryoticReads
ExcludedNonProkaryoticClassifiedReads
FunctionMappedReads
FunctionUnmappedEligibleReads
EligibleProkaryoticTaxa
FunctionMappedTaxa

ClassifiedReads + UnclassifiedReads == TotalReads
EligibleProkaryoticReads + ExcludedNonProkaryoticClassifiedReads == ClassifiedReads
FunctionMappedReads + FunctionUnmappedEligibleReads == EligibleProkaryoticReads
0 <= FunctionReadCount <= EligibleProkaryoticReads       [for each function]
```

`FunctionMappedReads` counts each read once if its taxon has at least one
FAPROTAX assignment. It is not the sum of `FunctionReadCount`, because one taxon
may contribute to multiple functions.

If no positive-count prokaryotic rows exist, write `faprotax_skipped.tsv` and
return structured status `skipped`. An installed package/API/database-version
mismatch is a failure, not a skip.

## 3. Rough implementation diffs

These blocks specify the intended interfaces and calculations; Luna should
adjust hunk locations to the current files and keep existing style.

### 3.1 Module-aware dependency checks

```diff
diff --git a/analysis/utils/dependencies.R b/analysis/utils/dependencies.R
--- a/analysis/utils/dependencies.R
+++ b/analysis/utils/dependencies.R
@@
 RUNTIME_PACKAGES <- c(
   ...
 )
+MODULE_PACKAGES <- list(faprotax = "microeco")
+MICROECO_MIN_VERSION <- "2.3.0"

+get_module_packages <- function(modules) {
+  unique(unlist(MODULE_PACKAGES[intersect(modules, names(MODULE_PACKAGES))],
+                use.names = FALSE))
+}
+
+check_module_dependencies <- function(modules) {
+  pkgs <- get_module_packages(modules)
+  check_dependencies(pkgs)
+  if ("faprotax" %in% modules &&
+      utils::packageVersion("microeco") < utils::package_version(MICROECO_MIN_VERSION)) {
+    stop(sprintf("Module 'faprotax' requires microeco >= %s.", MICROECO_MIN_VERSION),
+         call. = FALSE)
+  }
+  invisible(TRUE)
+}

-get_required_packages <- function(include_tests = FALSE) {
-  if (isTRUE(include_tests)) REQUIRED_PACKAGES else RUNTIME_PACKAGES
+get_required_packages <- function(include_tests = FALSE, include_modules = FALSE) {
+  ans <- if (isTRUE(include_tests)) REQUIRED_PACKAGES else RUNTIME_PACKAGES
+  if (isTRUE(include_modules)) ans <- unique(c(ans, unlist(MODULE_PACKAGES)))
+  ans
 }
```

```diff
diff --git a/analysis/install_packages.R b/analysis/install_packages.R
--- a/analysis/install_packages.R
+++ b/analysis/install_packages.R
@@
-REQUIRED_PACKAGES <- get_required_packages(include_tests = TRUE)
+REQUIRED_PACKAGES <- get_required_packages(include_tests = TRUE, include_modules = TRUE)
```

### 3.2 Configuration and registration

```diff
diff --git a/analysis/utils/config.R b/analysis/utils/config.R
--- a/analysis/utils/config.R
+++ b/analysis/utils/config.R
@@ validate_config
+  assert_scalar_number(cfg$faprotax$top_n_functions,
+                       "faprotax.top_n_functions", lower = 1, integer = TRUE)
@@ get_default_config
+    faprotax = list(
+      top_n_functions = 20L
+    ),
@@ output directories
     kreport = file.path(base_out, "07_Kreport"),
+    faprotax = file.path(base_out, "08_FAPROTAX")
```

```diff
diff --git a/config.example.yml b/config.example.yml
--- a/config.example.yml
+++ b/config.example.yml
@@
+faprotax:
+  # Plot selection only; all non-zero functions remain in the TSV output.
+  top_n_functions: 20
```

```diff
diff --git a/analysis/00_run_pipeline.R b/analysis/00_run_pipeline.R
--- a/analysis/00_run_pipeline.R
+++ b/analysis/00_run_pipeline.R
@@ source modules
 source(file.path(script_dir, "07_kreport_pavian.R"))
+source(file.path(script_dir, "08_faprotax.R"))
@@ module registry
   kreport     = run_kreport,
+  faprotax    = run_faprotax
 )
@@ after requested-module validation, before --validate-only/output mutation
+check_module_dependencies(requested_modules)
+if ("faprotax" %in% requested_modules) validate_faprotax_runtime()
@@ dependency provenance
-deps <- get_dependency_versions()
+deps <- get_dependency_versions(unique(c(
+  RUNTIME_PACKAGES,
+  get_module_packages(requested_modules)
+)))
```

Do not change the default module list in `load_config()` in this first patch.

### 3.3 New analysis module

```diff
diff --git a/analysis/08_faprotax.R b/analysis/08_faprotax.R
new file mode 100644
--- /dev/null
+++ b/analysis/08_faprotax.R
@@
+FAPROTAX_EXPECTED_VERSION <- "1.2.12"
+
+get_faprotax_database <- function() {
+  env <- new.env(parent = emptyenv())
+  utils::data("prok_func_FAPROTAX", package = "microeco", envir = env)
+  db <- env$prok_func_FAPROTAX
+  if (is.null(db) || is.null(db$ver)) {
+    stop("microeco did not expose embedded FAPROTAX version metadata.", call. = FALSE)
+  }
+  db
+}
+
+validate_faprotax_runtime <- function() {
+  if (!requireNamespace("microeco", quietly = TRUE)) {
+    stop("Module 'faprotax' requires the optional package 'microeco'.", call. = FALSE)
+  }
+  if (utils::packageVersion("microeco") < utils::package_version("2.3.0")) {
+    stop("Module 'faprotax' requires microeco >= 2.3.0.", call. = FALSE)
+  }
+  db <- get_faprotax_database()
+  observed <- as.character(db$ver)
+  if (!identical(observed, FAPROTAX_EXPECTED_VERSION)) {
+    stop(sprintf("Expected FAPROTAX %s, but microeco embeds %s.",
+                 FAPROTAX_EXPECTED_VERSION, observed), call. = FALSE)
+  }
+  list(
+    microeco_version = as.character(utils::packageVersion("microeco")),
+    faprotax_version = observed
+  )
+}
+
+prepare_faprotax_input <- function(context) {
+  classified <- setdiff(seq_len(nrow(context$count_matrix)), context$unclass_index)
+  tax <- context$taxonomy[classified, , drop = FALSE]
+  counts <- context$count_matrix[classified, , drop = FALSE]
+  eligible <- tax$superkingdom %in% c("Bacteria", "Archaea") & rowSums(counts) > 0
+  if (!any(eligible)) return(NULL)
+
+  tax <- tax[eligible, , drop = FALSE]
+  counts <- counts[eligible, , drop = FALSE]
+  ids <- sprintf("Taxon_%06d", classified[eligible])
+  rownames(counts) <- ids
+
+  # microeco expects Kingdom to contain Bacteria/Archaea. Deliberately omit
+  # the intervening NCBI kingdom rank from its collapsed taxonomy input.
+  tax_table <- data.frame(
+    Kingdom = tax$superkingdom,
+    Phylum = tax$phylum,
+    Class = tax$class,
+    Order = tax$order,
+    Family = tax$family,
+    Genus = tax$genus,
+    Species = tax$species,
+    check.names = FALSE,
+    stringsAsFactors = FALSE
+  )
+  tax_table[tax_table %in% c("Unknown", "Unclassified")] <- ""
+  rownames(tax_table) <- ids
+
+  sample_table <- context$metadata
+  if (is.null(sample_table)) {
+    sample_table <- data.frame(SampleID = context$samples,
+                               Group = context$samples,
+                               stringsAsFactors = FALSE)
+    rownames(sample_table) <- sample_table$SampleID
+  }
+  list(counts = counts, tax_table = tax_table, source_taxonomy = tax,
+       sample_table = sample_table, feature_ids = ids)
+}
+
+run_faprotax <- function(context) {
+  runtime <- validate_faprotax_runtime()
+  out_dir <- context$config$output$dirs$faprotax
+  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
+  prepared <- prepare_faprotax_input(context)
+  if (is.null(prepared)) {
+    path <- file.path(out_dir, "faprotax_skipped.tsv")
+    write.table(data.frame(Reason = "No positive-count Bacteria or Archaea rows"),
+                path, sep = "\t", row.names = FALSE, quote = FALSE)
+    return(list(status = "skipped", reason = "No eligible prokaryotic taxa",
+                outputs = path))
+  }
+
+  dataset <- microeco::microtable$new(
+    otu_table = as.data.frame(prepared$counts, check.names = FALSE),
+    sample_table = prepared$sample_table,
+    tax_table = prepared$tax_table,
+    auto_tidy = FALSE
+  )
+  predictor <- microeco::trans_func$new(dataset)
+  if (!identical(predictor$for_what, "prok")) {
+    stop("microeco did not recognize the prepared taxonomy as prokaryotic.", call. = FALSE)
+  }
+  predictor$cal_func(prok_database = "FAPROTAX")
+  binary <- as.matrix(predictor$res_func[prepared$feature_ids, , drop = FALSE])
+
+  if (anyNA(binary) || any(!binary %in% c(0, 1))) {
+    stop("microeco returned a non-binary or missing taxon-function matrix.", call. = FALSE)
+  }
+  function_counts <- t(binary) %*% as.matrix(prepared$counts)
+  mapped_taxa <- rowSums(binary) > 0
+  mapped_reads <- colSums(prepared$counts[mapped_taxa, , drop = FALSE])
+
+  # Build and write the three TSV contracts. Percentages use the explicit
+  # eligible, classified, and total denominators; retain only functions with
+  # a positive count in at least one selected sample.
+  # ... construct function_abundance, coverage, and positive long mappings ...
+  # ... assert all three read-accounting equations before writing ...
+  # ... select top_n_functions by summed FunctionReadCount and save bar/facets ...
+
+  provenance <- c(runtime, list(
+    engine = "microeco::trans_func$cal_func",
+    database = "FAPROTAX",
+    database_source = "microeco embedded data",
+    abundance_model = "sum of classified taxon counts carrying each function",
+    functions_are_nonexclusive = TRUE,
+    interpretation = paste(
+      "Taxon-based functional inference; not evidence of gene presence,",
+      "pathway completeness, expression, or activity."
+    )
+  ))
+  jsonlite::write_json(provenance,
+                       file.path(out_dir, "faprotax_provenance.json"),
+                       pretty = TRUE, auto_unbox = TRUE, null = "null")
+
+  list(status = "completed", outputs = c(/* exact five artifact paths */))
+}
```

Luna must replace the marked construction sections and placeholder output list
with ordinary R code; no ellipses or placeholders should remain in committed
implementation.

### 3.4 Tests and documentation

```diff
diff --git a/tests/testthat/test-config.R b/tests/testthat/test-config.R
--- a/tests/testthat/test-config.R
+++ b/tests/testthat/test-config.R
@@
+  expect_equal(cfg$output$dirs$faprotax,
+               file.path(cfg$output$base_dir, "08_FAPROTAX"))
+  expect_equal(cfg$faprotax$top_n_functions, 20L)
```

```diff
diff --git a/tests/testthat/test-faprotax.R b/tests/testthat/test-faprotax.R
new file mode 100644
--- /dev/null
+++ b/tests/testthat/test-faprotax.R
@@
+# Source dependencies/config/io/plotting and 08_faprotax.R.
+# Do not skip the release gate when microeco is missing.
+
+test_that("embedded FAPROTAX database is exactly 1.2.12", {
+  runtime <- validate_faprotax_runtime()
+  expect_identical(runtime$faprotax_version, "1.2.12")
+})
+
+test_that("FAPROTAX outputs conserve explicit read denominators", {
+  # Run the module on a controlled fixture containing:
+  #   one mapped prokaryote, one unmapped prokaryote,
+  #   one classified eukaryote, and the canonical unclassified row.
+  # Assert the three accounting equations, exact headers, integer counts,
+  # stable FeatureID/TaxonPath mapping, and all five output artifacts.
+})
+
+test_that("function percentages are not treated as a composition", {
+  # Assert each percentage is bounded [0,100], but do not assert that
+  # percentages across overlapping functions sum to 100.
+})
```

Update `README.md` to:

- add `08_faprotax.R` to the architecture tree;
- document `--modules faprotax` and the optional `microeco` dependency;
- reproduce the exact output/denominator definitions above;
- cite Louca et al. (2016), DOI
  [`10.1126/science.aaf4507`](https://doi.org/10.1126/science.aaf4507);
- state explicitly that this is FAPROTAX 1.2.12 via microeco, not FAPROTAX2;
- state that functions overlap and results are ecological inference rather than
  gene/pathway/activity measurements.

Add an `Unreleased` changelog entry only after deciding whether this lands
before 0.2.0 or as 0.3.0. Do not silently bump `VERSION` or `CITATION.cff` while
0.2.0 remains unreleased.

## 4. Acceptance gates

```bash
Rscript analysis/install_packages.R --install
Rscript -e "files <- list.files('analysis', pattern='[.]R$', recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))"
Rscript tests/testthat.R
Rscript analysis/00_run_pipeline.R --config config.yml --modules faprotax --validate-only
Rscript analysis/00_run_pipeline.R --config config.yml --modules faprotax --output-dir output_faprotax_reference
git diff --check
```

Before release, compare the tracked-fixture taxon-function matrix and coverage
against the official FAPROTAX 1.2.12 `collapse_table.py` workflow. An
[open microeco issue](https://github.com/ChiLiubio/microeco/issues/534) reports
differences between the two engines. Any observed difference must be explained
and recorded; do not describe microeco output as byte-for-byte equivalent to
`collapse_table.py` without evidence.

Manual review must also confirm:

1. no unclassified or non-prokaryotic lineage enters the microeco input;
2. the NCBI `superkingdom`/microeco `Kingdom` mapping is exact;
3. all read-accounting equations reconcile per sample;
4. the plot uses independent bars rather than a compositional stack;
5. the manifest records the actual `microeco` package version and the module's
   five returned artifacts; and
6. no producer-output fixture directories or unrelated plans are staged.
