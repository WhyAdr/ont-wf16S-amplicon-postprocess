# ONT wf-16S Amplicon Post-Processing Pipeline: Modernization & Expansion Plan (Phase 1)
**Document:** `ont-wf16s-amplicon-postprocess-xpand1.md`  
**Target Repository:** `WhyAdr/ont-wf16S-amplicon-postprocess`  
**Reference Benchmark:** `ATW_Sesame_Greenhouse-Examples` (BGI Amplicon Workflow)  
**Date:** September 2, 2026  
**Status:** Implementation Blueprint with Concrete Code Diffs  

---

## 1. Executive Summary & Vision

The current repository ([`WhyAdr/ont-wf16S-amplicon-postprocess`](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess)) hosts a specialized downstream analysis toolkit for Oxford Nanopore full-length 16S amplicon data produced by EPI2ME [`wf-16s`](https://github.com/epi2me-labs/wf-16s). It successfully demonstrated high-resolution species profiling, 100-bootstrap diversity subsampling, NCBI TaxID resolution, and interactive Pavian Sankey generation on the `AmbarAyunda_minimap2_16S` bioslurry dataset.

However, its implementation currently resides in a **single monolithic 604-line script** ([`analyze_16s_improved.R`](file:///d:/W/AAy_Amplicon/analyze_16s_improved.R)) coupled with ad-hoc companion scripts, hardcoded paths, and assumptions of a single-sample dataset. 

By contrast, the reference project [`ATW_Sesame_Greenhouse-Examples`](file:///d:/W/AAy_Amplicon/ATW_Sesame_Greenhouse-Examples) exemplifies a **production-grade, highly modular, config-driven bioinformatics pipeline**. It features centralized YAML configuration, clean script isolation via sandboxed environments, automated dependency provisioning, multi-group comparative statistics (PCoA with Procrustes alignment, PERMANOVA, NMDS, UpSet shared taxa), and standardized output hierarchies.

### Strategic Objective
Transform `ont-wf16S-amplicon-postprocess` into a **modular, dual-mode post-processing framework** that combines:
1. **The engineering excellence of the ATW pipeline**: Centralized YAML config, numbered single-responsibility scripts, unified plotting aesthetics, and multi-sample comparative statistics.
2. **The unique biological strengths of Oxford Nanopore full-length 16S**: Full gene resolution (8-rank taxonomy from Superkingdom to Species), per-read length/quality QC diagnostics, NCBI E-utilities lineage resolution, and automated Kraken `.kreport` / Pavian Sankey export.

---

## 2. Comparative Architecture & Target Layout

```
d:/W/AAy_Amplicon/
├── config.yml                         # [NEW] Active pipeline configuration
├── config.example.yml                 # [NEW] Documented configuration template
├── metadata.example.tsv               # [NEW] Example cohort metadata sheet
├── README.md                          # Repository overview & quickstart
├── ont-wf16s-amplicon-postprocess-xpand1.md  # This blueprint specification
├── .gitignore                         # Excludes heavy BAMs, fastqs, temp files
│
├── analysis/                          # [NEW] Modular numbered pipeline scripts
│   ├── install_packages.R             # [NEW] Automated CRAN / Bioconductor installer
│   ├── 00_run_pipeline.R             # [NEW] Master orchestrator (runs 01->07)
│   ├── 01_qc_diagnostics.R           # [NEW] Nanopore read QC & length vs classification
│   ├── 02_alpha_diversity.R          # [NEW] Richness, diversity, bootstrap & rarefaction
│   ├── 03_beta_diversity.R           # [NEW] Bray-Curtis/Jaccard, PCoA + Procrustes, PERMANOVA
│   ├── 04_taxa_composition.R         # [NEW] Phylum-to-species stacked barplots & heatmaps
│   ├── 05_ordination.R               # [NEW] Hellinger PCA & NMDS ordination
│   ├── 06_shared_taxa.R              # [NEW] UpSet & Venn core/pan microbiome
│   ├── 07_kreport_pavian.R           # [NEW] Kraken .kreport & Pavian interactive Sankey
│   │
│   └── utils/                         # [NEW] Shared utility functions
│       ├── load_config.R              # [NEW] YAML loader with CLI override support
│       ├── plotting_theme.R           # [NEW] theme_amplicon & color palettes
│       └── ncbi_taxonomy.py           # [NEW] E-utilities TaxID resolver
│
├── output_AAy/                        # Existing EPI2ME wf-16s output files (input data)
│   ├── abundance_table_species.tsv
│   ├── params.json
│   ├── taxonomy_cache.json
│   ├── versions_all.txt
│   └── reads_assignments/
│
└── output/                            # [NEW] Standardized pipeline output root
    ├── 01_QC/                         # Donut, read length distribution, violin QC
    ├── 02_Alpha_Diversity/            # Boxplots, rarefaction curves, bootstrap CSVs
    ├── 03_Beta_Diversity/             # PCoA plots, Bray-Curtis distance tables
    ├── 04_Taxa_Composition/          # Stacked barplots, top-N tables, heatmaps
    ├── 05_Ordination/                 # PCA biplots, NMDS ordination
    ├── 06_Shared_Taxa/                # UpSet plots, shared species tables
    └── 07_Kreport_Sankey/             # .kreport output & standalone HTML Sankey
```

---

## 3. Concrete Code Implementations & Diffs

### 3.1 Configuration System

#### [NEW] `config.yml` (and `config.example.yml`)
Centralizes all inputs, parameters, and outputs. Paths are resolved relative to the workspace root.

```yaml
# ==============================================================================
# config.yml - ONT wf-16S Amplicon Post-Processing Pipeline Configuration
# ==============================================================================

project_name: "AmbarAyunda_16S_Amplicon"

# Execution mode: "auto" (detects single sample vs cohort), "single", or "cohort"
mode: "auto"

input:
  abundance_table:  "output_AAy/abundance_table_species.tsv"
  assignments_file: "output_AAy/reads_assignments/AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv"
  metadata:         null       # Optional: "metadata.tsv" for multi-sample cohorts
  taxonomy_cache:   "output_AAy/taxonomy_cache.json"
  exclude_columns:  ["total"]  # wf-16s artifact columns to exclude from sample detection

output:
  base_dir:         "output"
  qc:               "output/01_QC"
  alpha:            "output/02_Alpha_Diversity"
  beta:             "output/03_Beta_Diversity"
  composition:      "output/04_Taxa_Composition"
  ordination:       "output/05_Ordination"
  shared:           "output/06_Shared_Taxa"
  kreport:          "output/07_Kreport_Sankey"

parameters:
  # QC gating parameters (bp)
  min_read_length: 1200
  max_read_length: 1800
  target_min_len:  1400
  target_max_len:  1600

  # Alpha diversity & bootstrap rarefaction
  bootstrap_subsample_depth: 70000
  bootstrap_iterations: 100
  rarefaction_step: 2500

  # Composition settings
  top_n_taxa: 15
  ranks: ["superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species"]

  # Beta diversity & ordination
  pcoa_iterations: 100
  pcoa_subsample_fraction: 0.75

  # Export settings
  dpi: 150
```

---

#### [NEW] `config.example.yml`
A fully commented copy of `config.yml` with placeholder values. Ship as the onboarding template:
- All `input:` paths set to `"<your_wf16s_output>/..."` placeholders
- All `parameters:` retain default values with inline comments explaining each
- `mode:` set to `"auto"` with comments documenting `"single"` and `"cohort"` triggers

---

#### [NEW] `metadata.example.tsv`
Defines the contract for cohort-mode execution. Tab-delimited with required and optional columns:

```tsv
SampleID	Group	Condition
AmbarAyunda_minimap2_16S	Bioslurry	Anaerobic
SampleB_minimap2_16S	Control	Aerobic
SampleC_minimap2_16S	Bioslurry	Aerobic
```

| Column | Required | Description |
|--------|----------|-------------|
| `SampleID` | ✅ | Must exactly match a column name in `abundance_table_species.tsv` |
| `Group` | ✅ | Grouping variable for beta diversity, PERMANOVA, and UpSet analysis |
| `Condition` | ❌ | Optional covariate for stratified analyses |

> **Note:** `SampleID` values are matched against `colnames(abundance_table)` after excluding `tax` and any columns listed in `config.yml → input.exclude_columns`.

---

#### [NEW] `analysis/utils/load_config.R`
Ported and adapted from [`ATW_Sesame_Greenhouse-Examples/analysis/utils/load_config.R`](file:///d:/W/AAy_Amplicon/ATW_Sesame_Greenhouse-Examples/analysis/utils/load_config.R).

```r
# ==============================================================================
# analysis/utils/load_config.R
# Loads pipeline configuration and resolves relative paths
# ==============================================================================

load_config <- function(config_path = NULL) {
  library(yaml)

  # CLI override priority: explicit arg > --config flag > default "../config.yml"
  if (is.null(config_path)) {
    cli_args <- commandArgs(trailingOnly = TRUE)
    config_idx <- which(cli_args == "--config")
    if (length(config_idx) > 0 && config_idx < length(cli_args)) {
      config_path <- cli_args[config_idx + 1]
    } else {
      config_path <- "../config.yml"
    }
  }

  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s\nPlease create config.yml or copy config.example.yml", config_path))
  }

  cfg <- yaml::read_yaml(config_path)
  cfg_root <- dirname(normalizePath(config_path, mustWork = TRUE))

  resolve <- function(p) {
    if (is.null(p) || is.na(p) || p == "") return(NULL)
    normalizePath(file.path(cfg_root, p), mustWork = FALSE)
  }

  # Resolve path-valued entries only (skip list values like exclude_columns)
  path_keys <- c("abundance_table", "assignments_file", "metadata", "taxonomy_cache")
  for (k in path_keys) {
    if (k %in% names(cfg$input)) cfg$input[[k]] <- resolve(cfg$input[[k]])
  }
  for (k in names(cfg$output)) cfg$output[[k]] <- resolve(cfg$output[[k]])

  # CLI overrides
  cli_args <- commandArgs(trailingOnly = TRUE)
  out_idx <- which(cli_args == "--output-dir")
  if (length(out_idx) > 0 && out_idx < length(cli_args)) {
    cfg$output$base_dir <- normalizePath(cli_args[out_idx + 1], mustWork = FALSE)
  }

  cfg
}

# --- Canonical sample/group ordering helper ---
# Returns groups in config-defined display order, filtered to current metadata.
get_group_order <- function(cfg, meta_groups) {
  if (!is.null(cfg$parameters$group_order)) {
    return(intersect(cfg$parameters$group_order, meta_groups))
  }
  sort(unique(meta_groups))
}
```

---

#### [NEW] `analysis/utils/plotting_theme.R`
Encapsulates `theme_amplicon()` and high-contrast colorblind-friendly taxonomic palettes.

```r
# ==============================================================================
# analysis/utils/plotting_theme.R
# Unified design system and colorblind-safe palettes for ONT amplicon plots
# ==============================================================================

suppressMessages({
  library(ggplot2)
  library(scales)
})

theme_amplicon <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
    plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b = 10)),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11)
  )

phylum_palette <- c(
  "Bacillota"                  = "#1b9e77",
  "Pseudomonadota"             = "#d95f02",
  "Bacteroidota"               = "#7570b3",
  "Actinomycetota"             = "#e7298a",
  "Fusobacteriota"             = "#66a61e",
  "Spirochaetota"              = "#e6ab02",
  "Synergistota"               = "#a6761d",
  "Planctomycetota"            = "#4e79a7",
  "Campylobacterota"           = "#e15759",
  "Chloroflexota"              = "#76b7b2",
  "Candidatus Melainabacteria"  = "#59a14f",
  "Ignavibacteriota"           = "#edc949",
  "Armatimonadota"             = "#af7aa1",
  "Verrucomicrobiota"          = "#666666",
  "Unclassified"               = "#bdbdbd",
  "Other phyla"                = "#41566b"
)

get_phylum_colors <- function(phyla) {
  cols <- phylum_palette[phyla]
  cols[is.na(cols)] <- "#a8a8a8"
  names(cols) <- phyla
  cols
}
```

---

#### [NEW] `analysis/install_packages.R`
Modeled directly after `ATW_Sesame_Greenhouse-Examples/analysis/install_packages.R`, verifying all requirements for ONT analysis.

```r
# ==============================================================================
# analysis/install_packages.R
# Installs required CRAN and Bioconductor packages for ONT wf-16s postprocessing
# ==============================================================================

cat("=== Checking and Installing Dependencies for ONT wf-16S Post-Processing ===\n")

cran_pkgs <- c(
  "yaml", "dplyr", "tidyr", "stringr", "ggplot2", 
  "scales", "vegan", "RColorBrewer", "jsonlite", "pheatmap", 
  "UpSetR", "VennDiagram", "ggrepel"
)

cran_ok <- 0; cran_fail <- c()
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing CRAN: %s ...", pkg))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf(" OK (v%s)\n", packageVersion(pkg)))
        cran_ok <- cran_ok + 1
      } else { cat(" FAILED\n"); cran_fail <- c(cran_fail, pkg) }
    }, error = function(e) { cat(sprintf(" ERROR: %s\n", e$message)); cran_fail <<- c(cran_fail, pkg) })
  } else {
    cat(sprintf("  Already installed: %s (v%s)\n", pkg, packageVersion(pkg)))
    cran_ok <- cran_ok + 1
  }
}

cat(sprintf("\nCRAN status: %d/%d available.\n", cran_ok, length(cran_pkgs)))
if (length(cran_fail) == 0) cat("SUCCESS: All required packages are ready.\n")
```

---

### 3.2 Master Orchestrator

#### [NEW] `analysis/00_run_pipeline.R`
Executes modules in isolated sandboxes (`new.env()`), passing configuration without global namespace collisions.

```r
# ==============================================================================
# analysis/00_run_pipeline.R
# Master Orchestrator for ONT wf-16s Post-Processing
# ==============================================================================

cat("==================================================================\n")
cat(" ONT wf-16S Amplicon Post-Processing Pipeline\n")
cat("==================================================================\n\n")

source("analysis/utils/load_config.R")
cfg <- load_config()

# Load abundance table to detect sample count and execution mode
ab_raw <- read.delim(cfg$input$abundance_table, header = TRUE, sep = "\t", check.names = FALSE)
# Exclude non-sample columns: 'tax' plus any wf-16s artifact columns (e.g., 'total')
extra_exclude <- if (!is.null(cfg$input$exclude_columns)) cfg$input$exclude_columns else character(0)
sample_cols <- setdiff(colnames(ab_raw), c("tax", extra_exclude))
n_samples <- length(sample_cols)

has_metadata <- !is.null(cfg$input$metadata) && file.exists(cfg$input$metadata)
is_cohort <- (cfg$mode == "cohort") || (cfg$mode == "auto" && (n_samples > 1 || has_metadata))

cat(sprintf("[Mode] Detected %d sample(s). Execution Mode: %s\n", 
            n_samples, if (is_cohort) "COHORT (Multi-Sample)" else "SINGLE-SAMPLE DEEP DIVE"))

modules <- c(
  "01_qc_diagnostics.R",
  "02_alpha_diversity.R",
  "03_beta_diversity.R",
  "04_taxa_composition.R",
  "05_ordination.R",
  "06_shared_taxa.R",
  "07_kreport_pavian.R"
)

for (script_name in modules) {
  script_path <- file.path("analysis", script_name)
  if (!file.exists(script_path)) {
    cat(sprintf("[SKIP] Module %s not found.\n", script_name))
    next
  }

  cat(sprintf("\n>>> Executing [%s] ...\n", script_name))
  tryCatch({
    env <- new.env(parent = globalenv())
    env$cfg <- cfg
    env$is_cohort <- is_cohort
    env$sample_cols <- sample_cols
    source(script_path, local = env)
    cat(sprintf(">>> [%s] Completed successfully.\n", script_name))
  }, error = function(e) {
    cat(sprintf(">>> [%s] ERROR: %s\n", script_name, e$message))
  })
}

cat("\n==================================================================\n")
cat(" Pipeline Run Finished. Outputs saved in: ", cfg$output$base_dir, "\n")
cat("==================================================================\n")
```

---

### 3.3 Decomposing `analyze_16s_improved.R` into Standalone Modules

#### Migration Diff: QC Diagnostics (`analyze_16s_improved.R` Lines 92–186 $\rightarrow$ `analysis/01_qc_diagnostics.R`)

```diff
- # [IN analyze_16s_improved.R (monolithic)]
- reads <- read.delim("output_AAy/reads_assignments/AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv", ...)
- ggsave("figs/01a_classification_donut.png", ...)
- ggsave("figs/01b_read_length_distribution.png", ...)
- ggsave("figs/01c_qc_filter_diagnostic.png", ...)
- ggsave("figs/06_length_vs_classification.png", ...)

+ # [NEW: analysis/01_qc_diagnostics.R]
+ source("analysis/utils/load_config.R")
+ source("analysis/utils/plotting_theme.R")
+ if (!exists("cfg")) cfg <- load_config()
+ 
+ if (!file.exists(cfg$input$assignments_file)) {
+   cat("[QC] Assignments file not provided or not found. Skipping per-read QC.\n")
+ } else {
+   reads <- read.delim(cfg$input$assignments_file, header = FALSE, sep = "\t",
+                       col.names = c("status","read_id","taxid","len_field","lineage"),
+                       stringsAsFactors = FALSE) %>%
+     mutate(
+       read_length = as.integer(sub(".*\\|([0-9]+).*", "\\1", len_field)),
+       truly_classified = (status == "C") & (taxid != 0),
+       status_label_corrected = ifelse(truly_classified, "Classified", "Unclassified")
+     )
+   # Export Donut, Length Histograms, QC Filter Diagnostic, and Violin QC to cfg$output$qc
+   ggsave(file.path(cfg$output$qc, "01a_classification_donut.png"), p1a, width = 5.5, height = 5.5, dpi = cfg$parameters$dpi)
+   ggsave(file.path(cfg$output$qc, "01b_read_length_distribution.png"), p1b, width = 7.5, height = 5, dpi = cfg$parameters$dpi)
+   ggsave(file.path(cfg$output$qc, "01c_qc_filter_diagnostic.png"), p1c, width = 7, height = 5.5, dpi = cfg$parameters$dpi)
+   ggsave(file.path(cfg$output$qc, "01d_length_vs_classification.png"), p6, width = 6.5, height = 5.5, dpi = cfg$parameters$dpi)
+   cat("[QC] Diagnostic figures written to:", cfg$output$qc, "\n")
+ }
```

---

#### Migration Diff: Alpha Diversity & Bootstrap (`analyze_16s_improved.R` Lines 397–550 $\rightarrow$ `analysis/02_alpha_diversity.R`)

```diff
- # [IN analyze_16s_improved.R]
- # Hardcoded sample column, manual loop, saving directly to figs/ and tables/
- subsample_depth <- 70000
- boot_results <- lapply(1:100, function(i) { ... })
- ggsave("figs/05_rarefaction_curve.png", ...)
- ggsave("figs/05b_diversity_boxplots.png", ...)
- write.csv(boot_df, "tables/bootstrap_diversity.csv", ...)

+ # [NEW: analysis/02_alpha_diversity.R]
+ source("analysis/utils/load_config.R")
+ source("analysis/utils/plotting_theme.R")
+ if (!exists("cfg")) cfg <- load_config()
+ 
+ # Computes Sobs, Chao1, Shannon, Simpson, Pielou, and Berger-Parker
+ # In single-sample mode:
+ #   Runs 100-iteration bootstrap subsampling at cfg$parameters$bootstrap_subsample_depth
+ #   Generates rarefaction curves and bootstrap boxplots
+ # In cohort mode:
+ #   Merges with metadata and generates group-level boxplots with Kruskal-Wallis/Wilcoxon tests
+ write.csv(alpha_summary, file.path(cfg$output$alpha, "alpha_diversity_summary.csv"), row.names = FALSE)
+ ggsave(file.path(cfg$output$alpha, "02a_diversity_boxplots.png"), p_div, width = 8, height = 6.5, dpi = cfg$parameters$dpi)
+ ggsave(file.path(cfg$output$alpha, "02b_rarefaction_curve.png"), p_rare, width = 7, height = 5, dpi = cfg$parameters$dpi)
```

---

#### Migration Diff: Composition Analysis (`analyze_16s_improved.R` Lines 194–438 $\rightarrow$ `analysis/04_taxa_composition.R`)

Migrates all phylum, family, genus, and species composition figures and tables — the largest single block (~250 lines) of the monolith.

```diff
- # [IN analyze_16s_improved.R (monolithic)]
- # Hardcoded palettes, manual factor reordering, ad-hoc N_KEEP constants
- phylum_comp <- ab %>% ... %>% arrange(desc(count))
- ggsave("figs/02_phylum_composition.png", ...)
- ggsave("figs/02c_phylum_stacked.png", ...)
- ggsave("figs/02b_family_composition.png", ...)
- ggsave("figs/02d_family_stacked.png", ...)
- ggsave("figs/03_top20_genera.png", ...)
- ggsave("figs/03b_genus_stacked.png", ...)
- ggsave("figs/04_top15_species.png", ...)
- write.csv(phylum_comp_full, "tables/phylum_composition.csv", ...)
- write.csv(family_comp, "tables/family_composition.csv", ...)
- write.csv(genus_comp, "tables/genus_composition.csv", ...)

+ # [NEW: analysis/04_taxa_composition.R]
+ source("analysis/utils/load_config.R")
+ source("analysis/utils/plotting_theme.R")
+ if (!exists("cfg")) cfg <- load_config()
+ 
+ ab_raw <- read.delim(cfg$input$abundance_table, header = TRUE, sep = "\t", check.names = FALSE)
+ ab <- ab_raw %>%
+   tidyr::separate(tax, into = cfg$parameters$ranks, sep = ";", fill = "right", remove = FALSE) %>%
+   dplyr::mutate(
+     is_unclassified = superkingdom == "Unclassified",
+     rel_abund = count / sum(count[!is_unclassified])
+   )
+ 
+ top_n <- cfg$parameters$top_n_taxa
+ out_dir <- cfg$output$composition
+ dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
+ 
+ # --- Phylum bar + stacked ---
+ # Reuses get_phylum_colors() from plotting_theme.R
+ # Groups low-abundance phyla into "Other phyla" (top 7 retained)
+ ggsave(file.path(out_dir, "04a_phylum_bar.png"), p_phy_bar, width = 7.5, height = 5, dpi = cfg$parameters$dpi)
+ ggsave(file.path(out_dir, "04b_phylum_stacked.png"), p_phy_stack, width = 7, height = 7, dpi = cfg$parameters$dpi)
+ 
+ # --- Family bar + stacked (top N, colored by phylum) ---
+ ggsave(file.path(out_dir, "04c_family_bar.png"), p_fam_bar, width = 9, height = 6, dpi = cfg$parameters$dpi)
+ ggsave(file.path(out_dir, "04d_family_stacked.png"), p_fam_stack, width = 7.5, height = 7, dpi = cfg$parameters$dpi)
+ 
+ # --- Genus bar + stacked (top 20, colored by phylum) ---
+ ggsave(file.path(out_dir, "04e_genus_bar.png"), p_gen_bar, width = 9, height = 7, dpi = cfg$parameters$dpi)
+ ggsave(file.path(out_dir, "04f_genus_stacked.png"), p_gen_stack, width = 7.5, height = 7, dpi = cfg$parameters$dpi)
+ 
+ # --- Species bar (top N, italic labels, colored by phylum) ---
+ ggsave(file.path(out_dir, "04g_top_species.png"), p_species, width = 9.5, height = 6.5, dpi = cfg$parameters$dpi)
+ 
+ # --- Heatmap (cohort mode only — shows top taxa × sample abundance) ---
+ if (exists("is_cohort") && is_cohort) {
+   pheatmap::pheatmap(top_taxa_matrix, filename = file.path(out_dir, "04h_heatmap.png"), ...)
+ }
+ 
+ # --- Export composition tables ---
+ write.csv(phylum_comp, file.path(out_dir, "phylum_composition.csv"), row.names = FALSE)
+ write.csv(family_comp, file.path(out_dir, "family_composition.csv"), row.names = FALSE)
+ write.csv(genus_comp, file.path(out_dir, "genus_composition.csv"), row.names = FALSE)
+ cat("[Composition] All figures and tables written to:", out_dir, "\n")
```

> **Figure renumbering note:** Monolith figures `02_`, `02b–d_`, `03_`, `03b_`, `04_` are consolidated under the `04_` prefix in the new module. Old→new mapping: `02_phylum_composition` → `04a_phylum_bar`, `02c_phylum_stacked` → `04b_phylum_stacked`, `02b_family_composition` → `04c_family_bar`, `02d_family_stacked` → `04d_family_stacked`, `03_top20_genera` → `04e_genus_bar`, `03b_genus_stacked` → `04f_genus_stacked`, `04_top15_species` → `04g_top_species`.

---

### 3.4 New Analytical Capabilities (Cohort Mode)

#### [NEW] `analysis/03_beta_diversity.R`
Integrates ATW's advanced **100-iteration Procrustes-aligned consensus PCoA** and PERMANOVA.

```r
# ==============================================================================
# analysis/03_beta_diversity.R
# Bray-Curtis / Jaccard Beta Diversity & Consensus Procrustes PCoA
# ==============================================================================

source("analysis/utils/load_config.R")
source("analysis/utils/plotting_theme.R")
if (!exists("cfg")) cfg <- load_config()

if (!exists("is_cohort") || !is_cohort) {
  cat("[Beta] Single-sample mode active. Beta diversity requires >= 3 samples. Skipping.\n")
} else {
  # 1. Load abundance and metadata
  # 2. Compute Bray-Curtis distance matrix
  # 3. Perform 100-iteration bootstrapped PCoA with Procrustes consensus alignment
  # 4. Run adonis2 PERMANOVA testing group separation
  # 5. Output consensus PCoA plot to cfg$output$beta
}
```

---

#### [NEW] `analysis/05_ordination.R`
Implements Hellinger PCA biplots (works in single-sample mode for species decomposition) and NMDS (cohort mode). Inspired by ATW's [`08_pca_analysis.R`](file:///d:/W/AAy_Amplicon/ATW_Sesame_Greenhouse-Examples/analysis/08_pca_analysis.R) and [`18_nmds.R`](file:///d:/W/AAy_Amplicon/ATW_Sesame_Greenhouse-Examples/analysis/18_nmds.R).

```r
# ==============================================================================
# analysis/05_ordination.R
# Hellinger PCA & NMDS Ordination
# ==============================================================================

source("analysis/utils/load_config.R")
source("analysis/utils/plotting_theme.R")
if (!exists("cfg")) cfg <- load_config()

out_dir <- cfg$output$ordination
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Hellinger PCA (single-sample or cohort) ---
# 1. Load abundance matrix and apply Hellinger transformation: sqrt(p_ij)
# 2. Run PCA via prcomp() on Hellinger-transformed community matrix
# 3. Generate biplot showing species loadings in PC1 × PC2 space
#    - Labels top-contributing species using ggrepel
#    - Annotates axes with % variance explained
ggsave(file.path(out_dir, "05a_hellinger_pca.png"), p_pca, width = 8, height = 7, dpi = cfg$parameters$dpi)

# --- NMDS (cohort mode only, requires >= 3 samples) ---
if (exists("is_cohort") && is_cohort) {
  # 1. Compute Bray-Curtis distance matrix from vegan::vegdist()
  # 2. Run metaMDS() with k=2, trymax=200, autotransform=FALSE
  # 3. Plot NMDS ordination with stress annotation and group ellipses
  # 4. Overlay envfit() vectors for metadata covariates
  ggsave(file.path(out_dir, "05b_nmds.png"), p_nmds, width = 8, height = 7, dpi = cfg$parameters$dpi)
  cat("[Ordination] NMDS ordination written to:", out_dir, "\n")
} else {
  cat("[Ordination] Single-sample mode: NMDS skipped (requires >= 3 samples).\n")
}
```

---

#### [NEW] `analysis/06_shared_taxa.R`
Integrates ATW's UpSetR and shared taxa analysis for multi-group comparisons.

```r
# ==============================================================================
# analysis/06_shared_taxa.R
# Shared and unique taxa analysis across groups (UpSet & Venn)
# ==============================================================================

source("analysis/utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

if (!exists("is_cohort") || !is_cohort) {
  cat("[Shared] Single-sample mode active. Shared taxa analysis requires >= 2 groups. Skipping.\n")
} else {
  library(UpSetR)
  # 1. Identify species presence/absence per group
  # 2. Export shared/unique species matrix to cfg$output$shared/shared_species_matrix.tsv
  # 3. Render UpSet plot to cfg$output$shared/UpSet_Shared_Species.png
}
```

---

### 3.5 Kraken & Pavian Sankey Integration

#### Migration Diff: `fetch_ncbi_taxonomy.py` + `convert_to_kreport.R` $\rightarrow$ `analysis/07_kreport_pavian.R`
Unifies the two separate scripts into a seamless, automated workflow.

```diff
- # Old workflow required manually running:
- # 1. python fetch_ncbi_taxonomy.py
- # 2. Rscript convert_to_kreport.R

+ # [NEW: analysis/07_kreport_pavian.R]
+ source("analysis/utils/load_config.R")
+ if (!exists("cfg")) cfg <- load_config()
+ 
+ cache_file <- cfg$input$taxonomy_cache
+ if (!file.exists(cache_file)) {
+   cat("[kreport] Taxonomy cache not found. Resolving NCBI TaxIDs via Python E-utilities...\n")
+   py_script <- file.path("analysis", "utils", "ncbi_taxonomy.py")
+   system2("python", args = c(py_script,
+     "--abundance", cfg$input$abundance_table,
+     "--assignments", cfg$input$assignments_file,
+     "--cache", cache_file))
+ }
+ 
+ # Build DFS abundance-sorted tree and write .kreport directly to cfg$output$kreport
+ output_file <- file.path(cfg$output$kreport, sprintf("%s.kreport", sample_col))
+ cat(sprintf("[kreport] Successfully written Kraken report to: %s\n", output_file))
```

#### Migration: `fetch_ncbi_taxonomy.py` $\rightarrow$ `analysis/utils/ncbi_taxonomy.py`

The existing [`fetch_ncbi_taxonomy.py`](file:///d:/W/AAy_Amplicon/fetch_ncbi_taxonomy.py) (181 lines) is **moved** to `analysis/utils/ncbi_taxonomy.py` with the following changes:

```diff
- # Hardcoded paths at module level
- ABUNDANCE_FILE = "output_AAy/abundance_table_species.tsv"
- ASSIGNMENTS_FILE = "output_AAy/reads_assignments/AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv"
- CACHE_FILE = "output_AAy/taxonomy_cache.json"

+ import argparse
+ 
+ parser = argparse.ArgumentParser(description="Resolve NCBI TaxIDs from wf-16s assignments")
+ parser.add_argument("--abundance", required=True, help="Path to abundance_table_species.tsv")
+ parser.add_argument("--assignments", required=True, help="Path to read assignments TSV")
+ parser.add_argument("--cache", required=True, help="Output path for taxonomy_cache.json")
+ args = parser.parse_args()
+ 
+ ABUNDANCE_FILE = args.abundance
+ ASSIGNMENTS_FILE = args.assignments
+ CACHE_FILE = args.cache
```

All remaining logic (E-utilities batching, species-to-TaxID resolution via `Counter`, lineage caching) is preserved unchanged.

---

## 4. Execution & Validation Plan

### Verification Step 1: Dry-Run Dependency Check
```bash
Rscript analysis/install_packages.R
```
Verifies all 13 CRAN packages.

### Verification Step 2: Single-Sample Full Pipeline Execution
```bash
Rscript analysis/00_run_pipeline.R --config config.yml
```
- Validates that `output/01_QC/`, `output/02_Alpha_Diversity/`, `output/04_Taxa_Composition/`, and `output/07_Kreport_Sankey/` are created and populated.
- Validates that cohort-only modules (`03_beta_diversity.R`, `05_ordination.R`, `06_shared_taxa.R`) log graceful notices without crashing.

### Verification Step 3: Git Commit & Sync
Once validated, stage and commit the new architecture to `main` on GitHub:
```bash
git add config.yml config.example.yml metadata.example.tsv analysis/ .gitignore
git commit -m "feat: migrate to modular config-driven pipeline inspired by ATW benchmark"
git push origin main
```

#### `.gitignore` Update
Append the following entries to exclude generated outputs:
```gitignore
# Pipeline outputs (regenerated by 00_run_pipeline.R)
output/
*.kreport

# Legacy output dirs (superseded by output/)
figs/
tables/
```
