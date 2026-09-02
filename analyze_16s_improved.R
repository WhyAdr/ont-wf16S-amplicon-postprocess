#!/usr/bin/env Rscript
# =============================================================================
# wf-16s amplicon results - downstream R analysis (Improved)
# Sample: AmbarAyunda_minimap2_16S  (bioslurry sample, minimap2 classifier)
# =============================================================================

# Suppress messages while loading only installed packages
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(vegan)
  library(RColorBrewer)
})

# Use the current directory as working directory to ensure portability.
# Create output subdirectories relative to the workspace root.
dir.create("figs", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

RANKS <- c("superkingdom","clade","phylum","class","order","family","genus","species")

# A premium, colorblind-friendly expanded palette for phyla that covers all 
# phyla in the top-abundance slots (no gray fallback for Planctomycetota or Campylobacterota)
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
  cols[is.na(cols)] <- "#a8a8a8"  # fallback grey for anything not in palette
  names(cols) <- phyla
  cols
}

theme_amplicon <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
    plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b = 10)),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# =============================================================================
# 1. LOAD ABUNDANCE TABLE (species-level)
# =============================================================================
# Using base R read.delim instead of read_tsv for robustness against missing 'readr'
ab_raw <- read.delim("output_AAy/abundance_table_species.tsv", 
                     header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
sample_col <- colnames(ab_raw)[2]   # "AmbarAyunda_minimap2_16S"
cat(sprintf("[1] Loaded abundance table: %d taxa, sample column = '%s'\n",
            nrow(ab_raw), sample_col))

ab <- ab_raw %>%
  rename(count = !!sample_col) %>%
  separate(tax, into = RANKS, sep = ";", fill = "right", remove = FALSE) %>%
  mutate(
    is_unclassified = superkingdom == "Unclassified",
    # rel_abund_total: proportion of ALL reads (including unclassified) — for QC plots
    rel_abund_total = count / sum(count),
    # rel_abund: proportion of CLASSIFIED reads only — for composition analysis
    rel_abund = count / sum(count[!is_unclassified])
  )

total_reads <- sum(ab$count)
n_unclassified <- sum(ab$count[ab$is_unclassified])
n_classified <- total_reads - n_unclassified

cat(sprintf("    Total reads: %d | Classified: %d (%.1f%%) | Unclassified: %d (%.1f%%)\n",
            total_reads, n_classified, 100*n_classified/total_reads,
            n_unclassified, 100*n_unclassified/total_reads))

# =============================================================================
# 2. LOAD PER-READ ASSIGNMENTS
# =============================================================================
# Corrected the input filename and location to match the workspace
reads <- read.delim(
  "output_AAy/reads_assignments/AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv",
  header = FALSE, sep = "\t", 
  col.names = c("status","read_id","taxid","len_field","lineage"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    read_length = as.integer(str_extract(len_field, "(?<=\\|)\\d+")),
    # Clean logic to account for post-alignment filters in Nanopore pipelines
    truly_classified = (status == "C") & (taxid != 0),
    status_label = ifelse(status == "C", "Classified", "Unclassified"),
    status_label_corrected = ifelse(truly_classified, "Classified", "Unclassified")
  )

n_status_C <- sum(reads$status == "C")
n_qc_failed <- sum(reads$status == "C" & reads$taxid == 0)
n_true_classified <- sum(reads$truly_classified)

cat(sprintf("[2] Loaded %d per-read assignments\n", nrow(reads)))
cat(sprintf("    Raw status=='C' (aligned at all):        %d\n", n_status_C))
cat(sprintf("    ...of which failed QC filter (taxid=0):   %d  <-- reclassified by coverage/identity filter\n", n_qc_failed))
cat(sprintf("    TRUE classified (status=C AND taxid!=0):  %d  <-- matches abundance_table_species.tsv\n", n_true_classified))
cat(sprintf("    Raw status=='U':                          %d\n", sum(reads$status=="U")))

# =============================================================================
# FIGURE 1 — Classification donut + read length distribution (combined QC panel)
# =============================================================================
donut_df <- reads %>%
  count(status_label_corrected) %>%
  mutate(
    frac = n / sum(n),
    ymax = cumsum(frac),
    ymin = c(0, head(ymax, -1)),
    label = sprintf("%s\n%s reads\n(%.1f%%)", status_label_corrected, comma(n), 100*frac)
  )

p1a <- ggplot(donut_df, aes(ymin = ymin, ymax = ymax, xmin = 3, xmax = 4, fill = status_label_corrected)) +
  geom_rect(color = "white", linewidth = 1.2) +
  coord_polar(theta = "y") +
  xlim(c(1, 4)) +
  scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
  annotate("text", x = 1, y = 0, label = sprintf("%s\nreads", comma(sum(donut_df$n))),
           size = 4.2, fontface = "bold") +
  theme_void() +
  labs(title = "Read Classification Outcome", fill = NULL) +
  theme(legend.position = "bottom", plot.title = element_text(face="bold", hjust=0.5, size=13)) +
  geom_text(aes(x = 3.5, y = (ymin+ymax)/2, label = label), inherit.aes = FALSE,
            data = donut_df, color = "black", size = 3.3, fontface = "bold")

p1b <- ggplot(reads, aes(x = read_length, fill = status_label_corrected)) +
  geom_histogram(binwidth = 10, alpha = 0.85, position = "identity") +
  scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
  # Using coord_cartesian for zooming is safer, though limits match the exact data range here
  coord_cartesian(xlim = c(1300, 1750)) +
  scale_y_continuous(labels = comma) +
  labs(title = "Read Length Distribution",
       subtitle = "ONT 16S amplicon - expected full-length peak ~1400-1550 bp",
       x = "Read length (bp)", y = "Read count", fill = NULL) +
  theme_amplicon +
  theme(legend.position = "bottom")

ggsave("figs/01a_classification_donut.png", p1a, width = 5.5, height = 5.5, dpi = 150)
ggsave("figs/01b_read_length_distribution.png", p1b, width = 7.5, height = 5, dpi = 150)
cat("[Figure 1] Classification donut + read length histogram saved.\n")

# =============================================================================
# FIGURE 1c — QC diagnostic: why raw status != true classification
# =============================================================================
qc_breakdown <- tibble(
  category = c("True classified\n(aligned + passed QC)",
               "QC-filtered\n(aligned but failed\ncoverage/identity)",
               "Never aligned\n(status = U)"),
  n = c(n_true_classified, n_qc_failed, sum(reads$status == "U"))
)
# Replaced 'fct_reorder' from missing 'forcats' package with base R 'reorder'
qc_breakdown$category <- reorder(factor(qc_breakdown$category), qc_breakdown$n)

p1c <- ggplot(qc_breakdown, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = comma(n)), vjust = -0.4, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = c(
    "True classified\n(aligned + passed QC)" = "#1b9e77",
    "QC-filtered\n(aligned but failed\ncoverage/identity)" = "#e6ab02",
    "Never aligned\n(status = U)" = "#bdbdbd"
  )) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Why Raw 'status' Overstates True Classification Rate",
       subtitle = sprintf(
         "min_ref_coverage / min_percent_identity filter reclassifies %s reads to taxid=0\nwithout changing their status letter \u2014 abundance_table_species.tsv correctly excludes them",
         comma(n_qc_failed)),
       x = NULL, y = "Read count") +
  theme_amplicon +
  theme(legend.position = "none", plot.subtitle = element_text(size = 8.5))

ggsave("figs/01c_qc_filter_diagnostic.png", p1c, width = 7, height = 5.5, dpi = 150)
cat("[Figure 1c] QC filter diagnostic saved.\n")

# =============================================================================
# FIGURE 2 — Phylum-level composition (bar)
# =============================================================================
phylum_comp_full <- ab %>%
  mutate(phylum_clean = ifelse(is_unclassified, "Unclassified", phylum)) %>%
  group_by(phylum_clean) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(rel = count / sum(count)) %>%
  arrange(desc(count))

write.csv(phylum_comp_full, "tables/phylum_composition.csv", row.names = FALSE, quote = FALSE)

# Keep top phyla and group lower ones into "Other phyla"
N_KEEP <- 7
phylum_comp <- phylum_comp_full %>%
  mutate(rank = row_number(),
         phylum_clean = ifelse(rank > N_KEEP, "Other phyla", phylum_clean)) %>%
  group_by(phylum_clean) %>%
  summarise(count = sum(count), rel = sum(rel), .groups = "drop") %>%
  arrange(desc(count))

# Reorder factor level by abundance using base R reorder instead of forcats
phylum_comp$phylum_clean <- reorder(factor(phylum_comp$phylum_clean), phylum_comp$count)

n_other <- nrow(phylum_comp_full) - N_KEEP

p2 <- ggplot(phylum_comp, aes(x = phylum_clean, y = rel, fill = phylum_clean)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", 100*rel)), hjust = -0.15, size = 3.3) +
  scale_fill_manual(values = c(get_phylum_colors(levels(phylum_comp$phylum_clean)),
                                "Other phyla" = "#41566b")) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(title = "Phylum-Level Community Composition",
       subtitle = sprintf("Bioslurry sample \u00b7 %s total reads \u00b7 'Other phyla' = %d phyla each <0.5%%",
                           comma(total_reads), n_other),
       x = NULL, y = "Relative abundance") +
  theme_amplicon +
  theme(legend.position = "none")

ggsave("figs/02_phylum_composition.png", p2, width = 7.5, height = 5, dpi = 150)
cat("[Figure 2] Phylum composition saved.\n")

# =============================================================================
# FIGURE 2c — Phylum-level composition (stacked bar)
# =============================================================================
phylum_stack <- phylum_comp
# Reverse factor levels so the most abundant phylum sits at the bottom of the stack
phylum_stack$phylum_clean <- factor(phylum_stack$phylum_clean,
                                    levels = rev(levels(phylum_stack$phylum_clean)))
phylum_stack$sample <- sample_col

p2c <- ggplot(phylum_stack, aes(x = sample, y = rel, fill = phylum_clean)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = ifelse(rel > 0.02, sprintf("%.1f%%", 100*rel), "")),
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = get_phylum_colors(levels(phylum_stack$phylum_clean)),
                    name = "Phylum") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Phylum-Level Community Composition (Stacked)",
       subtitle = sprintf("Bioslurry sample \u00b7 %s total reads \u00b7 includes unclassified fraction",
                           comma(total_reads)),
       x = NULL, y = "Relative abundance") +
  theme_amplicon +
  theme(axis.text.x = element_text(face = "bold", size = 10))

ggsave("figs/02c_phylum_stacked.png", p2c, width = 7, height = 7, dpi = 150)
cat("[Figure 2c] Phylum stacked barplot saved.\n")

# =============================================================================
# FIGURE 2b — Family-level composition (bar)
# =============================================================================
family_comp <- ab %>%
  filter(!is_unclassified) %>%
  group_by(phylum, family) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(desc(count))

top_family <- family_comp %>% slice_head(n = 15)
top_family$family <- reorder(factor(top_family$family), top_family$count)

p2b <- ggplot(top_family, aes(x = family, y = count, fill = phylum)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(count)), hjust = -0.15, size = 3) +
  scale_fill_manual(values = get_phylum_colors(unique(top_family$phylum)), name = "Phylum") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(title = "Top 15 Families by Read Count",
       subtitle = sprintf("Out of %d families detected \u00b7 colored by phylum",
                           n_distinct(family_comp$family)),
       x = NULL, y = "Read count") +
  theme_amplicon +
  theme(legend.position = "right")

ggsave("figs/02b_family_composition.png", p2b, width = 9, height = 6, dpi = 150)
cat("[Figure 2b] Family composition barplot saved.\n")

write.csv(family_comp, "tables/family_composition.csv", row.names = FALSE, quote = FALSE)

# =============================================================================
# FIGURE 2d — Family-level composition (stacked bar)
# =============================================================================
N_FAMILY_KEEP <- 10

family_stacked <- ab %>%
  filter(!is_unclassified) %>%
  group_by(family) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(
    rank = row_number(),
    family_label = ifelse(rank <= N_FAMILY_KEEP, family, "Other families")
  ) %>%
  group_by(family_label) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(rel = count / sum(count)) %>%
  arrange(desc(count)) %>%
  mutate(family_label = factor(family_label, levels = rev(family_label)))

family_stacked$sample <- sample_col

# Curated 11-color palette (top 10 taxa + grey for "Other")
family_stack_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#4e79a7", "#e15759", "#76b7b2",
  "#bdbdbd"
)
names(family_stack_colors) <- levels(family_stacked$family_label)

p2d <- ggplot(family_stacked, aes(x = sample, y = rel, fill = family_label)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = ifelse(rel > 0.02, sprintf("%.1f%%", 100*rel), "")),
            position = position_stack(vjust = 0.5), size = 2.8, fontface = "bold") +
  scale_fill_manual(values = family_stack_colors, name = "Family") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Family-Level Community Composition (Stacked)",
       subtitle = sprintf("Top %d families \u00b7 relative to %s classified reads",
                           N_FAMILY_KEEP, comma(n_classified)),
       x = NULL, y = "Relative abundance") +
  theme_amplicon +
  theme(axis.text.x = element_text(face = "bold", size = 10))

ggsave("figs/02d_family_stacked.png", p2d, width = 7.5, height = 7, dpi = 150)
cat("[Figure 2d] Family stacked barplot saved.\n")

# =============================================================================
# FIGURE 3 — Top 20 genera (stacked by phylum)
# =============================================================================
genus_comp <- ab %>%
  filter(!is_unclassified) %>%
  group_by(phylum, genus) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(desc(count))

top_genus <- genus_comp %>% slice_head(n = 20)
# Reorder by count using base R reorder
top_genus$genus <- reorder(factor(top_genus$genus), top_genus$count)

p3 <- ggplot(top_genus, aes(x = genus, y = count, fill = phylum)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(count)), hjust = -0.15, size = 3) +
  scale_fill_manual(values = get_phylum_colors(unique(top_genus$phylum)), name = "Phylum") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(title = "Top 20 Genera by Read Count",
       subtitle = sprintf("Out of %d genera detected \u00b7 colored by phylum",
                           n_distinct(genus_comp$genus)),
       x = NULL, y = "Read count") +
  theme_amplicon +
  theme(legend.position = "right")

ggsave("figs/03_top20_genera.png", p3, width = 9, height = 7, dpi = 150)
cat("[Figure 3] Top 20 genera saved.\n")

write.csv(genus_comp, "tables/genus_composition.csv", row.names = FALSE, quote = FALSE)

# =============================================================================
# FIGURE 3b — Genus-level composition (stacked bar)
# =============================================================================
N_GENUS_KEEP <- 10

genus_stacked <- ab %>%
  filter(!is_unclassified) %>%
  group_by(genus) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(
    rank = row_number(),
    genus_label = ifelse(rank <= N_GENUS_KEEP, genus, "Other genera")
  ) %>%
  group_by(genus_label) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(rel = count / sum(count)) %>%
  arrange(desc(count)) %>%
  mutate(genus_label = factor(genus_label, levels = rev(genus_label)))

genus_stacked$sample <- sample_col

genus_stack_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#4e79a7", "#e15759", "#76b7b2",
  "#bdbdbd"
)
names(genus_stack_colors) <- levels(genus_stacked$genus_label)

p3b <- ggplot(genus_stacked, aes(x = sample, y = rel, fill = genus_label)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = ifelse(rel > 0.02, sprintf("%.1f%%", 100*rel), "")),
            position = position_stack(vjust = 0.5), size = 2.8, fontface = "bold") +
  scale_fill_manual(values = genus_stack_colors, name = "Genus") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Genus-Level Community Composition (Stacked)",
       subtitle = sprintf("Top %d genera \u00b7 relative to %s classified reads",
                           N_GENUS_KEEP, comma(n_classified)),
       x = NULL, y = "Relative abundance") +
  theme_amplicon +
  theme(axis.text.x = element_text(face = "bold", size = 10))

ggsave("figs/03b_genus_stacked.png", p3b, width = 7.5, height = 7, dpi = 150)
cat("[Figure 3b] Genus stacked barplot saved.\n")

# =============================================================================
# FIGURE 4 — Top 15 species
# =============================================================================
top_species <- ab %>%
  filter(!is_unclassified) %>%
  arrange(desc(count)) %>%
  slice_head(n = 15)
# Reorder by count using base R reorder
top_species$species_lab <- reorder(factor(top_species$species), top_species$count)

p4 <- ggplot(top_species, aes(x = species_lab, y = rel_abund, fill = phylum)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%  (n=%s)", 100*rel_abund, comma(count))),
            hjust = -0.05, size = 3) +
  scale_fill_manual(values = get_phylum_colors(unique(top_species$phylum)), name = "Phylum") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.35))) +
  coord_flip() +
  labs(title = "Top 15 Species by Relative Abundance",
       subtitle = sprintf("Out of %d species-level taxa \u00b7 relative to %s classified reads",
                           sum(!ab$is_unclassified), comma(n_classified)),
       x = NULL, y = "Relative abundance (of classified reads)") +
  theme_amplicon +
  theme(legend.position = "right", axis.text.y = element_text(face = "italic"))

ggsave("figs/04_top15_species.png", p4, width = 9.5, height = 6.5, dpi = 150)
cat("[Figure 4] Top 15 species saved.\n")

# =============================================================================
# 5. ALPHA DIVERSITY (single-sample indices via vegan)
# =============================================================================
classified_counts <- ab %>% filter(!is_unclassified) %>% pull(count)

# Convert counts explicitly to integers to ensure perfect compatibility with vegan
classified_counts_int <- as.integer(round(classified_counts))

S <- length(classified_counts_int)
shannon <- diversity(classified_counts_int, index = "shannon")
simpson <- diversity(classified_counts_int, index = "simpson")
invsimpson <- diversity(classified_counts_int, index = "invsimpson")
pielou <- shannon / log(S)
ens <- exp(shannon)
bp_dominance <- max(classified_counts_int) / sum(classified_counts_int)
chao1 <- tryCatch(estimateR(classified_counts_int)["S.chao1"], error = function(e) NA)
fisher_alpha <- tryCatch(fisher.alpha(classified_counts_int), error = function(e) NA)

alpha_tbl <- tibble(
  Index = c("Observed species richness (S)", "Chao1 (estimated richness)",
            "Shannon (H)", "Effective number of species (e^H)",
            "Simpson's D (1-sum p^2)", "Inverse Simpson", "Pielou's evenness (J)",
            "Fisher's alpha", "Berger-Parker dominance"),
  Value = c(S, chao1, shannon, ens, simpson, invsimpson, pielou, fisher_alpha, bp_dominance)
) %>% mutate(Value = round(Value, 3))

write.csv(alpha_tbl, "tables/alpha_diversity.csv", row.names = FALSE, quote = FALSE)
cat("\n[5] Alpha diversity indices (species-level, classified reads only):\n")
print(alpha_tbl, n = Inf)

# =============================================================================
# FIGURE 5 — Within-sample rarefaction curve (richness vs sequencing depth)
# =============================================================================
# Replaced slow, noisy bootstrap/resampling loop with exact analytical rarefaction.
# This scales perfectly, runs in milliseconds, and avoids random noise.
depth_points <- unique(round(seq(100, sum(classified_counts_int), length.out = 25)))

# vegan::rarefy returns expected species count (row 1) and standard error (row 2)
rare_analytical <- rarefy(classified_counts_int, sample = depth_points, se = TRUE)

rare_results <- tibble(
  depth = depth_points,
  mean_richness = rare_analytical[1, ],
  sd_richness = rare_analytical[2, ]
)

p5 <- ggplot(rare_results, aes(x = depth, y = mean_richness)) +
  geom_ribbon(aes(ymin = mean_richness - sd_richness, ymax = mean_richness + sd_richness),
              fill = "#1b9e77", alpha = 0.2) +
  geom_line(color = "#1b9e77", linewidth = 1) +
  geom_point(color = "#1b9e77", size = 1.5) +
  geom_vline(xintercept = sum(classified_counts_int), linetype = "dashed", color = "grey50") +
  annotate("text", x = sum(classified_counts_int), y = 50,
           label = "actual\ndepth ", hjust = 1, vjust = 0, color = "grey40", size = 3) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0.01, 0.08))) +
  labs(title = "Within-Sample Rarefaction Curve",
       subtitle = "Species richness vs. sequencing depth (exact analytical calculation \u00b1 SE)",
       x = "Reads subsampled", y = "Observed species richness") +
  theme_amplicon

ggsave("figs/05_rarefaction_curve.png", p5, width = 8, height = 5.5, dpi = 150)
cat("\n[Figure 5] Rarefaction curve saved.\n")

# =============================================================================
# FIGURE 5b — Diversity distributions (Bootstrap Resampling Boxplots)
# =============================================================================
set.seed(42)  # For reproducibility of rrarefy
# Dynamic depth calculation to be robust on lower read counts
subsample_depth <- min(50000, round(sum(classified_counts_int) * 0.9, -2))
n_bootstrap <- 100

cat(sprintf("\n[5b] Performing %d bootstrap resamples at depth = %s reads...\n", n_bootstrap, comma(subsample_depth)))

boot_results <- lapply(seq_len(n_bootstrap), function(i) {
  sub_sample <- rrarefy(matrix(classified_counts_int, nrow = 1), subsample_depth)[1, ]
  sub_counts <- sub_sample[sub_sample > 0]
  
  S_boot <- length(sub_counts)
  shannon_boot <- diversity(sub_counts, index = "shannon")
  simpson_boot <- diversity(sub_counts, index = "simpson")
  chao1_boot <- as.numeric(tryCatch(estimateR(sub_counts)["S.chao1"], error = function(e) NA))
  
  c(Observed = S_boot, Shannon = shannon_boot, Simpson = simpson_boot, Chao1 = chao1_boot)
})

boot_df <- as.data.frame(do.call(rbind, boot_results))
write.csv(boot_df, "tables/bootstrap_diversity.csv", row.names = FALSE, quote = FALSE)

# Pivot to long format for ggplot faceting
boot_long <- boot_df %>%
  pivot_longer(cols = everything(), names_to = "Index", values_to = "Value")

p5b <- ggplot(boot_long, aes(x = "", y = Value, fill = Index)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6, outlier.alpha = 0.5, width = 0.5) +
  facet_wrap(~ Index, scales = "free", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Diversity Index Distributions (100 Bootstrap Resamples)",
       subtitle = sprintf("Subsampled at depth of %s reads \u00b7 variance reflects rarefaction noise, not biological replication",
                           comma(subsample_depth)),
       x = NULL, y = "Value") +
  theme_amplicon +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave("figs/05b_diversity_boxplots.png", p5b, width = 8, height = 6.5, dpi = 150)
cat("[Figure 5b] Diversity boxplots saved.\n")

# =============================================================================
# FIGURE 6 — Read length by classification status (violin/box) - QC diagnostic
# =============================================================================
p6 <- ggplot(reads, aes(x = status_label_corrected, y = read_length, fill = status_label_corrected)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.4, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
  coord_cartesian(ylim = c(1200, 1800)) +
  scale_y_continuous(labels = comma) +
  labs(title = "Read Length vs. Classification Outcome",
       subtitle = "Unclassified reads show wider length variance despite similar median",
       x = NULL, y = "Read length (bp)") +
  theme_amplicon +
  theme(legend.position = "none")

ggsave("figs/06_length_vs_classification.png", p6, width = 6.5, height = 5.5, dpi = 150)
cat("[Figure 6] Read length vs classification status saved.\n")

# Quick summary stat for that comparison
len_summary <- reads %>%
  group_by(status_label_corrected) %>%
  summarise(
    n = n(),
    median_length = median(read_length),
    mean_length = round(mean(read_length), 1),
    sd_length = round(sd(read_length), 1),
    .groups = "drop"
  )
write.csv(len_summary, "tables/read_length_by_status.csv", row.names = FALSE, quote = FALSE)
cat("\nRead length summary by status:\n")
print(len_summary)

# =============================================================================
# 7. TAXONOMIC SUMMARY TABLE (richness per rank)
# =============================================================================
richness_per_rank <- tibble(
  Rank = c("Phylum","Class","Order","Family","Genus","Species"),
  N_detected = c(
    n_distinct(ab$phylum[!ab$is_unclassified]),
    n_distinct(ab$class[!ab$is_unclassified]),
    n_distinct(ab$order[!ab$is_unclassified]),
    n_distinct(ab$family[!ab$is_unclassified]),
    n_distinct(ab$genus[!ab$is_unclassified]),
    n_distinct(ab$species[!ab$is_unclassified])
  )
)
write.csv(richness_per_rank, "tables/richness_per_rank.csv", row.names = FALSE, quote = FALSE)
cat("\n[7] Taxa detected per rank:\n")
print(richness_per_rank)

cat("\n=== DONE: all figures in figs/, all tables in tables/ ===\n")
cat("=== For Sankey diagrams, run convert_to_kreport.R and load the .kreport in Pavian ===\n")
