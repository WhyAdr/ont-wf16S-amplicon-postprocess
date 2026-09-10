#!/usr/bin/env Rscript
# =============================================================================
# Convert wf-16S abundance_table_species.tsv → Kraken-style report (.kreport)
# for interactive Sankey visualization in Pavian
#
# Usage:
#   Rscript convert_to_kreport.R
#   Then upload tables/AmbarAyunda_minimap2_16S.kreport to Pavian:
#   https://fbreitwieser.shinyapps.io/pavian/
#
# Input:  output_AAy/abundance_table_species.tsv  (semicolon-delimited 8-rank lineages)
# Output: tables/<sample_name>.kreport             (6-column Kraken report format)
# =============================================================================

suppressMessages(library(dplyr))
suppressMessages(library(jsonlite))

# --- Configuration ---
ABUNDANCE_FILE <- "output_AAy/abundance_table_species.tsv"
RANKS          <- c("superkingdom","clade","phylum","class","order","family","genus","species")
RANK_CODES     <- c("D",          "D1",   "P",     "C",    "O",    "F",     "G",    "S")

# --- Read abundance table ---
ab_raw <- read.delim(ABUNDANCE_FILE, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE)
sample_col  <- colnames(ab_raw)[2]
total_reads <- as.integer(sum(ab_raw[[sample_col]]))

cat(sprintf("[kreport] Sample: %s | Total reads: %s\n",
            sample_col, format(total_reads, big.mark = ",")))

# --- Parse lineages and build node tree ---
# Each node stores: name, depth, rank_code, reads_clade (this node + all
# descendants), reads_taxon (reads assigned directly to this node), and
# parent_path (for DFS traversal).
# An environment is used as a hash map for O(1) path lookups.

lineages <- strsplit(ab_raw$tax, ";")
counts   <- as.integer(ab_raw[[sample_col]])

node_env <- new.env(hash = TRUE, parent = emptyenv())

for (i in seq_along(lineages)) {
  lin <- lineages[[i]]
  cnt <- counts[i]

  # Skip the synthetic "Unclassified" lineage — handled separately below
  if (lin[1] == "Unclassified") next

  for (j in seq_along(lin)) {
    path <- paste(lin[1:j], collapse = ";")

    if (exists(path, envir = node_env)) {
      node <- get(path, envir = node_env)
      node$reads_clade <- node$reads_clade + cnt
      if (j == length(lin)) node$reads_taxon <- node$reads_taxon + cnt
      assign(path, node, envir = node_env)
    } else {
      assign(path, list(
        name        = lin[j],
        depth       = j,
        rank_code   = RANK_CODES[j],
        reads_clade = cnt,
        reads_taxon = if (j == length(lin)) cnt else 0L,
        parent_path = if (j > 1) paste(lin[1:(j - 1)], collapse = ";") else ""
      ), envir = node_env)
    }
  }
}

# Extract all nodes into a data frame
all_paths <- ls(node_env)
nodes <- do.call(rbind, lapply(all_paths, function(p) {
  n <- get(p, envir = node_env)
  data.frame(path = p, name = n$name, depth = n$depth,
             rank_code = n$rank_code, reads_clade = n$reads_clade,
             reads_taxon = n$reads_taxon, parent_path = n$parent_path,
             stringsAsFactors = FALSE)
}))

cat(sprintf("[kreport] Built tree: %d unique nodes from %d classified lineages\n",
            nrow(nodes), sum(sapply(lineages, function(x) x[1] != "Unclassified"))))

# --- Sort in depth-first order (most-abundant children first at each level) ---
# This ordering controls the visual layout of the Pavian Sankey — the most
# abundant taxa appear at the top of each rank column.

dfs_order <- function(parent, nodes_df) {
  children <- nodes_df[nodes_df$parent_path == parent, , drop = FALSE]
  if (nrow(children) == 0) return(character())

  children <- children[order(-children$reads_clade), ]
  result <- character()
  for (i in seq_len(nrow(children))) {
    result <- c(result, children$path[i])
    result <- c(result, dfs_order(children$path[i], nodes_df))
  }
  result
}

ordered_paths <- dfs_order("", nodes)
nodes_sorted  <- nodes[match(ordered_paths, nodes$path), ]

# --- Compute unclassified reads ---
uncl_idx   <- sapply(lineages, function(x) x[1] == "Unclassified")
uncl_reads <- sum(counts[uncl_idx])
cl_reads   <- total_reads - uncl_reads

# --- Write kreport ---
# Kraken report format (tab-separated, 6 columns):
#   1. Percentage of reads rooted at this node
#   2. Number of reads rooted at this node (clade)
#   3. Number of reads assigned directly to this node (taxon)
#   4. Rank code (U/R/D/D1/P/C/O/F/G/S)
#   5. NCBI taxonomy ID (0 = placeholder; Pavian does not require real taxids)
#   6. Scientific name (indented with 2 spaces per depth level)

dir.create("tables", showWarnings = FALSE)
output_file <- sprintf("tables/%s.kreport", sample_col)
con <- file(output_file, "w")

# Line 1: unclassified
writeLines(sprintf("%.2f\t%d\t%d\tU\t0\tunclassified",
                   100 * uncl_reads / total_reads, uncl_reads, uncl_reads), con)

# Line 2: root
writeLines(sprintf("%.2f\t%d\t%d\tR\t1\troot",
                   100 * cl_reads / total_reads, cl_reads, 0L), con)

# Load taxonomy cache (maps paths to actual NCBI TaxIDs)
cache_file <- "output_AAy/taxonomy_cache.json"
if (!file.exists(cache_file)) {
  cat("[kreport] Taxonomy cache not found. Resolving taxonomy via NCBI...\n")
  system("python fetch_ncbi_taxonomy.py")
}
cache <- jsonlite::fromJSON(cache_file)

# Classified tree nodes (depth-first, abundance-sorted)
for (i in seq_len(nrow(nodes_sorted))) {
  indent <- strrep("  ", nodes_sorted$depth[i])
  path <- nodes_sorted$path[i]
  taxid <- cache[[path]]
  if (is.null(taxid)) taxid <- 0
  
  writeLines(sprintf("%.2f\t%d\t%d\t%s\t%s\t%s%s",
                     100 * nodes_sorted$reads_clade[i] / total_reads,
                     nodes_sorted$reads_clade[i],
                     nodes_sorted$reads_taxon[i],
                     nodes_sorted$rank_code[i],
                     as.character(taxid),
                     indent,
                     nodes_sorted$name[i]), con)
}

close(con)

cat(sprintf("\n[kreport] Written: %s\n", output_file))
cat(sprintf("          %d nodes + 2 header lines = %d total lines\n",
            nrow(nodes_sorted), nrow(nodes_sorted) + 2))
cat(sprintf("          Classified: %s (%.1f%%) | Unclassified: %s (%.1f%%)\n",
            format(cl_reads, big.mark = ","), 100 * cl_reads / total_reads,
            format(uncl_reads, big.mark = ","), 100 * uncl_reads / total_reads))
cat("\n  Next step: upload this .kreport file to Pavian for interactive Sankey:\n")
cat("  https://fbreitwieser.shinyapps.io/pavian/\n")
