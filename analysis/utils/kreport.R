# =============================================================================
# Kraken Report (.kreport) Formatter and Tree Validator
# =============================================================================

suppressMessages(library(jsonlite))

RANKS_8 <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
RANK_CODES_8 <- c("D", "K", "P", "C", "O", "F", "G", "S")

build_kreport_tree <- function(lineages_str, counts) {
  lineages <- strsplit(lineages_str, ";")
  counts <- round(as.numeric(counts))

  node_env <- new.env(hash = TRUE, parent = emptyenv())

  for (i in seq_along(lineages)) {
    lin <- lineages[[i]]
    cnt <- counts[i]

    if (lin[1] == "Unclassified" || cnt == 0) next

    for (j in seq_along(lin)) {
      path <- paste(lin[1:j], collapse = ";")

      if (exists(path, envir = node_env)) {
        node <- get(path, envir = node_env)
        node$reads_clade <- node$reads_clade + cnt
        if (j == length(lin)) node$reads_taxon <- node$reads_taxon + cnt
        assign(path, node, envir = node_env)
      } else {
        assign(path, list(
          path        = path,
          name        = lin[j],
          depth       = j,
          rank_code   = RANK_CODES_8[j],
          reads_clade = cnt,
          reads_taxon = if (j == length(lin)) cnt else 0L,
          parent_path = if (j > 1) paste(lin[1:(j - 1)], collapse = ";") else ""
        ), envir = node_env)
      }
    }
  }

  all_paths <- ls(node_env)
  if (length(all_paths) == 0) {
    return(data.frame(
      path = character(0), name = character(0), depth = integer(0),
      rank_code = character(0), reads_clade = integer(0),
      reads_taxon = integer(0), parent_path = character(0),
      stringsAsFactors = FALSE
    ))
  }

  nodes_list <- lapply(all_paths, function(p) get(p, envir = node_env))
  nodes_df <- do.call(rbind.data.frame, c(nodes_list, stringsAsFactors = FALSE))

  # Depth-first search sorting with abundance tie-breaking
  dfs_order <- function(parent) {
    children <- nodes_df[nodes_df$parent_path == parent, , drop = FALSE]
    if (nrow(children) == 0) return(character(0))
    children <- children[order(-children$reads_clade, children$name), ]
    res <- character(0)
    for (idx in seq_len(nrow(children))) {
      c_path <- children$path[idx]
      res <- c(res, c_path, dfs_order(c_path))
    }
    res
  }

  ordered_paths <- dfs_order("")
  nodes_sorted <- nodes_df[match(ordered_paths, nodes_df$path), ]
  rownames(nodes_sorted) <- NULL
  nodes_sorted
}

validate_kreport_tree <- function(nodes_df, total_reads, uncl_reads) {
  cl_reads <- total_reads - uncl_reads

  # 1. Total reads check
  root_nodes <- nodes_df[nodes_df$parent_path == "", , drop = FALSE]
  sum_root_clade <- sum(root_nodes$reads_clade)

  if (sum_root_clade != cl_reads) {
    stop(sprintf("Kreport tree validation error: root clades sum (%d) != total classified reads (%d)",
                 sum_root_clade, cl_reads), call. = FALSE)
  }

  if (uncl_reads + sum_root_clade != total_reads) {
    stop(sprintf("Kreport tree validation error: unclassified (%d) + root (%d) != total reads (%d)",
                 uncl_reads, sum_root_clade, total_reads), call. = FALSE)
  }

  # 2. Clade = direct + sum(child clades) check
  for (i in seq_len(nrow(nodes_df))) {
    p <- nodes_df$path[i]
    clade_cnt <- nodes_df$reads_clade[i]
    direct_cnt <- nodes_df$reads_taxon[i]

    children <- nodes_df[nodes_df$parent_path == p, , drop = FALSE]
    child_sum <- if (nrow(children) > 0) sum(children$reads_clade) else 0L

    if (clade_cnt != (direct_cnt + child_sum)) {
      stop(sprintf("Kreport tree validation error at '%s': clade (%d) != direct (%d) + child sum (%d)",
                   p, clade_cnt, direct_cnt, child_sum), call. = FALSE)
    }
  }

  invisible(TRUE)
}

format_kreport_lines <- function(nodes_sorted, total_reads, uncl_reads, taxid_cache = list()) {
  cl_reads <- total_reads - uncl_reads

  lines <- character(nrow(nodes_sorted) + 2)

  # Line 1: unclassified
  lines[1] <- sprintf("%.2f\t%.0f\t%.0f\tU\t0\tunclassified",
                      100 * uncl_reads / total_reads, uncl_reads, uncl_reads)

  # Line 2: root
  lines[2] <- sprintf("%.2f\t%.0f\t%.0f\tR\t1\troot",
                      100 * cl_reads / total_reads, cl_reads, 0L)

  for (i in seq_len(nrow(nodes_sorted))) {
    indent <- strrep("  ", nodes_sorted$depth[i])
    p <- nodes_sorted$path[i]
    taxid <- taxid_cache[[p]]
    if (is.null(taxid)) taxid <- 0L

    pct <- 100 * nodes_sorted$reads_clade[i] / total_reads
    lines[i + 2] <- sprintf("%.2f\t%.0f\t%.0f\t%s\t%s\t%s%s",
                            pct,
                            nodes_sorted$reads_clade[i],
                            nodes_sorted$reads_taxon[i],
                            nodes_sorted$rank_code[i],
                            as.character(taxid),
                            indent,
                            nodes_sorted$name[i])
  }

  lines
}
