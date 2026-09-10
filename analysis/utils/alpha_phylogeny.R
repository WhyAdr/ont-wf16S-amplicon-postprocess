# =============================================================================
# Optional phylogenetic alpha-diversity contract
# =============================================================================

read_phylogenetic_tip_map <- function(path) {
  if (is.null(path)) return(NULL)
  header <- readLines(path, n = 1L, warn = FALSE)
  if (!length(header) || !identical(strsplit(header, "\t", fixed = TRUE)[[1]],
                                    c("TaxonPath", "TipLabel"))) {
    stop("input.phylogenetic_tip_map must have exactly: TaxonPath, TipLabel.",
         call. = FALSE)
  }
  map <- read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE,
                    colClasses = "character")
  if (!nrow(map) || anyNA(map) || any(!nzchar(map$TaxonPath)) ||
      any(!nzchar(map$TipLabel)) || anyDuplicated(map$TaxonPath) ||
      anyDuplicated(map$TipLabel)) {
    stop("Phylogenetic tip-map paths and labels must be non-empty and unique.",
         call. = FALSE)
  }
  if (any(grepl("^Unclassified(?:;|$)", map$TaxonPath, perl = TRUE))) {
    stop("Unclassified must not be mapped to a phylogenetic tip.", call. = FALSE)
  }
  map
}

prepare_alpha_phylogeny <- function(cfg, positive_taxon_paths) {
  tree_path <- cfg$input$phylogenetic_tree %||% NULL
  map_path <- cfg$input$phylogenetic_tip_map %||% NULL
  affected <- c("faith_pd", "psr", "pse")
  if (is.null(tree_path)) {
    if (!is.null(map_path)) {
      stop("input.phylogenetic_tip_map requires input.phylogenetic_tree.", call. = FALSE)
    }
    return(list(
      enabled = FALSE,
      skip = data.frame(
        Status = "Skipped", Reason = "input.phylogenetic_tree was not supplied",
        TreePath = "", TreeSHA256 = "", TipMapPath = "", TipMapSHA256 = "",
        AffectedMetrics = paste(affected, collapse = ","), stringsAsFactors = FALSE
      )
    ))
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("A configured phylogenetic tree requires the locked 'ape' package.", call. = FALSE)
  }
  tree <- tryCatch(ape::read.tree(tree_path), error = function(e) {
    stop(sprintf("Could not parse input.phylogenetic_tree: %s", conditionMessage(e)),
         call. = FALSE)
  })
  if (is.null(tree) || !inherits(tree, "phylo") || length(tree$tip.label) < 1L) {
    stop("input.phylogenetic_tree must contain one valid Newick tree.", call. = FALSE)
  }
  if (anyNA(tree$tip.label) || any(!nzchar(tree$tip.label)) || anyDuplicated(tree$tip.label)) {
    stop("Phylogenetic tree tip labels must be non-empty and unique.", call. = FALSE)
  }
  if (!ape::is.rooted(tree)) {
    stop("Phylogenetic tree must be rooted.", call. = FALSE)
  }
  if (is.null(tree$edge.length) || length(tree$edge.length) != nrow(tree$edge) ||
      anyNA(tree$edge.length) || any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0)) {
    stop("Phylogenetic tree must have finite, non-negative branch lengths.", call. = FALSE)
  }

  map <- read_phylogenetic_tip_map(map_path)
  if (is.null(map)) {
    map <- data.frame(TaxonPath = positive_taxon_paths,
                      TipLabel = positive_taxon_paths, stringsAsFactors = FALSE)
  }
  missing_paths <- setdiff(positive_taxon_paths, map$TaxonPath)
  if (length(missing_paths)) {
    stop(sprintf("Phylogenetic tip map is missing positive TaxonPath '%s'.",
                 missing_paths[[1]]), call. = FALSE)
  }
  mapped <- map$TipLabel[match(positive_taxon_paths, map$TaxonPath)]
  missing_tips <- setdiff(mapped, tree$tip.label)
  if (length(missing_tips)) {
    stop(sprintf("Mapped phylogenetic tip '%s' is absent from the tree.",
                 missing_tips[[1]]), call. = FALSE)
  }
  keep <- unique(mapped)
  extra <- setdiff(tree$tip.label, keep)
  mapping <- stats::setNames(map$TipLabel, map$TaxonPath)
  list(
    # Keep the original root and its stem paths. Physically dropping an
    # outgroup can re-root/collapse the retained topology and change rooted PD
    # and VCV-derived metrics. Extra tips are excluded by subsetting observed
    # names during each metric calculation instead.
    enabled = TRUE, tree = tree, mapping = mapping,
    tree_path = tree_path, tree_sha256 = compute_file_hash(tree_path),
    tip_map_path = map_path,
    tip_map_sha256 = if (is.null(map_path)) NA_character_ else compute_file_hash(map_path),
    original_tip_count = length(tree$tip.label), retained_tip_count = length(keep),
    pruned_tip_count = length(extra),
    pruning_method = "extra tips excluded from metric traversal; original root retained"
  )
}

faith_pd_rooted <- function(tree, observed_tips) {
  tip_indices <- match(observed_tips, tree$tip.label)
  active <- unique(tip_indices[!is.na(tip_indices)])
  if (!length(active)) return(NA_real_)
  included <- rep(FALSE, nrow(tree$edge))
  repeat {
    hit <- which(tree$edge[, 2] %in% active & !included)
    if (!length(hit)) break
    included[hit] <- TRUE
    active <- unique(c(active, tree$edge[hit, 1]))
  }
  sum(tree$edge.length[included])
}

calc_phylogenetic_alpha <- function(counts, phylogeny,
                                     registry = alpha_metric_registry()) {
  counts <- validate_alpha_counts(counts)
  if (is.null(names(counts)) || any(!nzchar(names(counts)))) {
    stop("Phylogenetic alpha counts must be named by TaxonPath.", call. = FALSE)
  }
  if (!isTRUE(phylogeny$enabled)) {
    return(do.call(rbind, lapply(c("faith_pd", "psr", "pse"), function(key) {
      alpha_metric_result(registry, key, NA_real_, FALSE, "phylogenetic module skipped")
    })))
  }
  tip_labels <- unname(phylogeny$mapping[names(counts)])
  if (anyNA(tip_labels)) stop("Positive alpha count lacks a phylogenetic mapping.", call. = FALSE)
  names(counts) <- tip_labels
  S <- length(counts)
  N <- sum(counts)
  pd <- faith_pd_rooted(phylogeny$tree, names(counts))
  psr <- pse <- NA_real_
  if (S > 1L) {
    covariance <- ape::vcv.phylo(phylogeny$tree, corr = TRUE)
    covariance <- covariance[names(counts), names(counts), drop = FALSE]
    psv <- (S * sum(diag(covariance)) - sum(covariance)) / (S * (S - 1))
    psr <- psv * S
    abundance <- as.numeric(counts)
    mean_abundance <- mean(abundance)
    numerator <- N * sum(diag(covariance) * abundance) -
      as.numeric(t(abundance) %*% covariance %*% abundance)
    denominator <- N^2 - N * mean_abundance
    if (is.finite(denominator) && denominator > 0) pse <- numerator / denominator
  }
  rbind(
    alpha_metric_result(registry, "faith_pd", pd, is.finite(pd),
                        "no mapped observed taxon"),
    alpha_metric_result(registry, "psr", psr, S > 1L && is.finite(psr),
                        if (S <= 1L) "richness <= 1" else "PSR was undefined"),
    alpha_metric_result(registry, "pse", pse, S > 1L && is.finite(pse),
                        if (S <= 1L) "richness <= 1" else "PSE was undefined")
  )
}

alpha_phylogeny_provenance <- function(phylogeny) {
  if (!isTRUE(phylogeny$enabled)) return(phylogeny$skip)
  data.frame(
    Status = "Completed", Reason = "",
    TreePath = phylogeny$tree_path, TreeSHA256 = phylogeny$tree_sha256,
    TipMapPath = phylogeny$tip_map_path %||% "",
    TipMapSHA256 = phylogeny$tip_map_sha256 %||% "",
    OriginalTipCount = phylogeny$original_tip_count,
    RetainedTipCount = phylogeny$retained_tip_count,
    PrunedTipCount = phylogeny$pruned_tip_count,
    PruningMethod = phylogeny$pruning_method,
    AffectedMetrics = "faith_pd,psr,pse", stringsAsFactors = FALSE
  )
}
