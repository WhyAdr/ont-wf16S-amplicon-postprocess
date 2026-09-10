# =============================================================================
# Mathematical & Alpha Diversity Metrics Utility
# =============================================================================

suppressMessages({
  library(vegan)
  library(dplyr)
})

derive_sample_seed <- function(seed, sample_id) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max ||
      abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
    stop("Base seed must be an integer in the set.seed() range.", call. = FALSE)
  }
  if (!is.character(sample_id) || length(sample_id) != 1L || is.na(sample_id) ||
      !nzchar(sample_id)) {
    stop("Sample ID for seed derivation must be one non-empty string.", call. = FALSE)
  }
  hash_part <- strtoi(substr(
    digest::digest(sample_id, algo = "xxhash32", serialize = FALSE), 1L, 7L
  ), base = 16L)
  as.integer((as.double(seed) + as.double(hash_part)) %%
               as.double(.Machine$integer.max))
}

validate_alpha_counts <- function(counts) {
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) ||
      any(counts < 0) || any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop("Alpha-diversity counts must be finite, non-negative integers.", call. = FALSE)
  }
  original_names <- names(counts)
  # Preserve double storage after validating integer-valued counts. Coercing to
  # R's 32-bit integer type silently turns otherwise valid large read counts
  # into NA and corrupts full-depth metrics.
  counts <- round(counts)
  names(counts) <- original_names
  counts[counts > 0L]
}

validate_alpha_orders <- function(orders, name) {
  if (!is.numeric(orders) || !length(orders) || anyNA(orders) ||
      any(!is.finite(orders)) || any(orders < 0) || anyDuplicated(orders)) {
    stop(sprintf("'%s' must contain unique, finite, non-negative q orders.", name),
         call. = FALSE)
  }
  as.numeric(orders)
}

alpha_q_token <- function(q) {
  token <- format(q, scientific = FALSE, trim = TRUE, digits = 15)
  token <- sub("[.]0+$", "", token)
  token <- sub("([.][0-9]*[1-9])0+$", "\\1", token)
  gsub("[.]", "p", token)
}

alpha_q_label <- function(q) format(q, scientific = FALSE, trim = TRUE, digits = 15)

alpha_metric_registry <- function(hill_orders = c(0, 1, 2), renyi_orders = 1) {
  hill_orders <- validate_alpha_orders(hill_orders, "alpha.hill_orders")
  renyi_orders <- validate_alpha_orders(renyi_orders, "alpha.renyi_orders")

  fixed <- data.frame(
    MetricKey = c(
      "richness", "shannon", "ens", "simpson", "invsimpson", "fisher_alpha",
      "chao1", "ace", "pielou", "berger_parker", "heip", "evar", "mcintosh"
    ),
    MetricLabel = c(
      "Richness (S)", "Shannon (H)", "ENS (e^H)", "Simpson (D)",
      "Inv. Simpson", "Fisher's alpha", "Chao1", "ACE",
      "Pielou's evenness (J)", "Berger-Parker index", "Heip's evenness",
      "Smith-Wilson Evar", "McIntosh diversity index"
    ),
    Parameter = rep(NA_real_, 13L),
    FigureKey = c(rep("02b", 6L), rep("02c", 4L), rep("02d", 3L)),
    FacetOrder = c(seq_len(6L), seq_len(4L), seq_len(3L)),
    Formula = c(
      "number of positive-count taxa",
      "-sum(p_i * log(p_i))", "exp(shannon)", "1 - sum(p_i^2)",
      "1 / sum(p_i^2)", "vegan::fisher.alpha(counts)",
      "vegan::estimateR(counts)['S.chao1']",
      "vegan::estimateR(counts)['S.ACE']", "shannon / log(richness)",
      "max(n_i) / N", "(exp(shannon) - 1) / (richness - 1)",
      "1 - (2/pi) * atan(mean((log(n_i) - mean(log(n_i)))^2))",
      "(N - sqrt(sum(n_i^2))) / (N - sqrt(N))"
    ),
    Method = c(
      "closed_form", "vegan::diversity", "closed_form", "vegan::diversity",
      "vegan::diversity", "vegan::fisher.alpha", "vegan::estimateR",
      "vegan::estimateR", rep("closed_form", 5L)
    ),
    ValidDomain = c(
      rep("N > 0", 8L), "richness > 1", "N > 0", "richness > 1",
      "richness > 1", "N > 1"
    ),
    UndefinedRule = c(
      rep("NA when no classified reads", 8L),
      "NA when richness <= 1", "NA when no classified reads",
      "NA when richness <= 1", "NA when richness <= 1", "NA when N <= 1"
    ),
    stringsAsFactors = FALSE
  )

  renyi <- do.call(rbind, lapply(seq_along(renyi_orders), function(index) {
    q <- renyi_orders[[index]]
    data.frame(
      MetricKey = paste0("renyi_q", alpha_q_token(q)),
      MetricLabel = sprintf("R\u00e9nyi entropy (q=%s)", alpha_q_label(q)),
      Parameter = q,
      FigureKey = if (identical(q, 1)) "02d" else "02f",
      FacetOrder = if (identical(q, 1)) 4L else index,
      Formula = "log(sum(p_i^q))/(1-q), with Shannon limit at q=1",
      Method = "closed_form",
      ValidDomain = "N > 0 and q >= 0",
      UndefinedRule = "NA when no classified reads",
      stringsAsFactors = FALSE
    )
  }))
  hill <- do.call(rbind, lapply(seq_along(hill_orders), function(index) {
    q <- hill_orders[[index]]
    data.frame(
      MetricKey = paste0("hill_q", alpha_q_token(q)),
      MetricLabel = sprintf("Hill number (q=%s)", alpha_q_label(q)),
      Parameter = q,
      FigureKey = "02e",
      FacetOrder = index,
      Formula = "exp(Renyi entropy at q)",
      Method = "closed_form",
      ValidDomain = "N > 0 and q >= 0",
      UndefinedRule = "NA when no classified reads",
      stringsAsFactors = FALSE
    )
  }))
  phylo <- data.frame(
    MetricKey = c("faith_pd", "psr", "pse"),
    MetricLabel = c("Faith PD", "PSR", "PSE"),
    Parameter = rep(NA_real_, 3L), FigureKey = rep("02e", 3L),
    FacetOrder = length(hill_orders) + seq_len(3L),
    Formula = c(
      "sum of branches in the rooted subtree connecting observed taxa",
      "PSV * richness", "abundance-weighted phylogenetic species evenness"
    ),
    Method = rep("ape VCV; equations equivalent to picante 1.8.2", 3L),
    ValidDomain = c("rooted tree with branch lengths; richness > 0",
                    "rooted tree with branch lengths; richness > 1",
                    "rooted tree with branch lengths; richness > 1"),
    UndefinedRule = c("NA when no mapped observed taxon", rep("NA when richness <= 1", 2L)),
    stringsAsFactors = FALSE
  )
  rbind(fixed, renyi, hill, phylo)
}

alpha_metric_result <- function(registry, key, value, valid = is.finite(value), reason = "") {
  row <- registry[registry$MetricKey == key, c("MetricKey", "MetricLabel", "Parameter"), drop = FALSE]
  if (nrow(row) != 1L) stop(sprintf("Unknown alpha metric key '%s'.", key), call. = FALSE)
  if (!isTRUE(valid) || !is.finite(value)) {
    value <- NA_real_
    valid <- FALSE
    if (!nzchar(reason)) reason <- "metric was undefined"
  }
  row$Value <- as.numeric(value)
  row$Valid <- isTRUE(valid)
  row$Reason <- if (isTRUE(valid)) "" else reason
  row
}

safe_vegan_scalar <- function(expr) {
  tryCatch(suppressWarnings(as.numeric(expr)), error = function(e) NA_real_)
}

renyi_entropy <- function(probabilities, q) {
  if (identical(as.numeric(q), 1)) return(-sum(probabilities * log(probabilities)))
  log(sum(probabilities^q)) / (1 - q)
}

calc_alpha_metric_long <- function(counts, hill_orders = c(0, 1, 2), renyi_orders = 1) {
  counts <- validate_alpha_counts(counts)
  registry <- alpha_metric_registry(hill_orders, renyi_orders)
  registry <- registry[!registry$MetricKey %in% c("faith_pd", "psr", "pse"), , drop = FALSE]
  N <- sum(counts)
  S <- length(counts)
  if (N == 0L) {
    rows <- lapply(registry$MetricKey, function(key) {
      alpha_metric_result(registry, key, NA_real_, FALSE, "no classified reads")
    })
    return(do.call(rbind, rows))
  }

  p <- counts / N
  shannon <- vegan::diversity(counts, index = "shannon")
  simpson <- vegan::diversity(counts, index = "simpson")
  invsimpson <- vegan::diversity(counts, index = "invsimpson")
  estimated <- tryCatch(
    suppressWarnings(vegan::estimateR(counts)),
    error = function(e) stats::setNames(rep(NA_real_, 5L),
      c("S.obs", "S.chao1", "se.chao1", "S.ACE", "se.ACE"))
  )
  chao1 <- unname(estimated["S.chao1"])
  ace <- unname(estimated["S.ACE"])
  fisher <- safe_vegan_scalar(vegan::fisher.alpha(counts))

  values <- c(
    richness = S, shannon = shannon, ens = exp(shannon), simpson = simpson,
    invsimpson = invsimpson, fisher_alpha = fisher, chao1 = chao1, ace = ace,
    pielou = if (S > 1L) shannon / log(S) else NA_real_,
    berger_parker = max(counts) / N,
    heip = if (S > 1L) (exp(shannon) - 1) / (S - 1) else NA_real_,
    evar = if (S > 1L) {
      1 - (2 / pi) * atan(mean((log(counts) - mean(log(counts)))^2))
    } else NA_real_,
    mcintosh = if (N > 1L) (N - sqrt(sum(counts^2))) / (N - sqrt(N)) else NA_real_
  )
  for (q in renyi_orders) {
    values[[paste0("renyi_q", alpha_q_token(q))]] <- renyi_entropy(p, q)
  }
  for (q in hill_orders) {
    values[[paste0("hill_q", alpha_q_token(q))]] <- exp(renyi_entropy(p, q))
  }

  rows <- lapply(registry$MetricKey, function(key) {
    reason <- if (key %in% c("pielou", "heip", "evar") && S <= 1L) {
      "richness <= 1"
    } else if (identical(key, "mcintosh") && N <= 1L) {
      "classified depth <= 1"
    } else if (!is.finite(values[[key]])) {
      "estimator returned no finite value"
    } else ""
    alpha_metric_result(registry, key, values[[key]], !nzchar(reason), reason)
  })
  do.call(rbind, rows)
}

calc_alpha_indices <- function(counts, hill_orders = c(0, 1, 2), renyi_orders = 1) {
  long <- calc_alpha_metric_long(counts, hill_orders, renyi_orders)
  full_labels <- c(
    richness = "Observed species richness (S)",
    chao1 = "Chao1 (estimated richness)", shannon = "Shannon (H)",
    ens = "Effective number of species (e^H)",
    simpson = "Simpson's D (1-sum p^2)", invsimpson = "Inverse Simpson",
    pielou = "Pielou's evenness (J)", fisher_alpha = "Fisher's alpha",
    berger_parker = "Berger-Parker dominance", ace = "ACE",
    heip = "Heip's evenness", evar = "Smith-Wilson Evar",
    mcintosh = "McIntosh diversity index"
  )
  legacy_first <- c("richness", "chao1", "shannon", "ens", "simpson", "invsimpson",
                    "pielou", "fisher_alpha", "berger_parker")
  remaining <- setdiff(long$MetricKey, legacy_first)
  ordered <- long[match(c(legacy_first, remaining), long$MetricKey), , drop = FALSE]
  labels <- unname(full_labels[ordered$MetricKey])
  labels[is.na(labels)] <- ordered$MetricLabel[is.na(labels)]
  data.frame(Metric = labels, Value = ordered$Value, stringsAsFactors = FALSE)
}

calc_analytical_rarefaction <- function(counts, n_points = 25) {
  counts <- validate_alpha_counts(counts)
  total_classified <- sum(counts)

  if (total_classified < 1) {
    return(data.frame(depth = integer(0), mean_richness = numeric(0), sd_richness = numeric(0)))
  }

  start_depth <- min(100L, total_classified)
  depth_points <- unique(round(seq(start_depth, total_classified, length.out = n_points)))
  # vegan warns when every taxon count is greater than one. That warning is a
  # sampling-design advisory, not an arithmetic failure for genuine read counts.
  rare_res <- suppressWarnings(vegan::rarefy(counts, sample = depth_points, se = TRUE))

  data.frame(
    depth = depth_points,
    mean_richness = as.numeric(rare_res[1, ]),
    sd_richness = as.numeric(rare_res[2, ])
  )
}

calc_rarefaction_resamples <- function(counts, subsample_depth, n_iterations = 100,
                                       seed = 42, hill_orders = c(0, 1, 2),
                                       renyi_orders = 1, phylogeny = NULL) {
  counts <- validate_alpha_counts(counts)
  if (!is.numeric(subsample_depth) || length(subsample_depth) != 1L ||
      is.na(subsample_depth) || subsample_depth != floor(subsample_depth) ||
      subsample_depth < 1L || subsample_depth > sum(counts)) {
    stop("Invalid rarefaction resampling depth for the supplied counts.", call. = FALSE)
  }
  if (!is.numeric(n_iterations) || length(n_iterations) != 1L || is.na(n_iterations) ||
      n_iterations != floor(n_iterations) || n_iterations < 1L) {
    stop("Rarefaction iteration count must be one positive integer.", call. = FALSE)
  }
  registry <- alpha_metric_registry(hill_orders, renyi_orders)
  non_phylo <- registry[!registry$MetricKey %in% c("faith_pd", "psr", "pse"), , drop = FALSE]
  set.seed(seed)
  count_mat <- matrix(counts, nrow = 1L, dimnames = list(NULL, names(counts)))

  wide_rows <- vector("list", n_iterations)
  long_rows <- vector("list", n_iterations)
  for (iteration in seq_len(n_iterations)) {
    sub <- suppressWarnings(vegan::rrarefy(
      count_mat, as.integer(subsample_depth)
    ))[1, ]
    calculated <- calc_alpha_metric_long(sub, hill_orders, renyi_orders)
    if (!is.null(phylogeny) && isTRUE(phylogeny$enabled)) {
      calculated <- rbind(calculated, calc_phylogenetic_alpha(sub, phylogeny, registry))
    }
    values <- stats::setNames(calculated$Value, calculated$MetricKey)
    legacy_first <- c("richness", "chao1", "shannon", "ens", "simpson", "invsimpson", "pielou")
    available <- registry$MetricKey[registry$MetricKey %in% calculated$MetricKey]
    new_keys <- setdiff(available, legacy_first)
    keys <- c(legacy_first, new_keys)
    wide_rows[[iteration]] <- data.frame(
      iteration = iteration,
      subsample_depth = as.integer(subsample_depth),
      as.list(values[keys]), check.names = FALSE
    )
    calculated$iteration <- iteration
    calculated$subsample_depth <- as.integer(subsample_depth)
    long_rows[[iteration]] <- calculated[, c(
      "iteration", "subsample_depth", "MetricKey", "MetricLabel", "Parameter",
      "Value", "Valid", "Reason"
    )]
  }
  wide <- do.call(rbind, wide_rows)
  attr(wide, "alpha_long") <- do.call(rbind, long_rows)
  wide
}
