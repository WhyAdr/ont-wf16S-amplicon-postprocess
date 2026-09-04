# =============================================================================
# Mathematical & Alpha Diversity Metrics Utility
# =============================================================================

suppressMessages({
  library(vegan)
  library(dplyr)
})

calc_alpha_indices <- function(counts) {
  counts <- as.integer(round(counts[counts > 0]))
  total_classified <- sum(counts)
  S <- length(counts)

  if (S == 0 || total_classified == 0) {
    return(data.frame(
      Metric = c("Observed species richness (S)", "Chao1 (estimated richness)",
                 "Shannon (H)", "Effective number of species (e^H)",
                 "Simpson's D (1-sum p^2)", "Inverse Simpson", "Pielou's evenness (J)",
                 "Fisher's alpha", "Berger-Parker dominance"),
      Value = rep(NA_real_, 9),
      stringsAsFactors = FALSE
    ))
  }

  shannon <- vegan::diversity(counts, index = "shannon")
  simpson <- vegan::diversity(counts, index = "simpson")
  invsimpson <- vegan::diversity(counts, index = "invsimpson")
  ens <- exp(shannon)
  pielou <- if (S > 1) shannon / log(S) else NA_real_
  bp_dominance <- max(counts) / total_classified

  chao1 <- tryCatch({
    as.numeric(vegan::estimateR(counts)["S.chao1"])
  }, error = function(e) NA_real_)

  fisher_alpha <- tryCatch({
    as.numeric(vegan::fisher.alpha(counts))
  }, error = function(e) NA_real_)

  data.frame(
    Metric = c("Observed species richness (S)", "Chao1 (estimated richness)",
               "Shannon (H)", "Effective number of species (e^H)",
               "Simpson's D (1-sum p^2)", "Inverse Simpson", "Pielou's evenness (J)",
               "Fisher's alpha", "Berger-Parker dominance"),
    Value = c(S, chao1, shannon, ens, simpson, invsimpson, pielou, fisher_alpha, bp_dominance),
    stringsAsFactors = FALSE
  )
}

calc_analytical_rarefaction <- function(counts, n_points = 25) {
  counts <- as.integer(round(counts[counts > 0]))
  total_classified <- sum(counts)

  if (total_classified < 1) {
    return(data.frame(depth = integer(0), mean_richness = numeric(0), sd_richness = numeric(0)))
  }

  start_depth <- min(100L, total_classified)
  depth_points <- unique(round(seq(start_depth, total_classified, length.out = n_points)))

  rare_res <- vegan::rarefy(counts, sample = depth_points, se = TRUE)

  data.frame(
    depth = depth_points,
    mean_richness = as.numeric(rare_res[1, ]),
    sd_richness = as.numeric(rare_res[2, ])
  )
}

calc_rarefaction_resamples <- function(counts, subsample_depth, n_iterations = 100, seed = 42) {
  counts <- as.integer(round(counts[counts > 0]))
  if (length(counts) == 0L || subsample_depth < 1L || subsample_depth > sum(counts)) {
    stop("Invalid rarefaction resampling depth for the supplied counts.", call. = FALSE)
  }
  set.seed(seed)

  count_mat <- matrix(counts, nrow = 1)

  res_list <- vector("list", n_iterations)
  for (i in seq_len(n_iterations)) {
    sub <- vegan::rrarefy(count_mat, subsample_depth)[1, ]
    sub_counts <- sub[sub > 0]
    S_sub <- length(sub_counts)
    shannon_sub <- vegan::diversity(sub_counts, index = "shannon")
    simpson_sub <- vegan::diversity(sub_counts, index = "simpson")
    invsimpson_sub <- vegan::diversity(sub_counts, index = "invsimpson")
    ens_sub <- exp(shannon_sub)
    pielou_sub <- if (S_sub > 1) shannon_sub / log(S_sub) else NA_real_
    chao1_sub <- tryCatch(as.numeric(vegan::estimateR(sub_counts)["S.chao1"]), error = function(e) NA_real_)

    res_list[[i]] <- data.frame(
      iteration = i,
      subsample_depth = subsample_depth,
      richness = S_sub,
      chao1 = chao1_sub,
      shannon = shannon_sub,
      ens = ens_sub,
      simpson = simpson_sub,
      invsimpson = invsimpson_sub,
      pielou = pielou_sub
    )
  }

  do.call(rbind, res_list)
}
