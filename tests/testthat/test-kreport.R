# =============================================================================
# Unit Tests: Kraken Report (.kreport) & Taxonomy Resolution
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "dependencies.R"))
source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "kreport.R"))
source(file.path("..", "..", "analysis", "07_kreport_pavian.R"))

test_that("kreport tree builder uses standard rank codes D, K, P, C, O, F, G, S", {
  ranks_template <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_subtilis"
  unclass_lineage <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"

  lineages <- c(unclass_lineage, ranks_template)
  counts <- c(10, 50)

  nodes <- build_kreport_tree(lineages, counts)

  # Standard rank codes: D, K, P, C, O, F, G, S
  expect_equal(nodes$rank_code, c("D", "K", "P", "C", "O", "F", "G", "S"))
  expect_false("D1" %in% nodes$rank_code) # Defect fixed!
  expect_equal(nodes$rank_code[2], "K")   # Kingdom uses K
})

test_that("kreport tree validates clade arithmetic and total read invariants", {
  lin1 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp1"
  lin2 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp2"
  uncl <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"

  lineages <- c(uncl, lin1, lin2)
  counts <- c(20, 30, 50)
  total <- sum(counts)

  nodes <- build_kreport_tree(lineages, counts)

  # Should validate without error
  expect_silent(validate_kreport_tree(nodes, total_reads = total, uncl_reads = 20))

  # Check failure when counts corrupted
  expect_error(
    validate_kreport_tree(nodes, total_reads = total + 10, uncl_reads = 20),
    "Kreport tree validation error"
  )
})

test_that("Real Ambar Ayunda fixture builds valid .kreport and runs offline", {
  ab_path <- file.path("..", "..", "output_AAy", "abundance_table_species.tsv")
  cache_path <- file.path("..", "..", "output_AAy", "taxonomy_cache.json")
  asgn_path <- file.path("..", "..", "output_AAy", "reads_assignments",
                         "AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv")

  skip_if_not(file.exists(ab_path), "Abundance table not found")
  skip_if_not(file.exists(cache_path), "Cache not found")

  tmp <- tempdir()
  out_dir <- file.path(tmp, "test_kreport_out")

  cfg <- get_default_config()
  cfg$config_dir <- normalizePath(file.path("..", ".."), winslash = "/")
  cfg$pipeline_root <- cfg$config_dir
  cfg$input$abundance_table <- ab_path
  cfg$taxonomy$cache <- cache_path
  cfg$input$assignments <- list(AmbarAyunda_minimap2_16S = asgn_path)
  cfg$output$base_dir <- out_dir
  cfg$output$dirs <- list(kreport = file.path(out_dir, "07_Kreport"))

  context <- build_context(cfg)

  res <- run_kreport(context)
  expect_equal(res$status, "completed")

  kreport_file <- file.path(cfg$output$dirs$kreport, "AmbarAyunda_minimap2_16S.kreport")
  expect_true(file.exists(kreport_file))

  lines <- readLines(kreport_file)
  expect_gt(length(lines), 100)

  # Line 1: unclassified
  expect_match(lines[1], "^[0-9.]+\t33500\t33500\tU\t0\tunclassified")
  # Line 2: root
  expect_match(lines[2], "^[0-9.]+\t80556\t0\tR\t1\troot")

  # Check kingdom row uses K
  expect_true(any(grepl("\tK\t", lines)))
  expect_false(any(grepl("\tD1\t", lines)))

  unresolved <- read.delim(file.path(cfg$output$dirs$kreport,
                                     "unresolved_taxids.tsv"),
                           check.names = FALSE)
  conflicts <- read.delim(file.path(cfg$output$dirs$kreport,
                                    "taxonomy_conflicts.tsv"),
                          check.names = FALSE)
  expect_equal(nrow(unresolved), 46L)
  expect_equal(nrow(conflicts), 26L)

  resolution <- read.delim(file.path(cfg$output$dirs$kreport,
                                     "taxonomy_resolution.tsv"),
                           check.names = FALSE)
  expect_true("ResolutionSource" %in% names(resolution))
  expect_true(all(resolution$ResolutionSource %in%
                    c("source_cache", "assignment", "ncbi_refresh", "unresolved")))
})

test_that("kreport resolver handles input and output paths containing spaces", {
  root <- tempfile("wf16s path with spaces ")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  lineage <- paste(c("Bacteria", "Bacillati", "Bacillota", "Bacilli",
                     "Bacillales", "Bacillaceae", "Bacillus",
                     "Bacillus subtilis"), collapse = ";")
  abundance <- file.path(root, "abundance table.tsv")
  writeLines(c(
    "tax\tS1\ttotal",
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t1\t1",
    paste(lineage, "2", "2", sep = "\t")
  ), abundance)

  assignments <- file.path(root, "read assignments.tsv")
  writeLines(c(
    "C\tread1\t1423\t0|1500\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis",
    "C\tread2\t1423\t1501\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis",
    "U\tread3\t0\t1490\tUnclassified"
  ), assignments)

  parts <- strsplit(lineage, ";", fixed = TRUE)[[1]]
  cache <- setNames(as.list(seq_len(7L)),
                    vapply(seq_len(7L), function(i) paste(parts[seq_len(i)], collapse = ";"),
                           character(1)))
  cache[[lineage]] <- 0L
  cache_file <- file.path(root, "taxonomy cache.json")
  jsonlite::write_json(cache, cache_file, auto_unbox = TRUE)

  cfg <- get_default_config()
  cfg$pipeline_root <- normalizePath(file.path("..", ".."), winslash = "/")
  cfg$input$abundance_table <- abundance
  cfg$input$assignments <- list(S1 = assignments)
  cfg$taxonomy$cache <- cache_file
  cfg$output$base_dir <- file.path(root, "output directory")
  cfg$output$dirs <- list(kreport = file.path(cfg$output$base_dir, "07_Kreport"))
  cfg$cli <- list(modules = "kreport", validate_only = FALSE)

  context <- build_context(cfg)
  expect_equal(run_kreport(context)$status, "completed")
  expect_true(file.exists(file.path(cfg$output$dirs$kreport, "S1.kreport")))
})
