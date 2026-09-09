# =============================================================================
# Unit Tests: Kraken Report (.kreport) & Taxonomy Resolution
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "dependencies.R"))
source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "version.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "manifest.R"))
source(file.path("..", "..", "analysis", "utils", "atomic_io.R"))
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

test_that("Krona lines use direct terminal counts without ancestor clade duplication", {
  lin1 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp1"
  lin2 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp2"
  uncl <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
  nodes <- build_kreport_tree(c(uncl, lin1, lin2), c(20, 30, 50))

  lines <- build_krona_lines(nodes, total_reads = 100, uncl_reads = 20)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  magnitudes <- as.numeric(vapply(fields, function(field) field[[1]], character(1)))

  expect_length(lines, 3L)
  expect_equal(fields[[1]], c("20", "Unclassified"))
  expect_equal(sum(magnitudes), 100)
  expect_equal(sum(magnitudes[-1]), 80)
  expect_true(all(vapply(fields[-1], function(field) length(field) == 9L, logical(1))))
  expect_equal(
    sort(vapply(fields[-1], function(field) field[[9]], character(1))),
    c("Bacillus_sp1", "Bacillus_sp2")
  )
})

test_that("Krona lines omit zero unclassified contributions and preserve exact sums", {
  lin1 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp1"
  lin2 <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus_sp2;Bacillus_sp2"
  nodes <- build_kreport_tree(c(lin1, lin2), c(30, 70))

  lines <- build_krona_lines(nodes, total_reads = 100, uncl_reads = 0)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  magnitudes <- as.numeric(vapply(fields, function(field) field[[1]], character(1)))

  expect_length(lines, 2L)
  expect_false(any(vapply(fields, function(field) identical(field[[2]], "Unclassified"), logical(1))))
  expect_equal(sum(magnitudes), 100)
})

test_that("Krona paths reject physical delimiters and normalize empty labels", {
  lineage <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp1"
  nodes <- build_kreport_tree(lineage, 10)

  bad_tab <- nodes
  bad_tab$path[8] <- "Bacteria;bad\tlabel"
  expect_error(build_krona_lines(bad_tab, 10, 0), "tab or newline")

  bad_newline <- nodes
  bad_newline$path[8] <- "Bacteria;bad\nlabel"
  expect_error(build_krona_lines(bad_newline, 10, 0), "tab or newline")

  empty_label <- nodes
  empty_label$path[8] <- "Bacteria;;Bacillota"
  lines <- build_krona_lines(empty_label, 10, 0)
  expect_true(grepl("Bacteria\tUnknown\tBacillota$", lines[1]))
})

test_that("Krona arithmetic rejects invalid and mismatched counts", {
  lineage <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp1"
  nodes <- build_kreport_tree(lineage, 10)

  expect_error(build_krona_lines(nodes, total_reads = 0, uncl_reads = 0), "greater than zero")
  expect_error(build_krona_lines(nodes, total_reads = 10.5, uncl_reads = 0), "finite")
  expect_error(build_krona_lines(nodes, total_reads = 10, uncl_reads = -1), "finite")
  expect_error(build_krona_lines(nodes, total_reads = 10, uncl_reads = 1), "classified direct")

  non_integer <- nodes
  non_integer$reads_taxon[8] <- 1.5
  expect_error(build_krona_lines(non_integer, 10, 0), "finite")

  negative <- nodes
  negative$reads_taxon[8] <- -1
  expect_error(build_krona_lines(negative, 10, 0), "finite")
})

test_that("write_krona_input creates a physical TSV under paths with spaces", {
  lineage <- "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_sp1"
  nodes <- build_kreport_tree(lineage, 10)
  output <- file.path(tempdir(), "krona path with spaces", "sample id.krona.tsv")

  expect_equal(write_krona_input(output, nodes, 10, 0), output)
  expect_true(file.exists(output))
  expect_equal(readLines(output, warn = FALSE),
               "10\tBacteria\tBacillati\tBacillota\tBacilli\tBacillales\tBacillaceae\tBacillus\tBacillus_sp1")
})

test_that("builtin Krona renderer is deterministic and strict external policy fails closed", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/")
  renderer <- resolve_krona_renderer(
    list(enabled = TRUE, render_html = TRUE, html_renderer = "builtin",
         executable = "ktImportText"),
    repo_root
  )
  expect_equal(renderer$provider, "builtin")
  expect_true(file.exists(renderer$builder))
  expect_true(file.exists(file.path(renderer$vendor_dir, "SOURCE.json")))

  root <- tempfile("builtin krona path with spaces ")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  input <- file.path(root, "sample input.tsv")
  output <- file.path(root, "sample output.html")
  writeLines(c("2\tBacteria\tFirmicutes", "1\tBacteria\tProteobacteria"), input)
  render_builtin_krona_html(
    find_python(), renderer$builder, renderer$vendor_dir, output,
    "sample id", input, 3
  )
  first <- readBin(output, "raw", n = file.info(output)$size)
  render_builtin_krona_html(
    find_python(), renderer$builder, renderer$vendor_dir, output,
    "sample id", input, 3
  )
  expect_identical(first, readBin(output, "raw", n = file.info(output)$size))
  expect_true(grepl("id=\"hiddenImage\" src=\"data:image/png", rawToChar(first), fixed = TRUE))
  expect_length(list.files(root, pattern = "[.]tmp-", full.names = TRUE), 0L)

  expect_error(
    resolve_krona_renderer(
      list(enabled = TRUE, render_html = TRUE, html_renderer = "kronatools",
           executable = file.path(root, "missing-ktImportText")),
      repo_root
    ),
    "was not found"
  )
})

test_that("external Krona renderer failures preserve prior output and clean partial HTML", {
  root <- tempfile("krona_renderer_failure_")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  input <- file.path(root, "input.tsv")
  output <- file.path(root, "output.html")
  writeLines("1\tBacteria", input)
  writeLines("prior HTML", output)
  executable <- if (identical(.Platform$OS.type, "windows")) {
    script <- file.path(root, "failing-renderer.cmd")
    writeLines(c(
      "@echo off",
      "echo renderer stderr 1>&2",
      "echo partial>\"%~2\"",
      "exit /b 7"
    ), script)
    script
  } else {
    script <- file.path(root, "failing-renderer.sh")
    writeLines(c(
      "#!/bin/sh",
      "printf partial > \"$2\"",
      "printf 'renderer stderr\\n' >&2",
      "exit 7"
    ), script)
    Sys.chmod(script, mode = "0755")
    script
  }

  expect_error(
    render_kronatools_html(executable, output, "S1", input),
    "renderer stderr"
  )
  expect_identical(readLines(output, warn = FALSE), "prior HTML")
  expect_length(list.files(root, pattern = "[.]tmp-", full.names = TRUE), 0L)
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
  cfg$input$params_json <- file.path("..", "..", "output_AAy", "params.json")
  cfg$taxonomy$cache <- cache_path
  cfg$input$assignments <- list(AmbarAyunda_minimap2_16S = asgn_path)
  cfg$output$base_dir <- out_dir
  cfg$output$dirs <- list(kreport = file.path(out_dir, "07_Kreport"))
  cfg$krona <- list(enabled = TRUE, render_html = FALSE, executable = "ktImportText")

  context <- build_context(cfg)

  expect_warning(
    res <- run_kreport(context),
    "26 lineage-to-TaxID conflict"
  )
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

  krona_file <- file.path(cfg$output$dirs$kreport, "krona",
                          "AmbarAyunda_minimap2_16S.krona.tsv")
  krona_provenance_file <- file.path(cfg$output$dirs$kreport, "krona",
                                     "krona_provenance.json")
  expect_true(file.exists(krona_file))
  expect_true(file.exists(krona_provenance_file))
  krona_magnitudes <- as.numeric(vapply(
    strsplit(readLines(krona_file, warn = FALSE), "\t", fixed = TRUE),
    function(field) field[[1]], character(1)
  ))
  expect_equal(sum(krona_magnitudes), 114056)
  expect_equal(sum(krona_magnitudes[-1]), 80556)
  krona_provenance <- jsonlite::fromJSON(krona_provenance_file, simplifyVector = FALSE)
  expect_equal(krona_provenance$schema_version, 1)
  expect_equal(krona_provenance$path_basis, "run_dir")
  expect_equal(krona_provenance$html_status, "not_requested")
  expect_equal(krona_provenance$samples[[1]]$emitted_magnitude_sum, 114056)
  expect_equal(krona_provenance$samples[[1]]$tsv_path,
               "07_Kreport/krona/AmbarAyunda_minimap2_16S.krona.tsv")
  expect_match(krona_provenance$samples[[1]]$tsv_sha256, "^[0-9a-f]{64}$")
  expect_null(krona_provenance$samples[[1]]$html_path)
  expect_null(krona_provenance$samples[[1]]$html_sha256)

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
                    c("source_cache", "assignment", "assignment_conflict",
                      "ncbi_refresh", "unresolved")))
  expect_gt(sum(resolution$Status == "Conflicted"), 0L)
  expect_true(all(resolution$Status[resolution$ResolutionSource == "assignment_conflict"] ==
                    "Conflicted"))
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
  cfg$input$params_json <- create_temp_params(root)
  cfg$input$assignments <- list(S1 = assignments)
  cfg$taxonomy$cache <- cache_file
  run_stage <- file.path(root, "run stage")
  dir.create(run_stage, recursive = TRUE, showWarnings = FALSE)
  module_stage <- prepare_module_staging(run_stage, "kreport")
  cfg$output$base_dir <- module_stage
  cfg$output$dirs <- list(kreport = file.path(module_stage, "07_Kreport"))
  cfg$krona <- list(enabled = TRUE, render_html = FALSE, executable = "ktImportText")
  cfg$cli <- list(modules = "kreport", validate_only = FALSE)

  context <- build_context(cfg)
  result <- run_kreport(context)
  expect_equal(result$status, "completed")
  expect_true(file.exists(file.path(module_stage, "07_Kreport", "S1.kreport")))
  publish_module_staging(module_stage, run_stage, result$outputs)
  expect_false(dir.exists(module_stage))
  expect_true(file.exists(file.path(run_stage, "07_Kreport", "S1.kreport")))
  provenance <- jsonlite::fromJSON(
    file.path(run_stage, "07_Kreport", "krona", "krona_provenance.json"),
    simplifyVector = FALSE
  )
  record <- provenance$samples[[1]]
  expect_equal(record$tsv_path, "07_Kreport/krona/S1.krona.tsv")
  expect_match(record$tsv_sha256, "^[0-9a-f]{64}$")
  expect_false(grepl("^([A-Za-z]:|/)|[.]module-stage", record$tsv_path))
  expect_identical(compute_file_hash(file.path(run_stage, record$tsv_path)), record$tsv_sha256)

  moved <- file.path(root, "moved completed run")
  dir.create(moved, recursive = TRUE, showWarnings = FALSE)
  expect_true(all(file.copy(list.files(run_stage, full.names = TRUE), moved, recursive = TRUE)))
  moved_record <- jsonlite::fromJSON(
    file.path(moved, "07_Kreport", "krona", "krona_provenance.json"),
    simplifyVector = FALSE
  )$samples[[1]]
  expect_identical(compute_file_hash(file.path(moved, moved_record$tsv_path)), moved_record$tsv_sha256)
})
