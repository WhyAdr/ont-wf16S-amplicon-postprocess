#!/usr/bin/env Rscript

# Manual performance harness for the bounded assignment parser.  It deliberately
# writes its synthetic input to a temporary directory and does not form part of
# the automated test suite.

file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_argument) == 1L) {
  sub("^--file=", "", file_argument)
} else {
  "tests/benchmark_assignment_parser.R"
}
script_path <- normalizePath(script_path, winslash = "/", mustWork = FALSE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
source(file.path(repo_root, "analysis", "utils", "config.R"))
source(file.path(repo_root, "analysis", "utils", "io.R"))

arguments <- commandArgs(trailingOnly = TRUE)
row_argument <- grep("^--rows=", arguments, value = TRUE)
row_count <- if (length(row_argument) == 1L) {
  as.integer(sub("^--rows=", "", row_argument))
} else {
  250000L
}
if (is.na(row_count) || row_count < 1L) {
  stop("Use a positive integer, for example --rows=250000.", call. = FALSE)
}

peak_memory_mb <- function() {
  if (file.exists("/proc/self/status")) {
    status <- readLines("/proc/self/status", warn = FALSE)
    hwm <- grep("^VmHWM:", status, value = TRUE)
    if (length(hwm) == 1L) {
      return(as.numeric(sub("^VmHWM:\\s+([0-9]+)\\s+kB$", "\\1", hwm)) / 1024)
    }
  }
  if (.Platform$OS.type == "windows" && getRversion() < "4.4.0" &&
      exists("memory.size", where = asNamespace("utils"))) {
    peak <- suppressWarnings(utils::memory.size(max = TRUE))
    if (is.finite(peak)) return(peak)
  }
  if (.Platform$OS.type == "windows" && requireNamespace("ps", quietly = TRUE)) {
    memory <- ps::ps_memory_info()
    if ("peak_wset" %in% names(memory) && is.finite(memory[["peak_wset"]])) {
      return(unname(memory[["peak_wset"]]) / 1024^2)
    }
  }
  NA_real_
}

benchmark_root <- tempfile("assignment_parser_benchmark_")
dir.create(benchmark_root)
on.exit(unlink(benchmark_root, recursive = TRUE, force = TRUE), add = TRUE)
assignment_path <- file.path(benchmark_root, "synthetic_assignments.tsv")
connection <- file(assignment_path, open = "wt")

chunk_size <- 100000L
for (start in seq.int(1L, row_count, by = chunk_size)) {
  end <- min(row_count, start + chunk_size - 1L)
  indices <- start:end
  classified <- indices %% 5L != 0L
  lines <- sprintf(
    "%s\tbenchmark_%09d\t%d\t0|%d\t%s",
    ifelse(classified, "C", "U"), indices, ifelse(classified, 123L, 0L),
    1400L + indices %% 200L, ifelse(classified, "Bacteria|Example", "Unclassified")
  )
  writeLines(lines, connection)
}
close(connection)

invisible(gc(reset = TRUE))
started <- proc.time()[["elapsed"]]
reads <- read_assignments_file(assignment_path, "benchmark", chunk_size = chunk_size)
elapsed_seconds <- proc.time()[["elapsed"]] - started
peak_memory <- peak_memory_mb()

results <- list(
  rows = nrow(reads),
  elapsed_seconds = unname(elapsed_seconds),
  final_object_size_mb = unname(as.numeric(object.size(reads)) / 1024^2),
  peak_memory_mb = peak_memory,
  peak_memory_method = if (file.exists("/proc/self/status")) {
    "Linux /proc/self/status VmHWM"
  } else if (.Platform$OS.type == "windows" && is.finite(peak_memory)) {
    "Windows ps::ps_memory_info() peak_wset"
  } else {
    "unavailable on this host"
  }
)
print(results)
