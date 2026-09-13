# =============================================================================
# Kraken Report (.kreport) Formatter and Tree Validator
# =============================================================================

suppressMessages(library(jsonlite))

RANKS_8 <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
RANK_CODES_8 <- c("D", "K", "P", "C", "O", "F", "G", "S")

kreport_sort_key <- function(value) {
  paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(value)))), collapse = "")
}

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
    children <- children[order(
      -children$reads_clade,
      vapply(children$name, kreport_sort_key, character(1)),
      method = "radix"
    ), ]
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

validate_krona_count <- function(value, name) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= 0 &&
    abs(value - round(value)) <= sqrt(.Machine$double.eps)
  if (!valid) {
    stop(sprintf("Krona %s must be a finite, non-negative integer.", name), call. = FALSE)
  }
  as.numeric(value)
}

validate_krona_count_vector <- function(values, name) {
  if (!is.numeric(values) || any(!is.finite(values)) || any(is.na(values)) ||
      any(values < 0) || any(abs(values - round(values)) > sqrt(.Machine$double.eps))) {
    stop(sprintf("Krona %s must contain only finite, non-negative integers.", name),
         call. = FALSE)
  }
  as.numeric(values)
}

normalize_krona_path <- function(path, index) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop(sprintf("Krona taxonomy path at row %d is empty or invalid.", index), call. = FALSE)
  }

  labels <- strsplit(path, ";", fixed = TRUE)[[1]]
  if (length(labels) == 0L || length(labels) > length(RANKS_8)) {
    stop(sprintf("Krona taxonomy path at row %d has an invalid rank depth.", index),
         call. = FALSE)
  }

  labels <- vapply(labels, function(label) {
    if (!is.na(label) && grepl("[\t\r\n]", label, perl = TRUE)) {
      stop(sprintf("Krona taxonomy label at row %d contains a tab or newline.", index),
           call. = FALSE)
    }
    if (is.na(label) || !nzchar(trimws(label))) "Unknown" else label
  }, character(1))

  paste(labels, collapse = "\t")
}

validate_kreport_label <- function(label, index, context = "Kreport") {
  if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label)) {
    stop(sprintf("%s taxonomy label at row %d is empty or invalid.", context, index),
         call. = FALSE)
  }
  if (grepl("[[:cntrl:]]", label, perl = TRUE) || grepl("^[ ]", label, perl = TRUE)) {
    stop(sprintf(
      "%s taxonomy label at row %d contains a control character or ambiguous leading indentation.",
      context, index), call. = FALSE)
  }
  label
}

format_krona_magnitude <- function(value) {
  trimws(formatC(value, format = "f", digits = 0))
}

build_krona_lines <- function(nodes_df, total_reads, uncl_reads) {
  total_num <- validate_krona_count(total_reads, "total_reads")
  uncl_num <- validate_krona_count(uncl_reads, "uncl_reads")
  if (total_num <= 0) {
    stop("Krona total_reads must be greater than zero.", call. = FALSE)
  }
  if (uncl_num > total_num) {
    stop("Krona unclassified reads cannot exceed total_reads.", call. = FALSE)
  }

  if (!is.data.frame(nodes_df) || !all(c("path", "reads_taxon") %in% names(nodes_df))) {
    stop("Krona nodes_df must contain 'path' and 'reads_taxon' columns.", call. = FALSE)
  }

  direct_counts <- validate_krona_count_vector(nodes_df$reads_taxon, "nodes_df$reads_taxon")
  classified_num <- total_num - uncl_num
  direct_sum <- sum(direct_counts)
  if (!isTRUE(direct_sum == classified_num)) {
    stop(sprintf(
      "Krona classified direct counts (%s) do not equal ClassifiedReads (%s).",
      format_krona_magnitude(direct_sum), format_krona_magnitude(classified_num)
    ), call. = FALSE)
  }

  positive_rows <- which(direct_counts > 0)
  lines <- character(0)
  if (uncl_num > 0) {
    lines <- paste(format_krona_magnitude(uncl_num), "Unclassified", sep = "\t")
  }

  if (length(positive_rows) > 0L) {
    classified_lines <- vapply(positive_rows, function(row_index) {
      path <- normalize_krona_path(nodes_df$path[row_index], row_index)
      paste(format_krona_magnitude(direct_counts[row_index]), path, sep = "\t")
    }, character(1))
    lines <- c(lines, classified_lines)
  }

  emitted_sum <- uncl_num + sum(direct_counts[positive_rows])
  if (!isTRUE(emitted_sum == total_num)) {
    stop(sprintf(
      "Krona emitted magnitude sum (%s) does not equal total_reads (%s).",
      format_krona_magnitude(emitted_sum), format_krona_magnitude(total_num)
    ), call. = FALSE)
  }

  if (length(lines) == 0L) {
    stop("Krona export produced no data lines.", call. = FALSE)
  }
  lines
}

write_krona_input <- function(output_path, nodes_df, total_reads, uncl_reads) {
  if (!is.character(output_path) || length(output_path) != 1L ||
      is.na(output_path) || !nzchar(trimws(output_path))) {
    stop("Krona output_path must be one non-empty path.", call. = FALSE)
  }

  parent_dir <- dirname(output_path)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent_dir)) {
    stop(sprintf("Could not create Krona output directory '%s'.", parent_dir), call. = FALSE)
  }

  lines <- build_krona_lines(nodes_df, total_reads, uncl_reads)
  writeLines(lines, output_path)
  if (!file.exists(output_path) || !isTRUE(file.info(output_path)$size > 0)) {
    stop(sprintf("Krona input was not written or is empty: '%s'.", output_path), call. = FALSE)
  }
  invisible(output_path)
}

render_kronatools_html <- function(executable, output_path, sample_id, input_path) {
  valid_executable <- is.character(executable) && length(executable) == 1L &&
    !is.na(executable) && nzchar(trimws(executable))
  valid_sample <- is.character(sample_id) && length(sample_id) == 1L &&
    !is.na(sample_id) && nzchar(sample_id)
  if (!valid_executable || !valid_sample) {
    stop("Krona HTML rendering requires a non-empty executable and sample ID.", call. = FALSE)
  }
  if (!file.exists(input_path) || !isTRUE(file.info(input_path)$size > 0)) {
    stop(sprintf("Krona input for sample '%s' is missing or empty.", sample_id), call. = FALSE)
  }

  parent_dir <- dirname(output_path)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent_dir)) {
    stop(sprintf("Could not create Krona HTML output directory '%s'.", parent_dir), call. = FALSE)
  }

  rendered_path <- tempfile(pattern = paste0(".", basename(output_path), ".tmp-"),
                            tmpdir = parent_dir)
  on.exit(if (file.exists(rendered_path)) unlink(rendered_path, force = TRUE), add = TRUE)
  args <- c("-o", rendered_path, "-n", sample_id, input_path)
  result <- tryCatch(
    processx::run(
      command = executable,
      args = args,
      error_on_status = FALSE
    ),
    error = function(e) {
      stop(sprintf(
        "Krona HTML rendering failed for sample '%s' using executable '%s': %s",
        sample_id, executable, e$message
      ), call. = FALSE)
    }
  )

  if (!isTRUE(result$status == 0L)) {
    status <- if (is.null(result$status)) "unknown" else as.character(result$status)
    detail <- trimws(paste(result$stderr %||% "", result$stdout %||% ""))
    stop(sprintf(
      "Krona HTML rendering failed for sample '%s' using executable '%s' (exit status %s): %s",
      sample_id, executable, status, detail
    ), call. = FALSE)
  }
  if (!file.exists(rendered_path) || !isTRUE(file.info(rendered_path)$size > 0)) {
    stop(sprintf(
      "Krona HTML renderer reported success but produced no output for sample '%s' using executable '%s'.",
      sample_id, executable
    ), call. = FALSE)
  }

  if (!exists("atomic_replace", mode = "function")) {
    stop("Krona HTML rendering requires the atomic output helper.", call. = FALSE)
  }
  atomic_replace(output_path, function(temp) {
    if (!file.copy(rendered_path, temp, overwrite = TRUE)) {
      stop(sprintf("Could not stage Krona HTML output for sample '%s'.", sample_id), call. = FALSE)
    }
  })

  invisible(output_path)
}

krona_vendor_directory <- function(repo_root) {
  file.path(repo_root, "analysis", "vendor", "krona-2.8.1")
}

krona_artifact_relpath <- function(path, module_root) {
  canonical <- function(value) {
    if (exists("canonicalize_root_path", mode = "function")) {
      canonicalize_root_path(value)
    } else {
      normalizePath(value, winslash = "/", mustWork = TRUE)
    }
  }
  normalized_path <- canonical(path)
  normalized_root <- canonical(module_root)
  path_key <- if (identical(.Platform$OS.type, "windows")) tolower(normalized_path) else normalized_path
  root_key <- if (identical(.Platform$OS.type, "windows")) tolower(normalized_root) else normalized_root
  if (!startsWith(path_key, paste0(root_key, "/"))) {
    stop("Krona artifact is outside the module output root.", call. = FALSE)
  }
  relative <- substring(normalized_path, nchar(normalized_root) + 2L)
  parts <- strsplit(relative, "/", fixed = TRUE)[[1]]
  if (!nzchar(relative) || any(!nzchar(parts)) || any(parts %in% c(".", "..")) ||
      grepl("^[A-Za-z]:|^/", relative)) {
    stop("Krona artifact has an unsafe run-relative path.", call. = FALSE)
  }
  relative
}

krona_vendor_manifest <- function(vendor_dir) {
  if (nzchar(Sys.readlink(vendor_dir)) || !dir.exists(vendor_dir)) {
    stop(sprintf("Krona vendor directory is missing: '%s'.", vendor_dir), call. = FALSE)
  }
  manifest_path <- file.path(vendor_dir, "SOURCE.json")
  if (!file.exists(manifest_path) || dir.exists(manifest_path)) {
    stop(sprintf("Krona vendor manifest is missing: '%s'.", manifest_path), call. = FALSE)
  }
  manifest <- tryCatch(
    jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
    error = function(e) stop(sprintf("Krona vendor manifest is invalid: %s", e$message), call. = FALSE)
  )
  if (!identical(manifest$upstream, "marbl/Krona") || !identical(manifest$tag, "v2.8.1")) {
    stop("Krona vendor manifest does not identify marbl/Krona v2.8.1.", call. = FALSE)
  }
  if (!is.character(manifest$commit) || length(manifest$commit) != 1L ||
      !grepl("^[0-9a-f]{40}$", manifest$commit)) {
    stop("Krona vendor manifest has an invalid upstream commit.", call. = FALSE)
  }
  entries <- manifest$files %||% list()
  if (!is.list(entries) || length(entries) == 0L) {
    stop("Krona vendor manifest has no pinned files.", call. = FALSE)
  }
  expected_files <- c(
    "LICENSE.txt", "src/krona-2.0.js", "img/favicon.ico", "img/hidden.png",
    "img/loading.gif", "img/logo-med.png"
  )
  declared <- character(0)
  declared_casefolded <- character(0)
  manifest_entries <- vapply(entries, function(entry) {
    if (!is.character(entry$path) || length(entry$path) != 1L ||
        !is.character(entry$sha256) || length(entry$sha256) != 1L) {
      stop("Krona vendor manifest contains an incomplete file entry.", call. = FALSE)
    }
    path <- as.character(entry$path)
    expected <- as.character(entry$sha256)
    if (length(path) != 1L || is.na(path) || !nzchar(path) || grepl("\\\\", path) ||
        grepl("^/|^[A-Za-z]:", path)) {
      stop(sprintf("Krona vendor manifest path is unsafe: '%s'.", path), call. = FALSE)
    }
    parts <- strsplit(path, "/", fixed = TRUE)[[1]]
    if (any(!nzchar(parts)) || any(parts %in% c(".", "..")) ||
        !identical(path, paste(parts, collapse = "/"))) {
      stop(sprintf("Krona vendor manifest path is unsafe: '%s'.", path), call. = FALSE)
    }
    folded <- tolower(path)
    if (path %in% declared || folded %in% declared_casefolded) {
      stop(sprintf("Krona vendor manifest contains duplicate or case-colliding path: '%s'.", path), call. = FALSE)
    }
    declared <<- c(declared, path)
    declared_casefolded <<- c(declared_casefolded, folded)
    if (length(expected) != 1L || is.na(expected) || !grepl("^[0-9a-f]{64}$", expected)) {
      stop(sprintf("Krona vendor manifest has an invalid SHA-256 for '%s'.", path), call. = FALSE)
    }
    actual_path <- file.path(vendor_dir, path)
    if (nzchar(Sys.readlink(actual_path)) || !file.exists(actual_path) || dir.exists(actual_path)) {
      stop(sprintf("Krona vendor file is missing: '%s'.", actual_path), call. = FALSE)
    }
    actual <- compute_file_hash(actual_path)
    if (!identical(actual, expected)) {
      stop(sprintf("Krona vendor SHA-256 mismatch for '%s'.", path), call. = FALSE)
    }
    paste(path, expected, sep = "\t")
  }, character(1))
  all_paths <- list.files(vendor_dir, recursive = TRUE, all.files = TRUE,
                          no.. = TRUE, include.dirs = TRUE, full.names = FALSE)
  symlinks <- all_paths[vapply(file.path(vendor_dir, all_paths), function(path) {
    nzchar(Sys.readlink(path))
  }, logical(1))]
  if (length(symlinks)) {
    stop(sprintf("Krona vendor tree contains a symlink: '%s'.", symlinks[1]), call. = FALSE)
  }
  actual <- all_paths[vapply(all_paths, function(relative) {
    path <- file.path(vendor_dir, relative)
    file.exists(path) && !dir.exists(path) && !identical(gsub("\\\\", "/", relative), "SOURCE.json")
  }, logical(1))]
  actual <- gsub("\\\\", "/", actual)
  if (anyDuplicated(tolower(actual))) {
    stop("Krona vendor tree contains case-colliding files.", call. = FALSE)
  }
  if (!setequal(actual, declared)) {
    stop("Declared and actual Krona vendor inventories differ.", call. = FALSE)
  }
  if (!setequal(declared, expected_files)) {
    stop("Krona vendor inventory differs from the renderer contract.", call. = FALSE)
  }
  list(
    tag = as.character(manifest$tag),
    commit = as.character(manifest$commit),
    source_url = as.character(manifest$source_url),
    entries = manifest_entries,
    manifest_sha256 = compute_file_hash(manifest_path)
  )
}

resolve_krona_renderer <- function(krona_cfg, repo_root) {
  policy <- tolower(as.character(krona_cfg$html_renderer %||% "builtin"))
  if (!policy %in% c("builtin", "kronatools", "auto")) {
    stop(sprintf("Unsupported Krona html_renderer '%s'; use builtin, kronatools, or auto.", policy),
         call. = FALSE)
  }

  requested <- as.character(krona_cfg$executable %||% "ktImportText")
  external <- find_krona_executable(requested)
  vendor_dir <- krona_vendor_directory(repo_root)
  builder <- file.path(repo_root, "analysis", "utils", "krona_builder.py")
  pipeline_version <- read_pipeline_version(file.path(repo_root, "VERSION"))
  builtin <- list(
    policy = policy,
    provider = "builtin",
    renderer = "builtin_krona_compatible",
    renderer_version = pipeline_version,
    krona_version = NULL,
    builder = normalizePath(builder, winslash = "/", mustWork = FALSE),
    vendor_dir = normalizePath(vendor_dir, winslash = "/", mustWork = FALSE),
    requested_executable = requested,
    resolved_executable = NULL,
    vendor_manifest_sha256 = NULL
  )

  if (identical(policy, "builtin") || (identical(policy, "auto") && is.na(external))) {
    if (!file.exists(builder) || dir.exists(builder)) {
      stop(sprintf("Builtin Krona renderer is missing: '%s'.", builder), call. = FALSE)
    }
    vendor <- krona_vendor_manifest(vendor_dir)
    builtin$vendor_manifest_sha256 <- vendor$manifest_sha256
    builtin$krona_version <- sub("^v", "", vendor$tag)
    return(builtin)
  }

  if (is.na(external)) {
    stop(sprintf("KronaTools executable '%s' was not found for html_renderer='%s'.",
                 requested, policy), call. = FALSE)
  }
  list(
    policy = policy,
    provider = "kronatools",
    renderer = "ktImportText",
    renderer_version = get_krona_version(external) %||% "unknown",
    krona_version = get_krona_version(external) %||% "unknown",
    builder = NULL,
    vendor_dir = NULL,
    requested_executable = requested,
    resolved_executable = external,
    vendor_manifest_sha256 = NULL
  )
}

render_builtin_krona_html <- function(python_cmd, builder_path, vendor_dir,
                                      output_path, sample_id, input_path, expected_total) {
  if (!file.exists(builder_path) || dir.exists(builder_path)) {
    stop(sprintf("Builtin Krona renderer is missing: '%s'.", builder_path), call. = FALSE)
  }
  if (!file.exists(input_path) || !isTRUE(file.info(input_path)$size > 0)) {
    stop(sprintf("Krona input for sample '%s' is missing or empty.", sample_id), call. = FALSE)
  }
  parent_dir <- dirname(output_path)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent_dir)) {
    stop(sprintf("Could not create Krona HTML output directory '%s'.", parent_dir), call. = FALSE)
  }
  args <- c(
    builder_path,
    "--input", input_path,
    "--output", output_path,
    "--dataset-name", sample_id,
    "--expected-total", format_krona_magnitude(expected_total),
    "--vendor-dir", vendor_dir
  )
  result <- tryCatch(
    processx::run(command = python_cmd, args = args, error_on_status = FALSE),
    error = function(e) stop(sprintf(
      "Builtin Krona HTML rendering failed for sample '%s': %s", sample_id, e$message
    ), call. = FALSE)
  )
  if (!identical(result$status, 0L)) {
    detail <- trimws(paste(result$stderr, result$stdout))
    stop(sprintf("Builtin Krona HTML rendering failed for sample '%s' (exit status %s): %s",
                 sample_id, as.character(result$status %||% "unknown"), detail), call. = FALSE)
  }
  if (!file.exists(output_path) || !isTRUE(file.info(output_path)$size > 0)) {
    stop(sprintf("Builtin Krona renderer reported success but produced no output for sample '%s'.",
                 sample_id), call. = FALSE)
  }
  temporary <- list.files(parent_dir, full.names = TRUE)
  temporary <- temporary[startsWith(
    basename(temporary), paste0(basename(output_path), ".tmp-"))]
  if (length(temporary) > 0L) {
    stop(sprintf("Builtin Krona renderer left temporary output(s) for sample '%s'.", sample_id),
         call. = FALSE)
  }
  invisible(output_path)
}

format_kreport_lines <- function(nodes_sorted, total_reads, uncl_reads, taxid_cache = list()) {
  if (!is.numeric(total_reads) || length(total_reads) != 1L || !is.finite(total_reads) ||
      total_reads <= 0 || total_reads != round(total_reads)) {
    stop("Kreport total_reads must be a positive finite integer.", call. = FALSE)
  }
  if (!is.numeric(uncl_reads) || length(uncl_reads) != 1L || !is.finite(uncl_reads) ||
      uncl_reads < 0 || uncl_reads != round(uncl_reads) || uncl_reads > total_reads) {
    stop("Kreport unclassified reads must be a finite non-negative integer not exceeding total_reads.",
         call. = FALSE)
  }
  cl_reads <- total_reads - uncl_reads

  lines <- character(nrow(nodes_sorted) + 2)

  # Line 1: unclassified
  lines[1] <- sprintf("%.2f\t%.0f\t%.0f\tU\t0\tunclassified",
                      100 * uncl_reads / total_reads, uncl_reads, uncl_reads)

  # Line 2: root
  lines[2] <- sprintf("%.2f\t%.0f\t%.0f\tR\t1\troot",
                      100 * cl_reads / total_reads, cl_reads, 0L)

  for (i in seq_len(nrow(nodes_sorted))) {
    validate_kreport_label(nodes_sorted$name[i], i + 2L)
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
