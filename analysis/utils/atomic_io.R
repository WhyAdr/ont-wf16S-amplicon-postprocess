# =============================================================================
# Atomic output helpers and whole-run publication
# =============================================================================

utc_timestamp_now <- function() {
  strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

atomic_replace <- function(path, writer, expected_sha256 = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  target_dir <- dirname(path)
  target_name <- basename(path)
  temp <- tempfile(pattern = paste0(".", target_name, ".tmp-"), tmpdir = target_dir)
  backup <- tempfile(pattern = paste0(".", target_name, ".bak-"), tmpdir = target_dir)
  had_target <- file.exists(path)
  orig_mode <- if (had_target) file.info(path)$mode else NULL
  replacement_verified <- FALSE

  restore_backup <- function() {
    if (!had_target || !file.exists(backup)) return(invisible(TRUE))
    if (file.exists(path) && unlink(path, force = TRUE) != 0L) {
      stop(sprintf("Could not remove failed replacement for '%s'; backup retained at '%s'.",
                   path, backup), call. = FALSE)
    }
    if (!file.rename(backup, path)) {
      stop(sprintf("Could not restore backup for '%s'; backup retained at '%s'.",
                   path, backup), call. = FALSE)
    }
    invisible(TRUE)
  }

  on.exit({
    if (file.exists(temp)) unlink(temp, force = TRUE)
    if (file.exists(backup)) {
      if (!replacement_verified) {
        tryCatch(restore_backup(), error = function(e) {
          warning(conditionMessage(e), call. = FALSE)
        })
      } else if (unlink(backup, force = TRUE) != 0L || file.exists(backup)) {
        warning(sprintf("Verified replacement retained an undeleted backup at '%s'.", backup),
                call. = FALSE)
      }
    }
  }, add = TRUE)

  writer(temp)
  if (!file.exists(temp) || dir.exists(temp)) {
    stop(sprintf("Atomic writer did not create '%s'.", path), call. = FALSE)
  }

  temp_size <- file.info(temp)$size
  hash_file <- function(candidate) {
    if (exists("compute_file_hash", mode = "function")) {
      compute_file_hash(candidate)
    } else {
      digest::digest(candidate, file = TRUE, algo = "sha256")
    }
  }
  temp_hash <- hash_file(temp)

  if (!is.null(expected_sha256) && !is.null(temp_hash) &&
      !identical(tolower(temp_hash), tolower(expected_sha256))) {
    stop(sprintf("Replacement '%s' failed expected SHA-256 verification.", path), call. = FALSE)
  }

  if (had_target) {
    if (!file.rename(path, backup)) {
      stop(sprintf("Could not stage backup for '%s'.", path), call. = FALSE)
    }
  }

  if (!file.rename(temp, path)) {
    if (had_target && file.exists(backup)) restore_backup()
    stop(sprintf("Could not publish replacement for '%s'.", path), call. = FALSE)
  }

  if (identical(Sys.getenv("WF16S_TEST_MODE"), "1") &&
      identical(Sys.getenv("WF16S_INJECT_ATOMIC_POST_PUBLISH_CORRUPTION"), basename(path))) {
    writeBin(charToRaw("injected-corruption"), path)
  }

  pub_info <- file.info(path)
  if (is.na(pub_info$size) || pub_info$size != temp_size) {
    if (had_target && file.exists(backup)) restore_backup()
    stop(sprintf("Published file '%s' failed size verification.", path), call. = FALSE)
  }
  pub_hash <- hash_file(path)
  if (!identical(pub_hash, temp_hash)) {
    if (had_target && file.exists(backup)) restore_backup()
    stop(sprintf("Published file '%s' failed hash verification.", path), call. = FALSE)
  }

  if (!is.null(orig_mode) && !is.na(orig_mode)) {
    tryCatch(Sys.chmod(path, mode = orig_mode), error = function(e) NULL)
  }

  replacement_verified <- TRUE

  if (had_target && file.exists(backup)) {
    if (unlink(backup, force = TRUE) != 0L || file.exists(backup)) {
      stop(sprintf("Verified replacement for '%s' could not remove backup '%s'.", path, backup),
           call. = FALSE)
    }
  }

  invisible(path)
}

atomic_write_lines <- function(text, path) atomic_replace(path, function(temp) writeLines(text, temp))
atomic_write_yaml <- function(value, path) atomic_replace(path, function(temp) yaml::write_yaml(value, temp))
atomic_write_json <- function(value, path, auto_unbox = TRUE) {
  atomic_replace(path, function(temp) jsonlite::write_json(value, temp, pretty = TRUE,
    auto_unbox = auto_unbox, null = "null"))
}

prepare_run_staging <- function(final_root) {
  parent <- dirname(final_root)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(pattern = paste0(".", basename(final_root), ".staging-"), tmpdir = parent)
  if (!dir.create(stage, recursive = FALSE, showWarnings = FALSE)) {
    stop(sprintf("Could not create run staging directory '%s'.", stage), call. = FALSE)
  }
  normalizePath(stage, winslash = "/", mustWork = TRUE)
}

snapshot_staging_directory <- function(stage) {
  parent <- dirname(stage)
  snapshot <- tempfile(pattern = paste0(".", basename(stage), ".snapshot-"), tmpdir = parent)
  if (!dir.create(snapshot, recursive = FALSE, showWarnings = FALSE)) {
    stop(sprintf("Could not create module staging snapshot '%s'.", snapshot), call. = FALSE)
  }
  entries <- list.files(stage, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                        full.names = FALSE)
  directories <- entries[dir.exists(file.path(stage, entries))]
  for (relative in directories) {
    dir.create(file.path(snapshot, relative), recursive = TRUE, showWarnings = FALSE)
  }
  files <- setdiff(entries, directories)
  for (relative in files) {
    target <- file.path(snapshot, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(stage, relative), target, overwrite = FALSE)) {
      unlink(snapshot, recursive = TRUE, force = TRUE)
      stop(sprintf("Could not snapshot staged path '%s'.", relative), call. = FALSE)
    }
  }
  normalizePath(snapshot, winslash = "/", mustWork = TRUE)
}

restore_staging_directory <- function(stage, snapshot) {
  if (!dir.exists(snapshot)) stop("Module staging snapshot is missing.", call. = FALSE)
  if (dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
  if (!dir.create(stage, recursive = FALSE, showWarnings = FALSE)) {
    stop(sprintf("Could not recreate module staging directory '%s'.", stage), call. = FALSE)
  }
  entries <- list.files(snapshot, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                        full.names = FALSE)
  directories <- entries[dir.exists(file.path(snapshot, entries))]
  for (relative in directories) {
    dir.create(file.path(stage, relative), recursive = TRUE, showWarnings = FALSE)
  }
  files <- setdiff(entries, directories)
  for (relative in files) {
    target <- file.path(stage, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(snapshot, relative), target, overwrite = FALSE)) {
      stop(sprintf("Could not restore staged path '%s'.", relative), call. = FALSE)
    }
  }
  invisible(stage)
}

validate_prior_output <- function(final_root, overwrite) {
  if (!dir.exists(final_root)) return(NULL)
  entries <- list.files(final_root, all.files = TRUE, no.. = TRUE)
  if (!length(entries)) return(NULL)
  manifest_path <- file.path(final_root, "run_manifest.json")
  if (!isTRUE(overwrite)) {
    stop(sprintf("Output directory '%s' is non-empty; pass --overwrite for an owned prior run.", final_root), call. = FALSE)
  }
  if (!file.exists(manifest_path)) {
    stop(sprintf("E_OUTPUT_UNOWNED: non-empty output directory lacks a valid prior manifest: '%s'", final_root), call. = FALSE)
  }
  prior <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE), error = function(e) NULL)
  root <- canonicalize_root_path(final_root)
  prior_root <- if (!is.null(prior$output_root)) {
    canonicalize_root_path(as.character(prior$output_root))
  } else {
    NA_character_
  }
  if (is.null(prior) || !identical(prior$pipeline, "ont-wf16s-postprocess") ||
      is.na(prior_root) || !paths_are_same(prior_root, root)) {
    stop(sprintf("E_OUTPUT_UNOWNED: prior manifest does not prove ownership of '%s'", final_root), call. = FALSE)
  }
  if (!prior$run_status %in% c("completed", "failed")) {
    stop(sprintf("E_OUTPUT_UNOWNED: prior output '%s' has an invalid run state.", final_root), call. = FALSE)
  }

  path_is_within <- function(path) {
    path_is_same_or_descendant(path, root)
  }
  relative_owned_path <- function(path) {
    candidate <- if (is_absolute_path(path)) path else file.path(root, path)
    if (!path_is_within(candidate)) {
      stop(sprintf("E_OUTPUT_UNOWNED: prior manifest declares path outside output root: '%s'", path),
           call. = FALSE)
    }
    normalized <- canonicalize_root_path(candidate)
    substring(normalized, nchar(root) + 2L)
  }
  prior_files <- list.files(final_root, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                            full.names = FALSE)
  prior_files <- prior_files[!dir.exists(file.path(final_root, prior_files))]
  current_contract <- identical(prior$schema_version, 2L) &&
    (identical(prior$schema_revision, 1L) || identical(prior$schema_revision, 2L) ||
     identical(prior$schema_revision, 3L)) &&
    !is.null(prior$owned_outputs) &&
    is.list(prior$environment) && !is.null(prior$environment$library_paths) &&
    !is.null(prior$environment$package_locations)
  if (current_contract) {
    tryCatch(validate_manifest_v2(prior, physical_root = final_root),
             error = function(e) stop(sprintf("E_OUTPUT_UNOWNED: invalid prior manifest: %s", e$message),
                                       call. = FALSE))
    owned <- as.character(unlist(prior$owned_outputs, use.names = FALSE))
    owned <- unique(vapply(owned, relative_owned_path, character(1)))
  } else {
    if (!identical(prior$run_status, "completed")) {
      stop("E_OUTPUT_UNOWNED: failed prior runs must use the current strict manifest contract.",
           call. = FALSE)
    }
    # v0.4.0 and earlier manifests had no ownership inventory. Derive one only
    # when every physical file is explicitly listed by a legacy module record.
    if (is.null(prior$modules) || !is.list(prior$modules)) {
      stop("E_OUTPUT_MIGRATION_REQUIRED: prior manifest predates ownership tracking and cannot be migrated safely.",
           call. = FALSE)
    }
    declared <- c(unlist(prior$owned_outputs %||% list(), use.names = FALSE),
                  unlist(lapply(prior$modules, function(record) record$outputs %||% character(0)),
                         use.names = FALSE))
    owned <- unique(c("run_manifest.json", "resolved_config.yml", "session_info.txt",
                      vapply(as.character(declared), relative_owned_path, character(1))))
    missing_declarations <- setdiff(prior_files, owned)
    if (length(missing_declarations)) {
      stop(sprintf("E_OUTPUT_MIGRATION_REQUIRED: legacy output contains undeclared file(s): %s",
                   paste(missing_declarations, collapse = ", ")), call. = FALSE)
    }
    prior$ownership_migrated <- TRUE
  }
  if (length(setdiff(owned, prior_files))) {
    stop(sprintf("E_OUTPUT_UNOWNED: prior manifest declares missing file(s): %s",
                 paste(setdiff(owned, prior_files), collapse = ", ")), call. = FALSE)
  }
  prior$owned_outputs <- if (exists("json_array", mode = "function")) json_array(owned) else owned
  prior
}

preserve_unowned_outputs <- function(final_root, stage, prior_manifest) {
  if (is.null(prior_manifest) || !dir.exists(final_root)) return(character(0))
  owned <- unique(c("run_manifest.json", "resolved_config.yml", "session_info.txt",
                    unlist(prior_manifest$owned_outputs %||% list(), use.names = FALSE)))
  files <- list.files(final_root, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  files <- files[!dir.exists(file.path(final_root, files))]
  unowned <- sort(setdiff(files, owned))
  for (relative in unowned) {
    target <- file.path(stage, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(final_root, relative), target, overwrite = FALSE)) {
      stop(sprintf("Could not preserve unowned output '%s'.", relative), call. = FALSE)
    }
  }
  unowned
}

prepare_module_staging <- function(stage, module_name) {
  parent <- dirname(stage)
  module_stage <- tempfile(pattern = paste0(".", basename(stage), ".", module_name, "-"), tmpdir = parent)
  if (!dir.create(module_stage, recursive = TRUE, showWarnings = FALSE)) {
    stop(sprintf("Could not create private module staging root '%s'.", module_stage), call. = FALSE)
  }
  normalizePath(module_stage, winslash = "/", mustWork = TRUE)
}

publish_module_staging <- function(module_stage, run_stage, declared_outputs) {
  norm_mod <- canonicalize_root_path(module_stage)
  norm_stage <- canonicalize_root_path(run_stage)
  move_plan <- lapply(declared_outputs, function(src_path) {
    norm_src <- canonicalize_root_path(src_path)
    if (!path_is_descendant(norm_src, norm_mod)) {
      stop(sprintf("Declared module output '%s' is not within module stage '%s'.", src_path, module_stage),
           call. = FALSE)
    }
    rel_path <- substring(norm_src, nchar(norm_mod) + 2L)
    dest_path <- file.path(norm_stage, rel_path)
    list(src = norm_src, relative = rel_path, destination = dest_path)
  })
  relative_paths <- vapply(move_plan, function(item) item$relative, character(1))
  if (anyDuplicated(tolower(relative_paths))) {
    stop("E_OUTPUT_COLLISION: module declared duplicate case-folded output paths.", call. = FALSE)
  }
  collisions <- relative_paths[vapply(move_plan, function(item) file.exists(item$destination) ||
                                        dir.exists(item$destination), logical(1))]
  if (length(collisions)) {
    stop(sprintf("E_OUTPUT_COLLISION: module output would overwrite an existing staged path: %s",
                 paste(collisions, collapse = ", ")), call. = FALSE)
  }
  for (item in move_plan) {
    norm_src <- item$src
    rel_path <- item$relative
    dest_path <- item$destination
    dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(norm_src, dest_path)) {
      if (!file.copy(norm_src, dest_path, overwrite = FALSE) || unlink(norm_src, force = TRUE) != 0L) {
        stop(sprintf("Could not move module output '%s' to run stage.", rel_path), call. = FALSE)
      }
    }
  }
  unlink(module_stage, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

cleanup_module_staging <- function(module_stage) {
  if (!is.null(module_stage) && dir.exists(module_stage)) {
    unlink(module_stage, recursive = TRUE, force = TRUE)
  }
  invisible(TRUE)
}

census_physical_files <- function(root_dir) {
  if (!dir.exists(root_dir)) return(character(0))
  all_paths <- list.files(root_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE, full.names = FALSE)
  is_file <- vapply(file.path(root_dir, all_paths), function(p) !dir.exists(p), logical(1))
  sort(all_paths[is_file])
}

verify_physical_file_census <- function(stage, owned_outputs, preserved_unowned_outputs) {
  physical <- census_physical_files(stage)
  expected_set <- sort(unique(c(setdiff(owned_outputs, "run_manifest.json"), preserved_unowned_outputs)))
  physical_set <- sort(unique(physical))

  if (anyDuplicated(tolower(physical_set))) {
    stop("Physical staged files contain case-folded collisions.", call. = FALSE)
  }
  if (anyDuplicated(tolower(expected_set))) {
    stop("Expected output census contains case-folded collisions.", call. = FALSE)
  }

  missing_from_physical <- setdiff(expected_set, physical_set)
  extra_in_physical <- setdiff(physical_set, expected_set)

  if (length(missing_from_physical) > 0L) {
    stop(sprintf("E_CENSUS_MISMATCH: expected outputs missing from physical stage: %s",
                 paste(missing_from_physical, collapse = ", ")), call. = FALSE)
  }
  if (length(extra_in_physical) > 0L) {
    stop(sprintf("E_CENSUS_MISMATCH: unowned physical files found in stage: %s",
                 paste(extra_in_physical, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

canonicalize_root_path <- function(path) {
  path <- gsub("\\\\", "/", path)
  parent <- dirname(path)
  if (dir.exists(parent)) {
    norm_parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
    file.path(norm_parent, basename(path))
  } else if (file.exists(path) || dir.exists(path)) {
    normalizePath(path, winslash = "/", mustWork = TRUE)
  } else {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }
}

path_identity_key <- function(path, os_type = .Platform$OS.type) {
  if (identical(os_type, "windows")) tolower(path) else path
}

paths_are_same <- function(left, right, os_type = .Platform$OS.type) {
  identical(path_identity_key(left, os_type), path_identity_key(right, os_type))
}

path_is_descendant <- function(path, root, os_type = .Platform$OS.type) {
  path_key <- path_identity_key(canonicalize_root_path(path), os_type)
  root_key <- path_identity_key(canonicalize_root_path(root), os_type)
  startsWith(path_key, paste0(root_key, "/"))
}

path_is_same_or_descendant <- function(path, root, os_type = .Platform$OS.type) {
  paths_are_same(canonicalize_root_path(path), canonicalize_root_path(root), os_type) ||
    path_is_descendant(path, root, os_type)
}

get_output_lock_path <- function(final_root) {
  canonical <- canonicalize_root_path(final_root)
  parent <- dirname(canonical)
  root_hash <- digest::digest(path_identity_key(canonical), algo = "sha256")
  file.path(parent, sprintf(".%s.wf16s_output.lock", root_hash))
}

.wf16s_active_output_locks <- new.env(parent = emptyenv())

acquire_output_lock <- function(final_root, timeout_ms = 10000) {
  lock_path <- get_output_lock_path(final_root)
  norm_lock_path <- path_identity_key(canonicalize_root_path(lock_path))
  if (exists(norm_lock_path, envir = .wf16s_active_output_locks, inherits = FALSE)) {
    stop(sprintf("E_OUTPUT_BUSY: output directory '%s' is locked by another process (lock '%s').",
                 final_root, lock_path), call. = FALSE)
  }
  dir.create(dirname(lock_path), recursive = TRUE, showWarnings = FALSE)
  lock_handle <- tryCatch(filelock::lock(lock_path, timeout = timeout_ms),
                          error = function(e) NULL)
  if (is.null(lock_handle)) {
    stop(sprintf("E_OUTPUT_BUSY: output directory '%s' is locked by another process (lock '%s').",
                 final_root, lock_path), call. = FALSE)
  }
  assign(norm_lock_path, lock_handle, envir = .wf16s_active_output_locks)
  lock_handle
}

release_output_lock <- function(lock_handle) {
  if (!is.null(lock_handle)) {
    for (name in ls(.wf16s_active_output_locks)) {
      if (identical(.wf16s_active_output_locks[[name]], lock_handle)) {
        rm(list = name, envir = .wf16s_active_output_locks)
        break
      }
    }
    tryCatch(filelock::unlock(lock_handle), error = function(e) NULL)
  }
  invisible(TRUE)
}

.wf16s_active_taxonomy_locks <- new.env(parent = emptyenv())

taxonomy_cache_identity <- function(cache_path) {
  if (!is.character(cache_path) || length(cache_path) != 1L || is.na(cache_path) ||
      !nzchar(trimws(cache_path)) || !file.exists(cache_path) || dir.exists(cache_path)) {
    stop(sprintf("Taxonomy cache must already exist as a regular file: '%s'.", cache_path),
         call. = FALSE)
  }
  normalizePath(cache_path, winslash = "/", mustWork = TRUE)
}

get_taxonomy_lock_path <- function(cache_path) {
  paste0(taxonomy_cache_identity(cache_path), ".lock")
}

acquire_taxonomy_lock <- function(cache_path, timeout_ms = 10000) {
  lock_path <- get_taxonomy_lock_path(cache_path)
  norm_lock_path <- path_identity_key(canonicalize_root_path(lock_path))
  if (exists(norm_lock_path, envir = .wf16s_active_taxonomy_locks, inherits = FALSE)) {
    stop(sprintf("E_TAXONOMY_CACHE_BUSY: taxonomy cache '%s' is already locked in this process.",
                 cache_path), call. = FALSE)
  }
  lock_handle <- tryCatch(filelock::lock(lock_path, timeout = timeout_ms), error = function(e) NULL)
  if (is.null(lock_handle)) {
    stop(sprintf("E_TAXONOMY_CACHE_BUSY: taxonomy cache '%s' is locked by another process.",
                 cache_path), call. = FALSE)
  }
  assign(norm_lock_path, lock_handle, envir = .wf16s_active_taxonomy_locks)
  lock_handle
}

release_taxonomy_lock <- function(lock_handle) {
  if (!is.null(lock_handle)) {
    for (name in ls(.wf16s_active_taxonomy_locks)) {
      if (identical(.wf16s_active_taxonomy_locks[[name]], lock_handle)) {
        rm(list = name, envir = .wf16s_active_taxonomy_locks)
        break
      }
    }
    tryCatch(filelock::unlock(lock_handle), error = function(e) NULL)
  }
  invisible(TRUE)
}

get_taxonomy_journal_path <- function(cache_path) {
  paste0(taxonomy_cache_identity(cache_path),
         ".wf16s_transaction.json")
}

new_transaction_id <- function() {
  entropy <- paste(
    utc_timestamp_now(), Sys.getpid(), tempfile("wf16s_txn_"),
    sep = "|"
  )
  paste0("tx-", digest::digest(entropy, algo = "sha256", serialize = FALSE))
}

valid_transaction_id <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^tx-[0-9a-f]{64}$", value)
}

assert_taxonomy_backup_path <- function(path, cache_identity) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!paths_are_same(dirname(normalized), dirname(cache_identity)) ||
      !startsWith(basename(normalized), ".wf16s_tax_backup_")) {
    stop(sprintf("Invalid taxonomy recovery backup path '%s'.", path), call. = FALSE)
  }
  invisible(normalized)
}

write_taxonomy_journal <- function(cache_path, backup, output_root, original_sha256,
                                   candidate_sha256, phase, transaction_id = NULL) {
  if (!valid_transaction_id(transaction_id)) {
    stop("Taxonomy journal requires a valid transaction_id.", call. = FALSE)
  }
  if (!is.character(phase) || length(phase) != 1L ||
      !phase %in% c("prepared", "candidate_committed", "output_published")) {
    stop("Invalid taxonomy journal phase.", call. = FALSE)
  }
  cache_identity <- taxonomy_cache_identity(cache_path)
  backup_path <- assert_taxonomy_backup_path(backup, cache_identity)
  payload <- list(
    journal_schema_version = 2L,
    cache_path = cache_identity,
    configured_cache_path = gsub("\\\\", "/", as.character(cache_path)),
    backup = backup_path,
    output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE),
    original_sha256 = original_sha256,
    candidate_sha256 = candidate_sha256,
    transaction_id = as.character(transaction_id),
    phase = phase,
    timestamp = utc_timestamp_now(),
    pid = Sys.getpid()
  )
  atomic_write_json(payload, get_taxonomy_journal_path(cache_path))
  invisible(payload)
}

remove_file_checked <- function(path, context) {
  if (!file.exists(path)) return(invisible(TRUE))
  status <- unlink(path, force = TRUE)
  if (status != 0L || file.exists(path)) {
    stop(sprintf("E_OUTPUT_CLEANUP_FAILED: could not remove %s '%s'.", context, path),
         call. = FALSE)
  }
  invisible(TRUE)
}

cleanup_taxonomy_journal <- function(cache_path, backup = NULL) {
  # Delete the expendable backup first. If that cleanup fails, retain the
  # journal so operators still have an authoritative recovery record.
  if (!is.null(backup)) remove_file_checked(backup, "taxonomy-cache backup")
  remove_file_checked(get_taxonomy_journal_path(cache_path), "taxonomy transaction journal")
  invisible(TRUE)
}

recover_taxonomy_journal <- function(cache_path) {
  journal_path <- get_taxonomy_journal_path(cache_path)
  if (!file.exists(journal_path)) return(invisible(FALSE))
  journal <- tryCatch(jsonlite::fromJSON(journal_path, simplifyVector = FALSE),
                      error = function(e) NULL)
  required <- c("cache_path", "backup", "output_root", "original_sha256",
                "candidate_sha256", "phase")
  if (is.null(journal) || length(setdiff(required, names(journal))) ||
      !journal$phase %in% c("prepared", "candidate_committed", "output_published")) {
    stop(sprintf("E_TAXONOMY_RECOVERY_REQUIRED: corrupt taxonomy journal at '%s'.",
                 journal_path), call. = FALSE)
  }
  legacy <- is.null(journal$journal_schema_version)
  if (!legacy && !identical(journal$journal_schema_version, 2L)) {
    stop(sprintf("E_TAXONOMY_RECOVERY_REQUIRED: unsupported taxonomy journal schema at '%s'.",
                 journal_path), call. = FALSE)
  }
  if (!legacy && !valid_transaction_id(journal$transaction_id)) {
    stop(sprintf("E_TAXONOMY_RECOVERY_REQUIRED: taxonomy journal lacks a valid transaction_id at '%s'.",
                 journal_path), call. = FALSE)
  }
  canonical_cache <- taxonomy_cache_identity(cache_path)
  recorded_cache <- normalizePath(as.character(journal$cache_path), winslash = "/",
                                  mustWork = FALSE)
  backup <- normalizePath(as.character(journal$backup), winslash = "/", mustWork = FALSE)
  if (!paths_are_same(canonical_cache, recorded_cache) ||
      !paths_are_same(dirname(canonical_cache), dirname(backup)) ||
      !startsWith(basename(backup), ".wf16s_tax_backup_") ||
      !grepl("^[0-9a-f]{64}$", journal$original_sha256) ||
      !grepl("^[0-9a-f]{64}$", journal$candidate_sha256 %||% "")) {
    stop(sprintf("E_TAXONOMY_RECOVERY_REQUIRED: invalid taxonomy journal at '%s'.",
                 journal_path), call. = FALSE)
  }
  if (!file.exists(cache_path) || !file.exists(backup) ||
      !identical(compute_file_hash(backup), journal$original_sha256)) {
    stop(sprintf("E_TAXONOMY_RECOVERY_REQUIRED: taxonomy recovery material is missing or invalid at '%s'.",
                 journal_path), call. = FALSE)
  }
  current_sha256 <- compute_file_hash(cache_path)
  final_valid <- manifest_is_valid_run(journal$output_root)
  final_manifest <- if (final_valid) tryCatch(jsonlite::fromJSON(
    file.path(journal$output_root, "run_manifest.json"), simplifyVector = FALSE
  ), error = function(e) NULL) else NULL
  if (legacy && final_valid) {
    stop(sprintf(
      "E_TAXONOMY_RECOVERY_REQUIRED: legacy taxonomy journal has ambiguous completed-output identity at '%s'.",
      journal_path), call. = FALSE)
  }
  transaction_matches <- legacy ||
    (!is.null(final_manifest) && identical(final_manifest$transaction_id,
                                           journal$transaction_id))
  if (final_valid && !transaction_matches) {
    stop(sprintf(
      "E_TAXONOMY_RECOVERY_REQUIRED: completed output transaction identity does not match journal '%s'.",
      journal_path), call. = FALSE)
  }
  candidate_file <- file.path(journal$output_root, "07_Kreport",
                              "resolved_taxonomy_cache.json")
  final_has_candidate <- final_valid && file.exists(candidate_file) &&
    identical(compute_file_hash(candidate_file), journal$candidate_sha256)
  if (identical(current_sha256, journal$original_sha256)) {
    if (final_valid && !final_has_candidate) {
      stop(sprintf(
        "E_TAXONOMY_RECOVERY_REQUIRED: completed output lacks the recorded taxonomy candidate for '%s'.",
        cache_path), call. = FALSE)
    }
    if (final_has_candidate) {
      atomic_replace(cache_path, function(temp) {
        if (!file.copy(candidate_file, temp, overwrite = TRUE)) {
          stop(sprintf("Could not restore committed taxonomy candidate from '%s'.",
                       candidate_file), call. = FALSE)
        }
      }, expected_sha256 = journal$candidate_sha256)
    }
    cleanup_taxonomy_journal(cache_path, backup)
    return(invisible(TRUE))
  }
  if (!identical(current_sha256, journal$candidate_sha256)) {
    stop(sprintf("E_TAXONOMY_CACHE_CHANGED: taxonomy cache '%s' changed outside the recorded transaction; journal retained.",
                 cache_path), call. = FALSE)
  }
  if (final_valid && !final_has_candidate) {
    stop(sprintf(
      "E_TAXONOMY_RECOVERY_REQUIRED: completed output lacks the recorded taxonomy candidate for '%s'.",
      cache_path), call. = FALSE)
  }
  if (final_has_candidate) {
    cleanup_taxonomy_journal(cache_path, backup)
    return(invisible(TRUE))
  }
  atomic_replace(cache_path, function(temp) {
    if (!file.copy(backup, temp, overwrite = TRUE)) {
      stop(sprintf("Could not copy taxonomy recovery backup '%s'.", backup), call. = FALSE)
    }
  }, expected_sha256 = journal$original_sha256)
  cleanup_taxonomy_journal(cache_path, backup)
  invisible(TRUE)
}

get_output_journal_path <- function(final_root) {
  canonical <- canonicalize_root_path(final_root)
  parent <- dirname(canonical)
  root_hash <- digest::digest(path_identity_key(canonical), algo = "sha256")
  file.path(parent, sprintf(".%s.wf16s_journal.json", root_hash))
}

assert_publication_temp_path <- function(path, final_root, prefix, nullable = TRUE) {
  if (is.null(path)) {
    if (!nullable) stop("Publication journal path is required.", call. = FALSE)
    return(invisible(NULL))
  }
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  canonical <- canonicalize_root_path(final_root)
  expected_prefix <- paste0(".", basename(canonical), prefix)
  if (!paths_are_same(dirname(normalized), dirname(canonical)) ||
      !startsWith(basename(normalized), expected_prefix)) {
    stop(sprintf("Invalid publication journal temporary path '%s'.", path), call. = FALSE)
  }
  invisible(normalized)
}

write_publication_journal <- function(final_root, stage = NULL, backup = NULL,
                                      phase = "prepared", transaction_id = NULL) {
  if (!valid_transaction_id(transaction_id)) {
    stop("Publication journal requires a valid transaction_id.", call. = FALSE)
  }
  if (!is.character(phase) || length(phase) != 1L ||
      !phase %in% c("prepared", "prior_moved", "stage_published")) {
    stop("Invalid publication journal phase.", call. = FALSE)
  }
  canonical_final <- canonicalize_root_path(final_root)
  stage_path <- assert_publication_temp_path(stage, final_root, ".staging-")
  backup_path <- assert_publication_temp_path(backup, final_root, ".previous-")
  journal_path <- get_output_journal_path(final_root)
  payload <- list(
    journal_schema_version = 2L,
    final_root = canonical_final,
    stage = stage_path,
    backup = backup_path,
    had_prior = !is.null(backup_path),
    transaction_id = as.character(transaction_id),
    phase = phase,
    timestamp = utc_timestamp_now(),
    pid = Sys.getpid()
  )
  atomic_write_json(payload, journal_path)
  invisible(journal_path)
}

remove_publication_journal <- function(final_root) {
  journal_path <- get_output_journal_path(final_root)
  if (file.exists(journal_path)) {
    status <- unlink(journal_path, force = TRUE)
    if (status != 0L || file.exists(journal_path)) {
      stop(sprintf(
        "E_OUTPUT_CLEANUP_FAILED: could not remove publication journal '%s'.",
        journal_path
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

manifest_is_valid_run <- function(dir_path) {
  if (!dir.exists(dir_path)) return(FALSE)
  manifest_file <- file.path(dir_path, "run_manifest.json")
  if (!file.exists(manifest_file)) return(FALSE)
  manifest <- tryCatch(jsonlite::fromJSON(manifest_file, simplifyVector = FALSE),
                       error = function(e) NULL)
  valid_identity <- !is.null(manifest) &&
    identical(manifest$pipeline, "ont-wf16s-postprocess") &&
    identical(manifest$run_status, "completed") &&
    identical(manifest$schema_version, 2L) &&
    length(manifest$schema_revision) == 1L &&
    isTRUE(manifest$schema_revision %in% c(1L, 2L, 3L))
  if (!isTRUE(valid_identity)) return(FALSE)
  isTRUE(tryCatch({
    validate_manifest_v2(manifest, physical_root = dir_path)
    TRUE
  }, error = function(e) FALSE))
}

remove_tree_checked <- function(path, context) {
  if (!dir.exists(path)) return(invisible(TRUE))
  status <- unlink(path, recursive = TRUE, force = TRUE)
  if (status != 0L || dir.exists(path)) {
    stop(sprintf("E_OUTPUT_CLEANUP_FAILED: could not remove %s '%s'; publication journal retained.",
                 context, path), call. = FALSE)
  }
  invisible(TRUE)
}

recover_publication_journal <- function(final_root) {
  journal_path <- get_output_journal_path(final_root)
  if (!file.exists(journal_path)) return(invisible(FALSE))

  journal <- tryCatch(jsonlite::fromJSON(journal_path, simplifyVector = FALSE),
                      error = function(e) NULL)
  if (is.null(journal) || is.null(journal$final_root)) {
    stop(sprintf("E_OUTPUT_RECOVERY_REQUIRED: corrupt publication journal at '%s'.", journal_path),
         call. = FALSE)
  }

  canonical_final <- canonicalize_root_path(final_root)
  journal_final <- canonicalize_root_path(as.character(journal$final_root))
  if (!paths_are_same(canonical_final, journal_final)) {
    return(invisible(FALSE))
  }

  legacy <- is.null(journal$journal_schema_version)
  if (!legacy) {
    if (!identical(journal$journal_schema_version, 2L) ||
        !valid_transaction_id(journal$transaction_id) ||
        !is.logical(journal$had_prior) || length(journal$had_prior) != 1L ||
        is.na(journal$had_prior) ||
        !journal$phase %in% c("prepared", "prior_moved", "stage_published")) {
      stop(sprintf("E_OUTPUT_RECOVERY_REQUIRED: invalid publication journal schema at '%s'.",
                   journal_path), call. = FALSE)
    }
    assert_publication_temp_path(journal$stage, final_root, ".staging-")
    assert_publication_temp_path(journal$backup, final_root, ".previous-",
                                 nullable = !isTRUE(journal$had_prior))
    if (!isTRUE(journal$had_prior) && !is.null(journal$backup)) {
      stop(sprintf("E_OUTPUT_RECOVERY_REQUIRED: journal backup disagrees with had_prior at '%s'.",
                   journal_path), call. = FALSE)
    }
  }

  backup <- journal$backup
  has_final <- dir.exists(final_root)
  final_valid <- has_final && manifest_is_valid_run(final_root)
  if (final_valid && !legacy) {
    final_manifest <- tryCatch(jsonlite::fromJSON(
      file.path(final_root, "run_manifest.json"), simplifyVector = FALSE
    ), error = function(e) NULL)
    if (is.null(final_manifest) ||
        !identical(final_manifest$transaction_id, journal$transaction_id)) {
      final_valid <- FALSE
    }
  }
  if (legacy && final_valid) {
    stop(sprintf("E_OUTPUT_RECOVERY_REQUIRED: legacy publication journal has ambiguous completed-output identity at '%s'.",
                 journal_path), call. = FALSE)
  }
  has_backup <- !is.null(backup) && dir.exists(backup)
  backup_valid <- has_backup && manifest_is_valid_run(backup)

  if (!has_final && backup_valid) {
    if (!file.rename(backup, final_root)) {
      stop(sprintf("E_OUTPUT_RECOVERY_REQUIRED: failed to restore backup '%s' to '%s'.",
                   backup, final_root), call. = FALSE)
    }
    remove_publication_journal(final_root)
    cat(sprintf("[RECOVERY] Restored prior completed run from crash backup '%s'.\n", backup),
        file = stderr())
    return(invisible(TRUE))
  } else if (final_valid) {
    if (has_backup) remove_tree_checked(backup, "crash backup")
    remove_publication_journal(final_root)
    cat(sprintf("[RECOVERY] Resolved prior publication state for '%s'.\n", final_root),
        file = stderr())
    return(invisible(TRUE))
  } else {
    stop(sprintf("E_OUTPUT_RECOVERY_REQUIRED: incomplete publication transaction at '%s'; manual recovery required.",
                 final_root), call. = FALSE)
  }
}

publish_staged_run <- function(stage, final_root, transaction_id = NULL) {
  parent <- dirname(final_root)
  backup <- tempfile(pattern = paste0(".", basename(final_root), ".previous-"), tmpdir = parent)
  had_prior <- dir.exists(final_root)

  write_publication_journal(final_root, stage = stage, backup = if (had_prior) backup else NULL,
                            phase = "prepared", transaction_id = transaction_id)

  committed <- FALSE
  on.exit({
    if (!committed && had_prior && dir.exists(backup) && !dir.exists(final_root)) {
      if (!file.rename(backup, final_root)) {
        warning(sprintf(
          "E_OUTPUT_RECOVERY_REQUIRED: failed to restore backup '%s' to '%s'; publication journal retained.",
          backup, final_root
        ), call. = FALSE)
      } else {
        tryCatch(remove_publication_journal(final_root),
                 error = function(e) warning(conditionMessage(e), call. = FALSE))
      }
    }
  }, add = TRUE)

  if (had_prior) {
    if (!file.rename(final_root, backup)) {
      stop("Could not preserve the prior completed run.", call. = FALSE)
    }
    write_publication_journal(final_root, stage = stage, backup = backup, phase = "prior_moved",
                              transaction_id = transaction_id)
  }

  if (!file.rename(stage, final_root)) {
    stop("Could not publish staged run.", call. = FALSE)
  }
  committed <- TRUE
  write_publication_journal(final_root, stage = stage, backup = if (had_prior) backup else NULL,
                            phase = "stage_published", transaction_id = transaction_id)

  if (had_prior && dir.exists(backup)) {
    remove_tree_checked(backup, "prior-run backup")
  }
  remove_publication_journal(final_root)
  invisible(final_root)
}
