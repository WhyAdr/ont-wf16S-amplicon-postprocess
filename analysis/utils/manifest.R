# =============================================================================
# Manifest schema v2, revision 1
# =============================================================================

json_array <- function(x) {
  if (is.null(x) || length(x) == 0L) return(list())
  unname(lapply(as.list(x), function(value) if (is.list(value)) value else jsonlite::unbox(value)))
}

is_json_array <- function(x) is.list(x) && is.null(names(x))

manifest_fail <- function(path, message) {
  stop(sprintf("Manifest schema v2 violation at '%s': %s", path, message), call. = FALSE)
}

assert_manifest_array <- function(value, path, element_type = NULL) {
  if (!is_json_array(value)) manifest_fail(path, "expected a JSON array.")
  if (!is.null(element_type) && length(value)) {
    valid <- vapply(value, function(item) {
      is.atomic(item) && length(item) == 1L && !is.na(item) && typeof(item) == element_type
    }, logical(1))
    if (!all(valid)) manifest_fail(path, sprintf("expected %s array elements.", element_type))
  }
  invisible(TRUE)
}

assert_manifest_scalar <- function(value, path, type, nullable = FALSE) {
  if (is.null(value) && nullable) return(invisible(TRUE))
  if (is.null(value) || length(value) != 1L || is.list(value) || is.na(value) || typeof(value) != type) {
    manifest_fail(path, sprintf("expected one %s%s.", type, if (nullable) " or null" else ""))
  }
  invisible(TRUE)
}

new_not_run_module_record <- function(reason = NULL) {
  list(status = "not_run", outputs = character(0), warnings = character(0), error = NULL,
       reason = reason, start_time = NULL, end_time = NULL, duration_seconds = NULL)
}

manifest_module_record <- function(record) {
  list(status = record$status, outputs = json_array(record$outputs %||% character(0)),
       warnings = json_array(record$warnings %||% character(0)), error = record$error %||% NULL,
       reason = record$reason %||% NULL, start_time = record$start_time %||% NULL,
       end_time = record$end_time %||% NULL, duration_seconds = record$duration_seconds %||% NULL)
}

is_utc_timestamp <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}([.][0-9]+)?Z$", x)
}

array_values <- function(x) unname(vapply(x, function(item) as.character(item), character(1)))

validate_manifest_fingerprint <- function(value, path) {
  if (!is.list(value) || is.null(names(value))) manifest_fail(path, "expected a fingerprint object.")
  for (field in c("path", "size_bytes", "mtime_utc", "sha256")) {
    if (is.null(value[[field]])) manifest_fail(paste0(path, ".", field), "field is required.")
  }
  assert_manifest_scalar(value$path, paste0(path, ".path"), "character")
  if (!is.numeric(value$size_bytes) || length(value$size_bytes) != 1L ||
      is.na(value$size_bytes) || is.list(value$size_bytes)) {
    manifest_fail(paste0(path, ".size_bytes"), "expected one numeric value.")
  }
  if (!is.finite(value$size_bytes) || value$size_bytes < 0 || value$size_bytes != floor(value$size_bytes)) {
    manifest_fail(paste0(path, ".size_bytes"), "expected a finite nonnegative integer.")
  }
  if (!is_utc_timestamp(value$mtime_utc)) manifest_fail(paste0(path, ".mtime_utc"), "expected UTC timestamp.")
  assert_manifest_scalar(value$sha256, paste0(path, ".sha256"), "character")
  if (!grepl("^[0-9a-f]{64}$", value$sha256)) manifest_fail(paste0(path, ".sha256"), "expected lowercase SHA-256.")
  invisible(TRUE)
}

validate_manifest_input_collection <- function(value, path) {
  assert_manifest_array(value, path)
  ids <- character(0)
  for (index in seq_along(value)) {
    record <- value[[index]]
    if (!is.list(record) || is.null(record$sample_id)) manifest_fail(path, "sample_id is required.")
    assert_manifest_scalar(record$sample_id, paste0(path, "[", index, "].sample_id"), "character")
    ids <- c(ids, record$sample_id)
    validate_manifest_fingerprint(record, paste0(path, "[", index, "]"))
  }
  if (anyDuplicated(ids)) manifest_fail(path, "sample IDs must be unique.")
  invisible(TRUE)
}

validate_manifest_v2 <- function(manifest, physical_root = NULL) {
  if (!is.list(manifest) || is.null(names(manifest))) manifest_fail("<root>", "expected an object.")
  required <- c("pipeline", "pipeline_version", "schema_version", "schema_revision",
    "config_schema_version", "run_status", "start_time", "end_time", "duration_seconds",
    "samples", "command", "cli", "inputs", "modules", "owned_outputs", "warnings",
    "package_versions", "environment", "output_root")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) manifest_fail("<root>", paste("missing required key(s):", paste(missing, collapse = ", ")))
  assert_manifest_scalar(manifest$pipeline, "pipeline", "character")
  assert_manifest_scalar(manifest$pipeline_version, "pipeline_version", "character")
  if (!grepl("^[0-9]+[.][0-9]+[.][0-9]+$", manifest$pipeline_version)) {
    manifest_fail("pipeline_version", "expected SemVer.")
  }
  if (!identical(manifest$schema_version, 2L)) manifest_fail("schema_version", "expected integer 2.")
  if (!identical(manifest$schema_revision, 1L)) manifest_fail("schema_revision", "expected integer 1.")
  if (!identical(manifest$config_schema_version, 1L)) manifest_fail("config_schema_version", "expected integer 1.")
  if (!manifest$run_status %in% c("completed", "failed")) manifest_fail("run_status", "invalid value.")
  if (!is_utc_timestamp(manifest$start_time) || !is_utc_timestamp(manifest$end_time)) {
    manifest_fail("start_time/end_time", "expected ISO-8601 UTC timestamps.")
  }
  if (!is.numeric(manifest$duration_seconds) || length(manifest$duration_seconds) != 1L ||
      !is.finite(manifest$duration_seconds) || manifest$duration_seconds < 0) {
    manifest_fail("duration_seconds", "expected one finite nonnegative number.")
  }
  assert_manifest_scalar(manifest$project_name, "project_name", "character")
  if (!is.character(manifest$mode) || length(manifest$mode) != 1L ||
      !manifest$mode %in% c("single", "cohort")) manifest_fail("mode", "invalid value.")
  if (!is.numeric(manifest$seed) || length(manifest$seed) != 1L || !is.finite(manifest$seed) ||
      manifest$seed != floor(manifest$seed) || manifest$seed < 0) manifest_fail("seed", "invalid integer.")
  if (!is.null(manifest$config_file)) assert_manifest_scalar(manifest$config_file, "config_file", "character", nullable = TRUE)
  if (!is.null(manifest$git_commit)) assert_manifest_scalar(manifest$git_commit, "git_commit", "character", nullable = TRUE)
  if (!is.null(manifest$git_commit) && !grepl("^[0-9a-f]{40}$", manifest$git_commit)) manifest_fail("git_commit", "expected SHA-1 or null.")
  if (!is.null(manifest$git_dirty) && (!is.logical(manifest$git_dirty) || length(manifest$git_dirty) != 1L || is.na(manifest$git_dirty))) {
    manifest_fail("git_dirty", "expected logical or null.")
  }
  assert_manifest_scalar(manifest$source_digest_sha256, "source_digest_sha256", "character")
  if (!grepl("^[0-9a-f]{64}$", manifest$source_digest_sha256)) manifest_fail("source_digest_sha256", "expected SHA-256.")
  for (field in c("samples", "owned_outputs", "warnings")) {
    assert_manifest_array(manifest[[field]], field, "character")
    values <- array_values(manifest[[field]])
    if (anyDuplicated(values)) manifest_fail(field, "array elements must be unique.")
  }
  assert_manifest_array(manifest$command, "command", "character")
  assert_manifest_array(manifest$package_versions, "package_versions")
  packages <- vapply(manifest$package_versions, function(x) as.character(x$package %||% NA), character(1))
  if (anyNA(packages) || anyDuplicated(packages)) manifest_fail("package_versions", "package names must be present and unique.")
  for (index in seq_along(manifest$package_versions)) {
    record <- manifest$package_versions[[index]]
    assert_manifest_scalar(record$package, paste0("package_versions[", index, "].package"), "character")
    assert_manifest_scalar(record$version, paste0("package_versions[", index, "].version"), "character")
  }
  if (!is.list(manifest$cli) || is.null(names(manifest$cli))) manifest_fail("cli", "expected an object.")
  assert_manifest_array(manifest$cli$modules, "cli.modules", "character")
  requested <- array_values(manifest$cli$modules)
  if (anyDuplicated(requested)) manifest_fail("cli.modules", "module names must be unique.")
  for (field in c("validate_only", "keep_going", "overwrite", "allow_unlocked", "allow_dirty",
                  "online_preflight", "refresh_taxonomy", "krona")) {
    if (!is.logical(manifest$cli[[field]]) || length(manifest$cli[[field]]) != 1L || is.na(manifest$cli[[field]])) {
      manifest_fail(paste0("cli.", field), "expected one logical value.")
    }
  }
  if (!is.list(manifest$inputs) || is.null(names(manifest$inputs))) manifest_fail("inputs", "expected an object.")
  for (field in c("abundance_table", "params_json", "metadata", "taxonomy_cache")) {
    if (!is.null(manifest$inputs[[field]])) validate_manifest_fingerprint(manifest$inputs[[field]], paste0("inputs.", field))
  }
  validate_manifest_input_collection(manifest$inputs$assignments, "inputs.assignments")
  validate_manifest_input_collection(manifest$inputs$bamstats, "inputs.bamstats")
  if (!is.list(manifest$modules) || is.null(names(manifest$modules)) || anyDuplicated(names(manifest$modules))) {
    manifest_fail("modules", "expected a uniquely named object.")
  }
  registry <- c("qc", "alpha", "beta", "composition", "ordination", "shared", "kreport", "faprotax")
  if (!identical(sort(names(manifest$modules)), sort(registry))) {
    manifest_fail("modules", "module keys must exactly match the maintained registry.")
  }
  if (!setequal(requested, names(manifest$modules)[names(manifest$modules) %in% requested])) {
    manifest_fail("modules", "requested module keys disagree with cli.modules.")
  }
  all_outputs <- character(0)
  for (module_name in names(manifest$modules)) {
    path <- paste0("modules.", module_name)
    record <- manifest$modules[[module_name]]
    fields <- c("status", "outputs", "warnings", "error", "reason", "start_time", "end_time", "duration_seconds")
    if (!is.list(record) || length(setdiff(fields, names(record)))) manifest_fail(path, "invalid result object.")
    if (!is.character(record$status) || length(record$status) != 1L ||
        !record$status %in% c("completed", "skipped", "failed", "not_run")) manifest_fail(paste0(path, ".status"), "invalid status.")
    if (!(module_name %in% requested) && record$status != "not_run") {
      manifest_fail(path, "unrequested modules must be not_run.")
    }
    if (module_name %in% requested && record$status == "not_run" && manifest$run_status != "failed") {
      manifest_fail(path, "requested modules cannot be not_run.")
    }
    assert_manifest_array(record$outputs, paste0(path, ".outputs"), "character")
    assert_manifest_array(record$warnings, paste0(path, ".warnings"), "character")
    outputs <- array_values(record$outputs)
    if (anyDuplicated(outputs)) manifest_fail(paste0(path, ".outputs"), "outputs must be unique.")
    all_outputs <- c(all_outputs, outputs)
    if (record$status == "completed" && (!is.null(record$error) || !is.null(record$reason))) manifest_fail(path, "completed result cannot have error/reason.")
    if (record$status == "failed" && (is.null(record$error) || !nzchar(record$error))) manifest_fail(path, "failed result requires error.")
    if (record$status %in% c("skipped", "not_run") && (is.null(record$reason) || !nzchar(record$reason))) manifest_fail(path, "skipped/not_run result requires reason.")
    if (record$status == "not_run") {
      if (length(outputs) || length(record$warnings) || !is.null(record$error) || !is.null(record$start_time) ||
          !is.null(record$end_time) || !is.null(record$duration_seconds)) manifest_fail(path, "invalid not_run state.")
    } else {
      if (!is_utc_timestamp(record$start_time) || !is_utc_timestamp(record$end_time) ||
          !is.numeric(record$duration_seconds) || length(record$duration_seconds) != 1L ||
          !is.finite(record$duration_seconds) || record$duration_seconds < 0) manifest_fail(path, "invalid execution timing.")
    }
  }
  if (anyDuplicated(all_outputs)) manifest_fail("modules.*.outputs", "an output is declared more than once.")
  statuses <- vapply(manifest$modules[requested], `[[`, character(1), "status")
  expected_run <- if (any(statuses == "failed")) "failed" else "completed"
  if (!identical(manifest$run_status, expected_run)) manifest_fail("run_status", "disagrees with module states.")
  root <- normalizePath(manifest$output_root, winslash = "/", mustWork = FALSE)
  check_root <- normalizePath(physical_root %||% root, winslash = "/", mustWork = FALSE)
  if (!is.character(manifest$output_root) || length(manifest$output_root) != 1L || is.na(manifest$output_root)) {
    manifest_fail("output_root", "expected one character path.")
  }
  owned <- array_values(manifest$owned_outputs)
  if (any(!nzchar(owned)) || any(grepl("(^|/)[.][.](/|$)", owned)) ||
      any(grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]", owned))) {
    manifest_fail("owned_outputs", "paths must be non-empty relative paths beneath output_root.")
  }
  if (anyDuplicated(owned)) manifest_fail("owned_outputs", "paths must be unique.")
  if (!all(c("resolved_config.yml", "session_info.txt") %in% owned)) {
    manifest_fail("owned_outputs", "resolved_config.yml and session_info.txt must be owned.")
  }
  if (!file.exists(file.path(check_root, "resolved_config.yml")) ||
      !file.exists(file.path(check_root, "session_info.txt"))) {
    manifest_fail("owned_outputs", "required pipeline metadata files are missing.")
  }
  owned_normalized <- tolower(gsub("\\\\", "/", owned))
  for (output in all_outputs) {
    normalized <- normalizePath(output, winslash = "/", mustWork = FALSE)
    if (!startsWith(tolower(normalized), paste0(tolower(root), "/"))) manifest_fail("modules.*.outputs", "output is outside output_root.")
    relative <- substring(normalized, nchar(root) + 2L)
    physical <- file.path(check_root, relative)
    if (!file.exists(physical) || dir.exists(physical)) manifest_fail("modules.*.outputs", paste("missing regular file", output))
    if (!(tolower(gsub("\\\\", "/", relative)) %in% owned_normalized)) {
      manifest_fail("owned_outputs", sprintf("module output is not declared as owned: %s", relative))
    }
  }
  for (relative in owned) {
    physical <- file.path(check_root, relative)
    if (identical(gsub("\\\\", "/", relative), "run_manifest.json") && !file.exists(physical)) next
    if (!file.exists(physical) || dir.exists(physical)) manifest_fail("owned_outputs", paste("missing regular file", relative))
  }
  env <- manifest$environment
  if (!is.list(env) || !env$lock_status %in% c("synchronized", "mismatch", "missing", "not_checked") ||
      !is.logical(env$locked) || length(env$locked) != 1L || is.na(env$locked) ||
      !identical(env$locked, identical(env$lock_status, "synchronized"))) manifest_fail("environment", "inconsistent lock status.")
  assert_manifest_array(env$r_discrepancies, "environment.r_discrepancies", "character")
  assert_manifest_array(env$package_discrepancies, "environment.package_discrepancies", "character")
  assert_manifest_array(env$library_discrepancies, "environment.library_discrepancies", "character")
  assert_manifest_array(env$library_paths, "environment.library_paths", "character")
  if (!is.null(env$project_library)) assert_manifest_scalar(env$project_library, "environment.project_library", "character", nullable = TRUE)
  assert_manifest_array(env$package_locations, "environment.package_locations")
  package_location_names <- vapply(env$package_locations,
    function(record) as.character(record$package %||% NA), character(1))
  if (anyNA(package_location_names) || anyDuplicated(package_location_names)) {
    manifest_fail("environment.package_locations", "package names must be present and unique.")
  }
  for (index in seq_along(env$package_locations)) {
    record <- env$package_locations[[index]]
    assert_manifest_scalar(record$package, paste0("environment.package_locations[", index, "].package"), "character")
    assert_manifest_scalar(record$version, paste0("environment.package_locations[", index, "].version"), "character")
    if (!is.null(record$path)) assert_manifest_scalar(record$path, paste0("environment.package_locations[", index, "].path"), "character", nullable = TRUE)
  }
  invisible(manifest)
}

write_manifest_v2 <- function(manifest, path, physical_root = NULL) {
  validate_manifest_v2(manifest, physical_root = physical_root)
  if (exists("atomic_write_json", mode = "function")) {
    atomic_write_json(manifest, path)
  } else {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temp <- tempfile(pattern = ".manifest-", tmpdir = dirname(path))
    on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
    jsonlite::write_json(manifest, temp, pretty = TRUE, auto_unbox = TRUE, null = "null")
    if (file.exists(path)) unlink(path, force = TRUE)
    if (!file.rename(temp, path)) stop("Could not publish manifest atomically.", call. = FALSE)
  }
  parsed <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  validate_manifest_v2(parsed, physical_root = physical_root)
  invisible(path)
}
