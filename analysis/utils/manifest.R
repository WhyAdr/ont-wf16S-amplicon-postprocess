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

validate_manifest_v2 <- function(manifest, physical_root = NULL) {
  if (!is.list(manifest) || is.null(names(manifest))) manifest_fail("<root>", "expected an object.")
  required <- c("pipeline", "pipeline_version", "schema_version", "schema_revision",
    "config_schema_version", "run_status", "start_time", "end_time", "duration_seconds",
    "samples", "command", "cli", "inputs", "modules", "owned_outputs", "warnings",
    "package_versions", "environment", "output_root")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) manifest_fail("<root>", paste("missing required key(s):", paste(missing, collapse = ", ")))
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
  for (field in c("samples", "command", "owned_outputs", "warnings")) {
    assert_manifest_array(manifest[[field]], field, "character")
    values <- array_values(manifest[[field]])
    if (anyDuplicated(values)) manifest_fail(field, "array elements must be unique.")
  }
  assert_manifest_array(manifest$package_versions, "package_versions")
  packages <- vapply(manifest$package_versions, function(x) as.character(x$package %||% NA), character(1))
  if (anyNA(packages) || anyDuplicated(packages)) manifest_fail("package_versions", "package names must be present and unique.")
  if (!is.list(manifest$cli) || is.null(names(manifest$cli))) manifest_fail("cli", "expected an object.")
  assert_manifest_array(manifest$cli$modules, "cli.modules", "character")
  requested <- array_values(manifest$cli$modules)
  if (anyDuplicated(requested)) manifest_fail("cli.modules", "module names must be unique.")
  if (!is.list(manifest$inputs) || is.null(names(manifest$inputs))) manifest_fail("inputs", "expected an object.")
  for (field in c("assignments", "bamstats")) assert_manifest_array(manifest$inputs[[field]], paste0("inputs.", field))
  for (field in c("assignments", "bamstats")) {
    ids <- vapply(manifest$inputs[[field]], function(x) as.character(x$sample_id %||% NA), character(1))
    if (anyNA(ids) || anyDuplicated(ids)) manifest_fail(paste0("inputs.", field), "sample IDs must be present and unique.")
  }
  if (!is.list(manifest$modules) || is.null(names(manifest$modules)) || anyDuplicated(names(manifest$modules))) {
    manifest_fail("modules", "expected a uniquely named object.")
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
  for (output in all_outputs) {
    normalized <- normalizePath(output, winslash = "/", mustWork = FALSE)
    if (!startsWith(tolower(normalized), paste0(tolower(root), "/"))) manifest_fail("modules.*.outputs", "output is outside output_root.")
    relative <- substring(normalized, nchar(root) + 2L)
    physical <- file.path(check_root, relative)
    if (!file.exists(physical) || dir.exists(physical)) manifest_fail("modules.*.outputs", paste("missing regular file", output))
  }
  env <- manifest$environment
  if (!is.list(env) || !env$lock_status %in% c("synchronized", "mismatch", "missing", "not_checked") ||
      !is.logical(env$locked) || length(env$locked) != 1L || is.na(env$locked) ||
      !identical(env$locked, identical(env$lock_status, "synchronized"))) manifest_fail("environment", "inconsistent lock status.")
  assert_manifest_array(env$r_discrepancies, "environment.r_discrepancies", "character")
  assert_manifest_array(env$package_discrepancies, "environment.package_discrepancies", "character")
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
