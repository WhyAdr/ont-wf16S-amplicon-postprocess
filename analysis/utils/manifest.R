# =============================================================================
# Manifest Schema v2 Helpers
# =============================================================================

json_array <- function(x) {
  if (is.null(x) || length(x) == 0L) return(list())
  unname(lapply(as.list(x), function(value) {
    if (is.list(value)) value else jsonlite::unbox(value)
  }))
}

is_json_array <- function(x) {
  is.list(x) && is.null(names(x))
}

assert_manifest_array <- function(value, path) {
  if (!is_json_array(value)) {
    stop(sprintf("Manifest schema v2 violation at '%s': expected a JSON array.", path),
         call. = FALSE)
  }
  invisible(TRUE)
}

assert_manifest_scalar_or_null <- function(value, path) {
  if (!is.null(value) && is.list(value)) {
    stop(sprintf("Manifest schema v2 violation at '%s': expected a scalar or null.", path),
         call. = FALSE)
  }
  invisible(TRUE)
}

new_not_run_module_record <- function(reason = NULL) {
  list(
    status = "not_run",
    outputs = character(0),
    warnings = character(0),
    error = NULL,
    reason = reason,
    start_time = NULL,
    end_time = NULL,
    duration_seconds = NULL
  )
}

manifest_module_record <- function(record) {
  list(
    status = record$status,
    outputs = json_array(record$outputs %||% character(0)),
    warnings = json_array(record$warnings %||% character(0)),
    error = record$error %||% NULL,
    reason = record$reason %||% NULL,
    start_time = record$start_time %||% NULL,
    end_time = record$end_time %||% NULL,
    duration_seconds = record$duration_seconds %||% NULL
  )
}

validate_manifest_v2 <- function(manifest) {
  if (!is.list(manifest) || is.null(names(manifest))) {
    stop("Manifest schema v2 violation at '<root>': expected an object.", call. = FALSE)
  }
  required <- c(
    "schema_version", "config_schema_version", "samples", "command", "cli",
    "inputs", "modules", "warnings", "package_versions"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0L) {
    stop(sprintf("Manifest schema v2 violation at '<root>': missing required key(s): %s.",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!identical(manifest$schema_version, 2L)) {
    stop("Manifest schema v2 violation at 'schema_version': expected 2.", call. = FALSE)
  }
  if (!identical(manifest$config_schema_version, 1L)) {
    stop("Manifest schema v2 violation at 'config_schema_version': expected 1.", call. = FALSE)
  }
  for (field in c("samples", "command", "warnings", "package_versions")) {
    assert_manifest_array(manifest[[field]], field)
  }
  if (!is.list(manifest$cli) || is.null(names(manifest$cli))) {
    stop("Manifest schema v2 violation at 'cli': expected an object.", call. = FALSE)
  }
  if (is.null(manifest$cli$modules)) {
    stop("Manifest schema v2 violation at 'cli.modules': missing required key.", call. = FALSE)
  }
  assert_manifest_array(manifest$cli$modules, "cli.modules")
  if (!is.list(manifest$inputs) || is.null(names(manifest$inputs))) {
    stop("Manifest schema v2 violation at 'inputs': expected an object.", call. = FALSE)
  }
  for (field in c("assignments", "bamstats")) {
    value <- manifest$inputs[[field]]
    if (!is.null(value)) assert_manifest_array(value, paste0("inputs.", field))
  }
  if (!is.list(manifest$modules) || is.null(names(manifest$modules))) {
    stop("Manifest schema v2 violation at 'modules': expected an object.", call. = FALSE)
  }
  valid_statuses <- c("completed", "skipped", "failed", "not_run")
  for (module_name in names(manifest$modules)) {
    path <- paste0("modules.", module_name)
    record <- manifest$modules[[module_name]]
    if (!is.list(record) || is.null(names(record))) {
      stop(sprintf("Manifest schema v2 violation at '%s': expected an object.", path), call. = FALSE)
    }
    fields <- c("status", "outputs", "warnings", "error", "reason", "start_time", "end_time", "duration_seconds")
    missing_fields <- setdiff(fields, names(record))
    if (length(missing_fields) > 0L) {
      stop(sprintf("Manifest schema v2 violation at '%s': missing key(s): %s.", path,
                   paste(missing_fields, collapse = ", ")), call. = FALSE)
    }
    if (!is.character(record$status) || length(record$status) != 1L ||
        is.na(record$status) || !record$status %in% valid_statuses) {
      stop(sprintf("Manifest schema v2 violation at '%s.status': invalid module status.", path),
           call. = FALSE)
    }
    assert_manifest_array(record$outputs, paste0(path, ".outputs"))
    assert_manifest_array(record$warnings, paste0(path, ".warnings"))
    for (field in c("error", "reason", "start_time", "end_time", "duration_seconds")) {
      assert_manifest_scalar_or_null(record[[field]], paste0(path, ".", field))
    }
    if (identical(record$status, "not_run")) {
      if (length(record$outputs) != 0L || length(record$warnings) != 0L ||
          !is.null(record$error) || !is.null(record$start_time) ||
          !is.null(record$end_time) || !is.null(record$duration_seconds)) {
        stop(sprintf("Manifest schema v2 violation at '%s': invalid not_run record.", path),
             call. = FALSE)
      }
    }
  }
  invisible(manifest)
}

write_manifest_v2 <- function(manifest, path) {
  validate_manifest_v2(manifest)
  jsonlite::write_json(manifest, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(path)
}
