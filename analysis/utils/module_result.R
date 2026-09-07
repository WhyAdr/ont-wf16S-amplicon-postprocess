# =============================================================================
# Module result-envelope validation
# =============================================================================

validate_module_result <- function(result, module_name, output_root) {
  if (!is.list(result) || is.null(names(result))) {
    stop(sprintf("Module '%s' returned a non-object result.", module_name), call. = FALSE)
  }
  status <- result$status
  if (!is.character(status) || length(status) != 1L || is.na(status) ||
      !status %in% c("completed", "skipped", "failed")) {
    stop(sprintf("Module '%s' returned an invalid status.", module_name), call. = FALSE)
  }
  outputs <- result$outputs %||% character(0)
  if (!is.character(outputs) || anyNA(outputs) || anyDuplicated(outputs)) {
    stop(sprintf("Module '%s' returned an invalid output list.", module_name), call. = FALSE)
  }
  root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
  if (length(outputs)) {
    normalized <- normalizePath(outputs, winslash = "/", mustWork = FALSE)
    contained <- startsWith(tolower(normalized), paste0(tolower(root), "/"))
    valid <- file.exists(normalized) & !dir.exists(normalized) & contained
    if (!all(valid)) stop(sprintf("Module '%s' declared missing or unowned output '%s'.",
                                  module_name, outputs[which(!valid)[1]]), call. = FALSE)
    result$outputs <- normalized
  }
  if (status == "completed" && !is.null(result$error)) stop(sprintf("Module '%s' completed with an error.", module_name), call. = FALSE)
  if (status == "failed" && (is.null(result$error) || !nzchar(result$error))) stop(sprintf("Module '%s' failed without an error.", module_name), call. = FALSE)
  if (status == "skipped" && (is.null(result$reason) || !nzchar(result$reason))) stop(sprintf("Module '%s' skipped without a reason.", module_name), call. = FALSE)
  result
}
