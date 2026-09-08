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
  if (status == "failed" && length(outputs) > 0L) {
    stop(sprintf("Module '%s' failed but declared outputs.", module_name), call. = FALSE)
  }
  warnings <- result$warnings %||% character(0)
  if (!is.character(warnings) || anyNA(warnings)) {
    stop(sprintf("Module '%s' returned invalid warnings: must be a character vector with no NA values.", module_name),
         call. = FALSE)
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
  # Physical stage census check: require physical files to equal declared outputs
  physical_files <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  physical_files <- physical_files[!dir.exists(physical_files)]
  if (length(physical_files)) {
    norm_phys <- tolower(normalizePath(physical_files, winslash = "/", mustWork = FALSE))
    norm_decl <- if (length(outputs)) tolower(normalizePath(outputs, winslash = "/", mustWork = FALSE)) else character(0)
    extra_files <- setdiff(norm_phys, norm_decl)
    if (length(extra_files)) {
      stop(sprintf("Module '%s' produced undeclared output(s) in staging directory: %s",
                   module_name, paste(extra_files, collapse = ", ")), call. = FALSE)
    }
  }
  if (status == "completed" && !is.null(result$error)) stop(sprintf("Module '%s' completed with an error.", module_name), call. = FALSE)
  if (status == "failed" && (is.null(result$error) || !nzchar(result$error))) stop(sprintf("Module '%s' failed without an error.", module_name), call. = FALSE)
  if (status == "skipped" && (is.null(result$reason) || !nzchar(result$reason))) stop(sprintf("Module '%s' skipped without a reason.", module_name), call. = FALSE)
  result
}
