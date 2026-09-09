# =============================================================================
# Byte-strict pipeline version helpers
# =============================================================================

read_strict_semver <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !file.exists(path) ||
      dir.exists(path)) {
    stop(sprintf("VERSION is missing: '%s'.", path), call. = FALSE)
  }
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  text <- rawToChar(bytes)
  if (!grepl("^[0-9]+[.]([0-9]+)[.]([0-9]+)\n$", text)) {
    stop("VERSION must contain exactly one newline-terminated SemVer value.", call. = FALSE)
  }
  sub("\n$", "", text)
}

read_pipeline_version <- read_strict_semver
