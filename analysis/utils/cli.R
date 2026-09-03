# =============================================================================
# CLI Parser Utility
# =============================================================================

suppressMessages(library(optparse))

get_cli_parser <- function() {
  option_list <- list(
    optparse::make_option(
      c("-c", "--config"),
      type = "character",
      default = "config.yml",
      dest = "config",
      help = "Path to configuration YAML file [default: %default]"
    ),
    optparse::make_option(
      c("-o", "--output-dir"),
      type = "character",
      default = NULL,
      dest = "output_dir",
      help = "Override output base directory"
    ),
    optparse::make_option(
      c("-m", "--modules"),
      type = "character",
      default = NULL,
      dest = "modules",
      help = "Comma-separated list of modules to run (e.g. qc,alpha,composition,kreport)"
    ),
    optparse::make_option(
      c("--validate-only"),
      action = "store_true",
      default = FALSE,
      dest = "validate_only",
      help = "Perform config, dependency, schema, and reconciliation checks without writing outputs"
    ),
    optparse::make_option(
      c("--keep-going"),
      action = "store_true",
      default = FALSE,
      dest = "keep_going",
      help = "Continue executing remaining modules if an individual module fails"
    ),
    optparse::make_option(
      c("--refresh-taxonomy"),
      action = "store_true",
      default = FALSE,
      dest = "refresh_taxonomy",
      help = "Allow online NCBI queries to resolve missing TaxIDs and refresh cache"
    ),
    optparse::make_option(
      c("--overwrite"),
      action = "store_true",
      default = FALSE,
      dest = "overwrite",
      help = "Allow overwriting existing output files"
    )
  )
  
  optparse::OptionParser(
    usage = "%prog [options]",
    description = "ONT wf-16s Post-Processing Analytical Pipeline",
    option_list = option_list
  )
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  parser <- get_cli_parser()
  optparse::parse_args(parser, args = args)
}
