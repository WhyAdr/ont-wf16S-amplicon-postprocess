# =============================================================================
# Plotting Helpers, Theme, and Color Palettes
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
})

theme_amplicon <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(color = "grey40", size = base_size - 2, margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size)
    )
}

# Explicitly bounded palette for the new classified-read composition figures.
# The order is part of the output contract: the first named taxon receives the
# first colour and the mapping is reused across its figures and sidecars.
composition_palette <- c(
  "#FF8C5A", "#E6B66E", "#8DB17E", "#5C9462", "#2F6B34",
  "#34C79B", "#63C4C7", "#3699C5", "#5F65C8", "#9060C2",
  "#BF5BC4", "#ED54B8", "#F184B9", "#C69090", "#996666"
)

composition_colors <- function(display_paths, display_labels = NULL) {
  if (is.null(display_labels)) display_labels <- display_paths
  if (!is.character(display_paths) || !is.character(display_labels) ||
      length(display_paths) != length(display_labels) ||
      anyNA(display_paths) || anyNA(display_labels) ||
      anyDuplicated(display_paths) || anyDuplicated(display_labels)) {
    stop("Composition colour mapping requires unique path and label vectors.", call. = FALSE)
  }
  named <- display_paths != "__OTHER__"
  if (sum(named) > length(composition_palette)) {
    stop(sprintf("Composition figures support at most %d named taxa.",
                 length(composition_palette)), call. = FALSE)
  }
  values <- character(length(display_paths))
  values[named] <- composition_palette[seq_len(sum(named))]
  values[!named] <- "#F3C58F"
  stats::setNames(values, display_labels)
}

phylum_palette <- c(
  "Bacillota"                  = "#1b9e77",
  "Pseudomonadota"             = "#d95f02",
  "Bacteroidota"               = "#7570b3",
  "Actinomycetota"             = "#e7298a",
  "Fusobacteriota"             = "#66a61e",
  "Spirochaetota"              = "#e6ab02",
  "Synergistota"               = "#a6761d",
  "Planctomycetota"            = "#4e79a7",
  "Campylobacterota"           = "#e15759",
  "Chloroflexota"              = "#76b7b2",
  "Candidatus Melainabacteria"  = "#59a14f",
  "Ignavibacteriota"           = "#edc949",
  "Armatimonadota"             = "#af7aa1",
  "Verrucomicrobiota"          = "#666666",
  "Unclassified"               = "#bdbdbd",
  "Other phyla"                = "#41566b",
  "Other"                      = "#8c96a0"
)

get_phylum_colors <- function(phyla) {
  cols <- phylum_palette[phyla]
  missing_idx <- is.na(cols)
  if (any(missing_idx)) {
    # Generate fallback distinct colors
    n_miss <- sum(missing_idx)
    cols[missing_idx] <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_miss)
  }
  names(cols) <- phyla
  cols
}

save_plot <- function(filename, plot, width = 7, height = 5, dpi = 150) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
}

with_png_device <- function(filename, draw, width = 7, height = 5, dpi = 150) {
  if (!is.function(draw)) stop("'draw' must be a function.", call. = FALSE)
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
  device_id <- grDevices::dev.cur()
  completed <- FALSE
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(which = device_id)
    }
    if (!completed && file.exists(filename)) unlink(filename)
  }, add = TRUE)

  draw()
  grDevices::dev.off(which = device_id)
  if (!file.exists(filename) || !isTRUE(file.info(filename)$size > 0)) {
    stop(sprintf("PNG output was not written or is empty: '%s'.", filename), call. = FALSE)
  }
  completed <- TRUE
  invisible(filename)
}
