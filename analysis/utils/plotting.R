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
