#!/usr/bin/env Rscript
# Plot themes and styling utilities for Fig3-6 analysis pipeline
# Provides standardized ggplot2 themes and color schemes

library(ggplot2)
library(extrafont)

#' Get standard publication-ready ggplot2 theme
#' @param base_size Base font size (default: 14)
#' @return ggplot2 theme object
get_standard_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Arial"),
      axis.title.x = element_text(size = base_size + 4),
      axis.title.y = element_text(size = base_size + 4),
      axis.text.x = element_text(size = base_size),
      axis.text.y = element_text(size = base_size),
      plot.title = element_text(size = base_size + 6, hjust = 0.5),
      legend.title = element_text(size = base_size + 2),
      legend.text = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = base_size + 2)
    )
}

#' Get standard color palette for splicing events
#' @return Named character vector of colors
get_event_colors <- function() {
  colors <- c("#44ED8B", "#FF2FC1", "#3FAFFF", "#FFE270", "#FF604E")
  names(colors) <- c("A3SS", "A5SS", "MXE", "RI", "SE")
  return(colors)
}

#' Get color palette for KO genes
#' @return Character vector of colors for different KO genes
get_ko_colors <- function() {
  RColorBrewer::brewer.pal(11, "Spectral")
}

#' Get standard plot dimensions
#' @param plot_type Type of plot ("barplot", "violin", "heatmap", "venn")
#' @return List with width and height
get_plot_dimensions <- function(plot_type = "standard") {
  dimensions <- switch(plot_type,
    "barplot" = list(width = 8, height = 6),
    "violin" = list(width = 10, height = 8),
    "heatmap" = list(width = 12, height = 10),
    "venn" = list(width = 8, height = 8),
    "circos" = list(width = 10, height = 10),
    "wide_violin" = list(width = 24, height = 6),
    list(width = 8, height = 6)  # default
  )
  return(dimensions)
}

#' Save plots in standard formats
#' @param plot ggplot object
#' @param filename Base filename without extension
#' @param output_dir Output directory
#' @param plot_type Plot type for dimension selection
#' @param formats Vector of formats to save ("pdf", "jpg", "svg")
save_standard_plot <- function(plot, filename, output_dir, 
                              plot_type = "standard", 
                              formats = c("pdf", "jpg")) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get standard dimensions
  dims <- get_plot_dimensions(plot_type)
  
  # Save in requested formats
  for (format in formats) {
    output_path <- file.path(output_dir, paste0(filename, ".", format))
    
    if (format == "pdf") {
      ggsave(output_path, plot, 
             width = dims$width, height = dims$height, 
             family = "Arial", device = cairo_pdf)
    } else if (format == "svg") {
      ggsave(output_path, plot, 
             width = dims$width, height = dims$height,
             device = "svg")
    } else {
      ggsave(output_path, plot, 
             width = dims$width, height = dims$height)
    }
  }
  
  message("Saved plot: ", filename, " in formats: ", paste(formats, collapse = ", "))
}

#' Create standard violin plot
#' @param data Data frame
#' @param x_var X-axis variable name
#' @param y_var Y-axis variable name
#' @param fill_var Fill variable name (optional)
#' @param title Plot title
#' @return ggplot object
create_violin_plot <- function(data, x_var, y_var, fill_var = NULL, title = "") {
  p <- ggplot(data, aes_string(x = x_var, y = y_var))
  
  if (!is.null(fill_var)) {
    p <- p + aes_string(fill = fill_var)
  }
  
  p <- p +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(title = title) +
    get_standard_theme()
  
  return(p)
}