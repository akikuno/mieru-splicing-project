#!/usr/bin/env Rscript
# Configuration loading utilities
# Provides functions to load and validate configuration parameters

library(yaml)

#' Load configuration from YAML file
#' @param config_path Path to configuration file
#' @return List with configuration parameters
load_config <- function(config_path = "scripts/config/parameters.yaml") {
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  
  config <- read_yaml(config_path)
  
  # Validate required sections
  required_sections <- c("thresholds", "plotting", "paths", "analysis")
  missing_sections <- setdiff(required_sections, names(config))
  
  if (length(missing_sections) > 0) {
    stop("Missing required configuration sections: ", paste(missing_sections, collapse = ", "))
  }
  
  return(config)
}

#' Get output directory for a specific figure
#' @param fig_number Figure number (3, 4, 5, 6, 7, or "SFig")
#' @param config Configuration list (optional)
#' @return Character string with output directory path
get_output_dir <- function(fig_number, config = NULL) {
  if (is.null(config)) {
    config <- load_config()
  }
  
  fig_key <- paste0("fig", tolower(fig_number))
  
  if (fig_key %in% names(config$paths$reports)) {
    return(config$paths$reports[[fig_key]])
  } else {
    stop("Unknown figure number: ", fig_number)
  }
}

#' Create all output directories
#' @param config Configuration list (optional)
create_output_dirs <- function(config = NULL) {
  if (is.null(config)) {
    config <- load_config()
  }
  
  # Create all report directories
  for (dir_path in config$paths$reports) {
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  }
  
  cat("Created output directories:\n")
  cat(paste(unlist(config$paths$reports), collapse = "\n"))
  cat("\n")
}