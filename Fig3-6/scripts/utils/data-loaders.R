#!/usr/bin/env Rscript
# Data loading utilities for Fig3-6 analysis pipeline
# Provides standardized data loading functions with consistent filtering

library(tidyverse)

#' Load rMATS differential splicing data with standard filtering
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param dpsi_threshold Absolute delta PSI threshold (default: 0.1)
#' @return Filtered tibble of splicing events
load_rmats_data <- function(fdr_threshold = 0.05, dpsi_threshold = 0.1) {
  data_path <- "data/rmats/all_events_ko_target_fdr_dpsi.csv"
  
  if (!file.exists(data_path)) {
    stop("rMATS data file not found: ", data_path)
  }
  
  read_csv(data_path, show_col_types = FALSE) %>%
    filter(fdr < fdr_threshold, abs(dpsi) > dpsi_threshold)
}

#' Load differential expression data
#' @return Tibble of DEG results
load_deg_data <- function() {
  data_path <- "data/degs/ko_vs_mieru.csv"
  
  if (!file.exists(data_path)) {
    stop("DEG data file not found: ", data_path)
  }
  
  read_csv(data_path, show_col_types = FALSE)
}

#' Load normalized DESeq2 counts
#' @return Matrix of normalized counts
load_normalized_counts <- function() {
  data_path <- "data/degs/deseq2_normalized_counts.csv"
  
  if (!file.exists(data_path)) {
    stop("Normalized counts file not found: ", data_path)
  }
  
  read_csv(data_path, show_col_types = FALSE)
}

#' Load exon characteristics data
#' @return Tibble of exon features
load_exon_data <- function() {
  length_gc_path <- "data/exon_characteristics/exon_length_gc.tsv"
  conservation_path <- "data/exon_characteristics/exon_conservation.tsv"
  
  exon_data <- list()
  
  if (file.exists(length_gc_path)) {
    exon_data$length_gc <- read_tsv(length_gc_path, show_col_types = FALSE)
  }
  
  if (file.exists(conservation_path)) {
    exon_data$conservation <- read_tsv(conservation_path, show_col_types = FALSE)
  }
  
  return(exon_data)
}

#' Get list of RBP knockout genes
#' @return Character vector of KO gene symbols
get_ko_genes <- function() {
  c("Cd2bp2", "Qk", "Rbm24", "Rpl22l1", "Spen", "Strap", 
    "Tra2b", "Trim71", "Ubr5", "Wt1", "Ybx1")
}

#' Get event type labels
#' @return Named character vector of event types
get_event_types <- function() {
  c(
    "A3SS" = "Alternative 3' Splice Site",
    "A5SS" = "Alternative 5' Splice Site", 
    "MXE" = "Mutually Exclusive Exons",
    "RI" = "Retained Intron",
    "SE" = "Skipped Exon"
  )
}