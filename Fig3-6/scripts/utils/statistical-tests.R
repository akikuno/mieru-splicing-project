#!/usr/bin/env Rscript
# Statistical testing utilities for Fig3-6 analysis pipeline
# Provides standardized statistical analysis functions

library(tidyverse)

#' Perform Fisher's exact test for enrichment analysis
#' @param genes_of_interest Character vector of genes of interest
#' @param genes_in_category Character vector of genes in category
#' @param background_genes Character vector of background genes
#' @return List with test results
fisher_enrichment_test <- function(genes_of_interest, genes_in_category, background_genes) {
  
  # Create contingency table
  genes_in_both <- intersect(genes_of_interest, genes_in_category)
  genes_interest_only <- setdiff(genes_of_interest, genes_in_category)
  genes_category_only <- setdiff(genes_in_category, genes_of_interest)
  genes_neither <- setdiff(background_genes, union(genes_of_interest, genes_in_category))
  
  # Contingency table:
  #              In Category    Not in Category
  # Of Interest      a               b
  # Not Interest     c               d
  
  a <- length(genes_in_both)
  b <- length(genes_interest_only)
  c <- length(genes_category_only)
  d <- length(genes_neither)
  
  # Create matrix for fisher.test
  contingency_matrix <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
  
  # Calculate odds ratio and confidence interval
  odds_ratio <- fisher_result$estimate
  ci_lower <- fisher_result$conf.int[1]
  ci_upper <- fisher_result$conf.int[2]
  p_value <- fisher_result$p.value
  
  # Return results
  list(
    genes_in_both = genes_in_both,
    n_in_both = a,
    n_interest = a + b,
    n_category = a + c,
    n_background = a + b + c + d,
    odds_ratio = odds_ratio,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    contingency_matrix = contingency_matrix
  )
}

#' Perform multiple Fisher's tests with FDR correction
#' @param gene_lists Named list of gene vectors to test
#' @param category_genes Character vector of genes in category
#' @param background_genes Character vector of background genes
#' @return Tibble with test results
multiple_fisher_tests <- function(gene_lists, category_genes, background_genes) {
  
  # Perform Fisher's tests for each gene list
  results <- map_dfr(names(gene_lists), function(list_name) {
    test_result <- fisher_enrichment_test(
      gene_lists[[list_name]], 
      category_genes, 
      background_genes
    )
    
    tibble(
      comparison = list_name,
      n_overlap = test_result$n_in_both,
      n_interest = test_result$n_interest,
      n_category = test_result$n_category,
      odds_ratio = test_result$odds_ratio,
      ci_lower = test_result$ci_lower,
      ci_upper = test_result$ci_upper,
      p_value = test_result$p_value
    )
  })
  
  # Add FDR correction
  results$fdr <- p.adjust(results$p_value, method = "fdr")
  
  return(results)
}

#' Calculate percentage overlap between gene sets
#' @param gene_sets Named list of gene vectors
#' @return Matrix of percentage overlaps
calculate_overlap_matrix <- function(gene_sets) {
  
  set_names <- names(gene_sets)
  n_sets <- length(gene_sets)
  
  # Initialize matrix
  overlap_matrix <- matrix(0, nrow = n_sets, ncol = n_sets)
  rownames(overlap_matrix) <- set_names
  colnames(overlap_matrix) <- set_names
  
  # Calculate overlaps
  for (i in 1:n_sets) {
    for (j in 1:n_sets) {
      if (i == j) {
        overlap_matrix[i, j] <- 100  # Self-overlap is 100%
      } else {
        overlap_size <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
        overlap_matrix[i, j] <- (overlap_size / length(gene_sets[[i]])) * 100
      }
    }
  }
  
  return(overlap_matrix)
}

#' Perform hypergeometric test for pathway enrichment
#' @param genes_of_interest Character vector of genes of interest
#' @param pathway_genes Character vector of genes in pathway
#' @param universe_size Total number of genes in universe
#' @return List with test results
hypergeometric_test <- function(genes_of_interest, pathway_genes, universe_size) {
  
  # Calculate parameters
  k <- length(intersect(genes_of_interest, pathway_genes))  # successes in sample
  m <- length(pathway_genes)  # successes in universe
  n <- universe_size - m  # failures in universe  
  q <- length(genes_of_interest)  # sample size
  
  # Perform hypergeometric test
  p_value <- phyper(k - 1, m, n, q, lower.tail = FALSE)
  
  # Calculate fold enrichment
  expected <- (length(genes_of_interest) * length(pathway_genes)) / universe_size
  fold_enrichment <- k / expected
  
  list(
    observed = k,
    expected = expected,
    fold_enrichment = fold_enrichment,
    p_value = p_value,
    genes_in_pathway = intersect(genes_of_interest, pathway_genes)
  )
}