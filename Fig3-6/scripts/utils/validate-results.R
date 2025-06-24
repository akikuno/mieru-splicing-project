#!/usr/bin/env Rscript
# Validate analysis results for reproducibility and quality control

library(tidyverse)
library(yaml)

# Load configuration
source("scripts/utils/data-loaders.R")
config <- read_yaml("scripts/config/parameters.yaml")

cat("======================================\n")
cat("Validating Analysis Results\n")
cat("======================================\n")

# Create validation report directory
dir.create("reports/validation", showWarnings = FALSE, recursive = TRUE)

# Initialize validation results
validation_results <- list()

# 1. Validate splicing data
cat("\n1. Validating splicing event data...\n")
if (file.exists("data/rmats/all_events_ko_target_fdr_dpsi.csv")) {
  splicing_data <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv", show_col_types = FALSE)
  
  # Basic statistics
  total_events <- nrow(splicing_data)
  significant_events <- splicing_data %>% 
    filter(fdr < config$thresholds$fdr, abs(dpsi) > config$thresholds$dpsi) %>%
    nrow()
  
  # Validate event types
  expected_events <- config$analysis$event_types
  observed_events <- unique(splicing_data$event)
  missing_events <- setdiff(expected_events, observed_events)
  unexpected_events <- setdiff(observed_events, expected_events)
  
  # Validate KO genes
  expected_kos <- config$analysis$ko_genes
  observed_kos <- unique(splicing_data$ko_symbol)
  missing_kos <- setdiff(expected_kos, observed_kos)
  unexpected_kos <- setdiff(observed_kos, expected_kos)
  
  validation_results$splicing <- list(
    total_events = total_events,
    significant_events = significant_events,
    significance_rate = significant_events / total_events * 100,
    missing_event_types = missing_events,
    unexpected_event_types = unexpected_events,
    missing_ko_genes = missing_kos,
    unexpected_ko_genes = unexpected_kos,
    status = "CHECKED"
  )
  
  cat(sprintf("  ✓ Total events: %d\n", total_events))
  cat(sprintf("  ✓ Significant events: %d (%.1f%%)\n", significant_events, significant_events / total_events * 100))
  
  if (length(missing_events) > 0) {
    cat(sprintf("  ⚠ Missing event types: %s\n", paste(missing_events, collapse = ", ")))
  }
  if (length(missing_kos) > 0) {
    cat(sprintf("  ⚠ Missing KO genes: %s\n", paste(missing_kos, collapse = ", ")))
  }
} else {
  validation_results$splicing <- list(status = "NOT_FOUND")
  cat("  ❌ Splicing data file not found\n")
}

# 2. Validate DEG data
cat("\n2. Validating differential expression data...\n")
if (file.exists("data/degs/ko_vs_mieru.csv")) {
  deg_data <- read_csv("data/degs/ko_vs_mieru.csv", show_col_types = FALSE)
  
  total_genes <- nrow(deg_data)
  significant_degs <- deg_data %>% 
    filter(padj < config$thresholds$padj) %>%
    nrow()
  
  validation_results$degs <- list(
    total_genes = total_genes,
    significant_degs = significant_degs,
    significance_rate = significant_degs / total_genes * 100,
    status = "CHECKED"
  )
  
  cat(sprintf("  ✓ Total genes tested: %d\n", total_genes))
  cat(sprintf("  ✓ Significant DEGs: %d (%.1f%%)\n", significant_degs, significant_degs / total_genes * 100))
} else {
  validation_results$degs <- list(status = "NOT_FOUND")
  cat("  ❌ DEG data file not found\n")
}

# 3. Validate QC metrics
cat("\n3. Validating quality control metrics...\n")
if (file.exists("reports/QC/qc_summary_statistics.csv")) {
  qc_data <- read_csv("reports/QC/qc_summary_statistics.csv", show_col_types = FALSE)
  
  # Check mapping rates
  mean_mapping_rate <- qc_data$mean_mapping_rate
  samples_below_threshold <- qc_data$samples_below_80pct
  
  validation_results$qc <- list(
    mean_mapping_rate = mean_mapping_rate,
    samples_below_80pct = samples_below_threshold,
    mapping_quality = ifelse(mean_mapping_rate >= 80, "GOOD", "WARNING"),
    status = "CHECKED"
  )
  
  cat(sprintf("  ✓ Mean mapping rate: %.1f%%\n", mean_mapping_rate))
  if (samples_below_threshold > 0) {
    cat(sprintf("  ⚠ %d samples below 80%% mapping rate\n", samples_below_threshold))
  }
} else {
  validation_results$qc <- list(status = "NOT_FOUND")
  cat("  ❌ QC summary file not found\n")
}

# 4. Validate output files
cat("\n4. Validating output files...\n")
expected_outputs <- c(
  "reports/Fig3/01-event-frequency-barplot.pdf",
  "reports/Fig3/02-dpsi-distribution-violin.pdf",
  "reports/Fig4/025-venn_overlap_dsg_deg.pdf",
  "reports/QC/mapping-rates.pdf"
)

output_validation <- map_lgl(expected_outputs, file.exists)
names(output_validation) <- expected_outputs

validation_results$outputs <- list(
  expected_files = length(expected_outputs),
  found_files = sum(output_validation),
  missing_files = names(output_validation)[!output_validation],
  status = ifelse(all(output_validation), "COMPLETE", "INCOMPLETE")
)

cat(sprintf("  ✓ Expected output files: %d\n", length(expected_outputs)))
cat(sprintf("  ✓ Found output files: %d\n", sum(output_validation)))

if (any(!output_validation)) {
  cat("  ⚠ Missing output files:\n")
  for (missing_file in names(output_validation)[!output_validation]) {
    cat(sprintf("    - %s\n", missing_file))
  }
}

# 5. Statistical validation
cat("\n5. Performing statistical validation...\n")
if (validation_results$splicing$status == "CHECKED" && validation_results$degs$status == "CHECKED") {
  
  # Check if results are within expected ranges
  expected_significance_rate <- 5  # 5% is typical for FDR 0.05
  actual_significance_rate <- validation_results$splicing$significance_rate
  
  statistical_checks <- list(
    significance_rate_reasonable = actual_significance_rate >= 1 && actual_significance_rate <= 20,
    sufficient_events = validation_results$splicing$total_events >= 1000,
    sufficient_genes = validation_results$degs$total_genes >= 10000
  )
  
  validation_results$statistics <- statistical_checks
  
  cat(sprintf("  ✓ Significance rate check: %s (%.1f%%)\n", 
              ifelse(statistical_checks$significance_rate_reasonable, "PASS", "WARNING"),
              actual_significance_rate))
  cat(sprintf("  ✓ Event count check: %s (%d events)\n",
              ifelse(statistical_checks$sufficient_events, "PASS", "WARNING"),
              validation_results$splicing$total_events))
  cat(sprintf("  ✓ Gene count check: %s (%d genes)\n",
              ifelse(statistical_checks$sufficient_genes, "PASS", "WARNING"),
              validation_results$degs$total_genes))
}

# 6. Generate validation report
cat("\n6. Generating validation report...\n")

# Overall status
overall_status <- "PASS"
if (validation_results$splicing$status == "NOT_FOUND" || 
    validation_results$degs$status == "NOT_FOUND") {
  overall_status <- "FAIL"
} else if (validation_results$outputs$status == "INCOMPLETE" ||
           validation_results$qc$status == "NOT_FOUND") {
  overall_status <- "WARNING"
}

# Save validation results
validation_summary <- list(
  timestamp = Sys.time(),
  overall_status = overall_status,
  validation_results = validation_results,
  thresholds_used = config$thresholds,
  session_info = sessionInfo()
)

# Save as YAML for human readability
write_yaml(validation_summary, "reports/validation/validation_summary.yml")

# Save as RDS for R compatibility
saveRDS(validation_summary, "reports/validation/validation_summary.rds")

# Generate human-readable report
report_lines <- c(
  "Analysis Validation Report",
  "=========================",
  paste("Generated:", Sys.time()),
  paste("Overall Status:", overall_status),
  "",
  "Summary:",
  sprintf("- Splicing events: %s", validation_results$splicing$status),
  sprintf("- DEG analysis: %s", validation_results$degs$status),
  sprintf("- Quality control: %s", validation_results$qc$status),
  sprintf("- Output files: %s", validation_results$outputs$status),
  "",
  "Details available in validation_summary.yml"
)

writeLines(report_lines, "reports/validation/validation_report.txt")

cat("\n======================================\n")
cat("Validation Complete\n")
cat("======================================\n")
cat(sprintf("Overall Status: %s\n", overall_status))
cat("Reports saved to: reports/validation/\n")

if (overall_status == "FAIL") {
  cat("\n❌ Validation FAILED - Critical issues found\n")
  quit(status = 1)
} else if (overall_status == "WARNING") {
  cat("\n⚠ Validation completed with WARNINGS\n")
} else {
  cat("\n✅ Validation PASSED - All checks successful\n")
}