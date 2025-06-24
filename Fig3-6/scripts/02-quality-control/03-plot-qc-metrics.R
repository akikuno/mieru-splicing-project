#!/usr/bin/env Rscript
# Generate quality control plots and summary statistics

# Load required libraries
source("scripts/utils/data-loaders.R")
source("scripts/utils/plot-themes.R")
library(yaml)

# Load configuration
config <- read_yaml("scripts/config/parameters.yaml")

cat("Generating QC plots and metrics...\n")

# Create output directory
qc_dir <- "reports/QC"
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# Check if mapping summary exists
mapping_summary_path <- file.path(qc_dir, "mapping", "mapping_summary.txt")

if (!file.exists(mapping_summary_path)) {
  cat("Warning: Mapping summary not found. Run mapping stats first.\n")
} else {
  # Load mapping statistics
  mapping_stats <- read_csv(mapping_summary_path, show_col_types = FALSE)
  
  # Clean mapping rate column (remove % sign)
  mapping_stats$Mapping_Rate_Numeric <- as.numeric(gsub("%", "", mapping_stats$Mapping_Rate))
  
  # Create mapping rate plot
  p_mapping <- ggplot(mapping_stats, aes(x = reorder(Sample, Mapping_Rate_Numeric), y = Mapping_Rate_Numeric)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = "Mapping Rate by Sample",
      x = "Sample",
      y = "Mapping Rate (%)",
      caption = "Red dashed line indicates 80% threshold"
    ) +
    get_standard_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save mapping rate plot
  save_standard_plot(
    plot = p_mapping,
    filename = "mapping-rates",
    output_dir = qc_dir,
    plot_type = "barplot"
  )
  
  # Create total reads plot
  p_reads <- ggplot(mapping_stats, aes(x = reorder(Sample, Total_Reads), y = Total_Reads / 1e6)) +
    geom_col(fill = "darkgreen", alpha = 0.7) +
    labs(
      title = "Total Reads by Sample",
      x = "Sample", 
      y = "Total Reads (Millions)"
    ) +
    get_standard_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save total reads plot
  save_standard_plot(
    plot = p_reads,
    filename = "total-reads",
    output_dir = qc_dir,
    plot_type = "barplot"
  )
  
  # Generate summary statistics
  summary_stats <- mapping_stats %>%
    summarise(
      n_samples = n(),
      mean_mapping_rate = mean(Mapping_Rate_Numeric, na.rm = TRUE),
      median_mapping_rate = median(Mapping_Rate_Numeric, na.rm = TRUE),
      min_mapping_rate = min(Mapping_Rate_Numeric, na.rm = TRUE),
      max_mapping_rate = max(Mapping_Rate_Numeric, na.rm = TRUE),
      mean_total_reads = mean(Total_Reads, na.rm = TRUE),
      median_total_reads = median(Total_Reads, na.rm = TRUE),
      samples_below_80pct = sum(Mapping_Rate_Numeric < 80, na.rm = TRUE)
    )
  
  # Save summary statistics
  write_csv(summary_stats, file.path(qc_dir, "qc_summary_statistics.csv"))
  
  cat("✓ QC plots generated:\n")
  cat("  - Mapping rates:", file.path(qc_dir, "mapping-rates.pdf"), "\n")
  cat("  - Total reads:", file.path(qc_dir, "total-reads.pdf"), "\n")
  cat("  - Summary stats:", file.path(qc_dir, "qc_summary_statistics.csv"), "\n")
}

# Check if count data exists for additional QC
count_file <- "data/counts/featurecounts_gene_name.tsv.gz"

if (file.exists(count_file)) {
  cat("\nGenerating count-based QC metrics...\n")
  
  # Load count data
  count_data <- read_tsv(count_file, show_col_types = FALSE)
  count_matrix <- as.matrix(select(count_data, -c(Geneid, Length)))
  
  # Calculate library sizes
  library_sizes <- colSums(count_matrix)
  
  # Create library size plot
  library_df <- data.frame(
    Sample = names(library_sizes),
    Library_Size = library_sizes
  )
  
  p_library <- ggplot(library_df, aes(x = reorder(Sample, Library_Size), y = Library_Size / 1e6)) +
    geom_col(fill = "purple", alpha = 0.7) +
    labs(
      title = "Library Size by Sample",
      x = "Sample",
      y = "Library Size (Millions of Reads)"
    ) +
    get_standard_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save library size plot
  save_standard_plot(
    plot = p_library,
    filename = "library-sizes",
    output_dir = qc_dir,
    plot_type = "barplot"
  )
  
  # Calculate percentage of genes with zero counts
  zero_counts <- rowSums(count_matrix == 0) / ncol(count_matrix) * 100
  
  # Create zero counts distribution plot
  zero_df <- data.frame(Zero_Percentage = zero_counts)
  
  p_zeros <- ggplot(zero_df, aes(x = Zero_Percentage)) +
    geom_histogram(bins = 50, fill = "orange", alpha = 0.7) +
    labs(
      title = "Distribution of Zero Count Percentages",
      x = "Percentage of Samples with Zero Counts",
      y = "Number of Genes"
    ) +
    get_standard_theme()
  
  # Save zero counts plot
  save_standard_plot(
    plot = p_zeros,
    filename = "zero-counts-distribution",
    output_dir = qc_dir,
    plot_type = "standard"
  )
  
  cat("✓ Count-based QC plots generated:\n")
  cat("  - Library sizes:", file.path(qc_dir, "library-sizes.pdf"), "\n")
  cat("  - Zero counts:", file.path(qc_dir, "zero-counts-distribution.pdf"), "\n")
}

cat("\n✓ Quality control analysis completed successfully!\n")