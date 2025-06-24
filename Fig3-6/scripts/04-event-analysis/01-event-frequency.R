#!/usr/bin/env Rscript
# Analysis of alternative splicing event frequency across RBP knockouts
# Generates barplot showing percentage distribution of event types per KO

# Load required libraries and utilities
source("scripts/utils/data-loaders.R")
source("scripts/utils/plot-themes.R")
library(yaml)

# Load configuration
config <- read_yaml("scripts/config/parameters.yaml")

# Load and filter rMATS data
cat("Loading rMATS data...\n")
df_rmats <- load_rmats_data(
  fdr_threshold = config$thresholds$fdr,
  dpsi_threshold = config$thresholds$dpsi
)

cat("Filtered data:", nrow(df_rmats), "significant events\n")

# Calculate percentage of each event type for each KO
df_percentage <- df_rmats %>%
  group_by(ko_symbol, event) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ko_symbol) %>%
  mutate(perc = n / sum(n) * 100) %>%
  ungroup()

# Calculate total event counts for labeling
ko_symbol_counts <- df_rmats %>%
  group_by(ko_symbol) %>%
  count(ko_symbol) %>%
  mutate(ko_symbol_counts = str_glue("{ko_symbol} ({n})")) %>%
  select(ko_symbol, ko_symbol_counts)

# Prepare plotting data
df_plot <- df_percentage %>%
  left_join(ko_symbol_counts, by = "ko_symbol")

# Set factor levels for ordering
df_plot$ko_symbol_counts <- factor(df_plot$ko_symbol_counts)
df_plot$event <- factor(df_plot$event)
levels(df_plot$ko_symbol_counts) <- rev(sort(levels(df_plot$ko_symbol_counts)))

# Create barplot
colors <- get_event_colors()

g_barplot <- ggplot(df_plot, aes(x = perc, y = ko_symbol_counts, fill = event)) + 
  geom_bar(stat = "identity", color = "black") +
  labs(
    title = str_glue("FDR < {config$thresholds$fdr} and abs(ΔPSI) > {config$thresholds$dpsi}"),
    x = "Percentage", 
    y = "KO (No. of AS)", 
    fill = "Event"
  ) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(reverse = TRUE)) +
  get_standard_theme()

# Save plot
save_standard_plot(
  plot = g_barplot,
  filename = "01-event-frequency-barplot",
  output_dir = config$paths$reports$fig3,
  plot_type = "barplot",
  formats = config$plotting$formats
)

cat("Event frequency analysis completed successfully\n")