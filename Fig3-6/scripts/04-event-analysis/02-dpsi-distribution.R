#!/usr/bin/env Rscript
# Analysis of delta PSI distribution by event type
# Generates violin plots showing ΔPSI distribution for each event type

# Load required libraries and utilities
source("scripts/utils/data-loaders.R")
source("scripts/utils/plot-themes.R")
library(yaml)
library(patchwork)

# Load configuration
config <- read_yaml("scripts/config/parameters.yaml")

# Load and filter rMATS data
cat("Loading rMATS data...\n")
df_rmats <- load_rmats_data(
  fdr_threshold = config$thresholds$fdr,
  dpsi_threshold = config$thresholds$dpsi
)

# Calculate event counts for labeling
event_counts <- df_rmats %>%
  group_by(ko_symbol, event) %>%
  count(event) %>%
  ungroup() %>%
  mutate(event_counts = str_glue("{ko_symbol} ({n})")) %>%
  select(!n)

# Prepare plotting data
df_plot <- df_rmats %>%
  left_join(event_counts, by = c("ko_symbol", "event"))

# Get colors and order events
colors <- get_event_colors()
event_ordered <- rev(sort(unique(df_plot$event)))
len_ko_symbol <- length(unique(df_plot$ko_symbol))

# Create violin plots for each event type
violin_list <- map(event_ordered, function(target_event) {
  df_plot_event <- df_plot %>% filter(event == target_event)
  df_plot_event$event_counts <- factor(df_plot_event$event_counts)
  levels(df_plot_event$event_counts) <- rev(sort(levels(df_plot_event$event_counts)))
  
  ggplot(df_plot_event, aes(x = dpsi, y = event_counts, fill = event)) +
    geom_violin() +
    labs(title = target_event, x = "ΔPSI", y = "", fill = "Event") +
    scale_fill_manual(values = rep(colors[[target_event]], len_ko_symbol)) +
    get_standard_theme() +
    theme(
      legend.position = 'none',
      axis.text.y = element_text(size = 16),
      strip.text = element_text(size = 18)
    )
})

# Combine plots
g_violinplot <- wrap_plots(violin_list, nrow = 1)

# Save plot
save_standard_plot(
  plot = g_violinplot,
  filename = "02-dpsi-distribution-violin",
  output_dir = config$paths$reports$fig3,
  plot_type = "wide_violin",
  formats = config$plotting$formats
)

cat("ΔPSI distribution analysis completed successfully\n")