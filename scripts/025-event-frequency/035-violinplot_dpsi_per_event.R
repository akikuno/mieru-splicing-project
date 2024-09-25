library(tidyverse)
library(extrafont)
library(patchwork)

df <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

df_filter <- df %>% filter(fdr < 0.05, abs(dpsi) > 0.1)

###########################################################
# Calculate the occurrence count of each `event`
###########################################################
event_counts <- df_filter %>% 
    group_by(ko_symbol, event) %>%
    count(event) %>%
    ungroup() %>%
    mutate(event_counts = str_glue("{ko_symbol} ({n})")) %>%
    select(!n)

###########################################################
# Add the event_counts
###########################################################
df_plot <-
    df_filter %>%
    left_join(event_counts, by = c("ko_symbol", "event"))

###########################################################
# Plot
###########################################################

colors <- c("#44ED8B", "#FF2FC1","#3FAFFF","#FFE270","#FF604E")
names(colors) <- c("A3SS", "A5SS", "MXE", "RI", "SE")

target_event <- "SE"
len_ko_symbol <- length(unique(df_plot$ko_symbol))
df_plot_event <- df_plot %>% filter(event == target_event)
df_plot_event$event_counts <- factor(df_plot_event$event_counts)
levels(df_plot_event$event_counts) <- rev(sort(levels(df_plot_event$event_counts)))

event_ordered <- rev(sort(unique(df_plot$event)))
violin_list <-
    map(event_ordered, ~{
        df_plot_event <- df_plot %>% filter(event == .x)
        df_plot_event$event_counts <- factor(df_plot_event$event_counts)
        levels(df_plot_event$event_counts) <- rev(sort(levels(df_plot_event$event_counts)))

        ggplot(df_plot_event, aes(x = dpsi, y = event_counts, fill = event)) +
            geom_violin() +
            labs(title=.x, x = "ΔPSI", y = "", fill = "Event") +
            scale_fill_manual(values = rep(colors[[.x]], len_ko_symbol)) +
            theme_minimal() +
            theme(
                text = element_text(family = "Arial"),
                legend.position = 'none',
                axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 16),
                strip.text = element_text(size = 18))
    })

g_violinplot <- wrap_plots(violin_list, nrow = 1)

###########################################################
# Save the plot
###########################################################

dir.create("reports/figure/", showWarnings = FALSE)
ggsave("reports/figure/025-violinplot.pdf", g_violinplot, width = 24, height = 6, family = "Arial", device = cairo_pdf)
ggsave("reports/figure/025-violinplot.jpg", g_violinplot, width = 24, height = 6)
