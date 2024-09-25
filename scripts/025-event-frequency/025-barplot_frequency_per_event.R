library(tidyverse)
library(extrafont)

df <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

df_filter <- df %>% filter(fdr < 0.05, abs(dpsi) > 0.1)

###########################################################
# Calculate the percentage of each `KO_Symbol` for each `Event`
###########################################################
df_percentage <- df_filter %>%
    group_by(ko_symbol, event) %>%
    summarise(n = n()) %>%
    mutate(perc = n / sum(n) * 100)

###########################################################
# Calculate the occurrence count of each `KO_Symbol` (the number of AS events)
###########################################################
ko_symbol_counts <- df_filter %>%
    group_by(ko_symbol) %>%
    count(ko_symbol) %>%
    mutate(ko_symbol_counts = str_glue("{ko_symbol} ({n})")) %>%
    select(ko_symbol, ko_symbol_counts)

###########################################################
# Add the ko_symbol_counts to the df_percentage
###########################################################
df_plot <-
    df_percentage %>%
    left_join(ko_symbol_counts, by = c("ko_symbol"))

###########################################################
# Order of the ko_symbol_counts and event
###########################################################

df_plot$ko_symbol_counts <- factor(df_plot$ko_symbol_counts)
df_plot$event <- factor(df_plot$event)
levels(df_plot$ko_symbol_counts) <- rev(sort(levels(df_plot$ko_symbol_counts)))

# levels(df_plot$ko_symbol_counts)
# levels(df_plot$event)


###########################################################
# Plot
###########################################################

colors <- c("#44ED8B", "#FF2FC1","#3FAFFF","#FFE270","#FF604E")

g_barplot <-
    ggplot(df_plot, aes(x = perc, y = ko_symbol_counts, fill = event)) + 
    geom_bar(stat = "identity", color = "black") +
    labs(title="FDR < 0.05 and abs(ΔPSI) > 0.1", x = "Percentage", y = "KO (No. of AS)", fill = "Event") +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_minimal() +
    theme(
        text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

###########################################################
# Save the plot
###########################################################

dir.create("reports/figure/", showWarnings = FALSE)
ggsave("reports/figure/015-barplot.pdf", g_barplot, width = 8, height = 6, family = "Arial", device = cairo_pdf)
ggsave("reports/figure/015-barplot.jpg", g_barplot, width = 8, height = 6)
