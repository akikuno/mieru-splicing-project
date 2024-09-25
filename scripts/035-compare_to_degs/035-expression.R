library(tidyverse)
library(extrafont)

df_deg <- read_csv("data/degs/ko_vs_mieru.csv")

df_counts <- df_deg %>%
    group_by(ko_symbol) %>%
    summarize(n = n()) %>%
    mutate(ko_symbol_counts = str_glue("{ko_symbol} ({n})")) %>%
    select(ko_symbol, ko_symbol_counts)

df_plot <- df_deg %>% left_join(df_counts, by = c("ko_symbol"))

g_violin <-
    ggplot(df_plot, aes(x = ko_symbol_counts, y = log2FoldChange)) +
    geom_violin() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14))  +
    labs(title="FDR < 0.05", x = "KO (No. of DEGs)", y = "log2 FC")

###########################################################
# Save the plot
###########################################################

width <- 18
height <- 10

dir.create("reports/figure/", showWarnings = FALSE)
ggsave("reports/figure/045-violin_expression.pdf", g_violin, width = width, height = height, family = "Arial", device = cairo_pdf)
ggsave("reports/figure/045-violin_expression.jpg", g_violin, width = width, height = height)
