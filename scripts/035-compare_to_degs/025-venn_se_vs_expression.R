library(tidyverse)
library(extrafont)
library(ggVennDiagram)

df_se <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv") %>% filter(fdr < 0.05, abs(dpsi) > 0.1) %>% filter(event == "SE")

df_deg <- read_csv("data/degs/ko_vs_mieru.csv")

df_se_deg <- df_se %>% left_join(df_deg, by = c("ko_symbol", "target_symbol"))

ko_symbols <- df_se_deg$ko_symbol %>% unique()

input_ko_symbol <- ko_symbols[1]

venn_list <-
    map(ko_symbols, function(input_ko_symbol) {
        df_deg_ko_symbol <- df_deg %>% filter(ko_symbol == input_ko_symbol) %>% select(target_symbol)
        df_se_ko_symbol <- df_se %>% filter(ko_symbol == input_ko_symbol) %>% select(target_symbol)

        gene_list <- list(DEG = df_deg_ko_symbol$target_symbol, SE = df_se_ko_symbol$target_symbol)
        ggVennDiagram(gene_list) + labs(title = input_ko_symbol) + scale_fill_gradient(low="#EEE",high = "#FF604E")
    })

g_venn <- wrap_plots(venn_list, nrow = 3)

###########################################################
# Save the plot
###########################################################

width <- 18
height <- 18

dir.create("reports/figure/", showWarnings = FALSE)
ggsave("reports/figure/035-venn.pdf", g_venn, width = width, height = height, family = "Arial", device = cairo_pdf)
ggsave("reports/figure/035-venn.jpg", g_venn, width = width, height = height)
