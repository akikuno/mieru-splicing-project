library(tidyverse)
library(extrafont)
library(patchwork)
library(ggVennDiagram)

df <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv") %>% filter(fdr < 0.05, abs(dpsi) > 0.1)

df_deg <- read_csv("data/degs/ko_vs_mieru.csv")

ko_symbols <- df_deg$ko_symbol %>% unique()
events <- df$event %>% unique()

input_ko_symbol <- ko_symbols[1]
input_event <- events[1]

###########################################################
# Plot Venn Diagram
###########################################################

width <- 20
height <- 18

for (input_event in events) {
    venn_list <-
        map(ko_symbols, function(input_ko_symbol) {
            df_deg_ko_symbol <- df_deg %>% filter(ko_symbol == input_ko_symbol) %>% select(target_symbol)
            df_rmats_ko_symbol <- df %>% filter(ko_symbol == input_ko_symbol, event == input_event) %>% select(target_symbol)
            list_names <- c("DEG", input_event)
            gene_list <- list(df_deg_ko_symbol$target_symbol, df_rmats_ko_symbol$target_symbol)
            gene_list <- setNames(gene_list, list_names)
            ggVennDiagram(gene_list) + labs(title = input_ko_symbol) + scale_fill_gradient(low="#EEE",high = "#FF604E")
        })

    g_venn <- wrap_plots(venn_list, nrow = 3)

    # Save the plot
    ggsave(str_glue("reports/figure/037-venn-{input_event}.pdf"), g_venn, width = width, height = height, family = "Arial", device = cairo_pdf)
    ggsave(str_glue("reports/figure/037-venn-{input_event}.jpg"), g_venn, width = width, height = height)
}

###########################################################
# Extract intersect
###########################################################

genes_intersect <- tibble()

for (input_event in events) {
    for (input_ko_symbol in ko_symbols) {
        df_deg_ko_symbol <- df_deg %>% filter(ko_symbol == input_ko_symbol) %>% select(target_symbol)
        df_se_ko_symbol <- df %>% filter(ko_symbol == input_ko_symbol, event == input_event) %>% select(target_symbol)
        gene_intersect <- intersect(df_deg_ko_symbol$target_symbol, df_se_ko_symbol$target_symbol)
        results <- tibble(ko_symbol = input_ko_symbol, event = input_event, se_deg_overlapped = gene_intersect)
        genes_intersect <- bind_rows(genes_intersect, results)
    }
}

write_csv(genes_intersect, "reports/all_events_deg_overlapped.csv")

