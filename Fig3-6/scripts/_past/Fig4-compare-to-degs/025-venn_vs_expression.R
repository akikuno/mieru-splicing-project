library(tidyverse)
library(extrafont)
library(patchwork)
library(ggVennDiagram)

df_spliced <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv") %>% filter(fdr < 0.05, abs(dpsi) > 0.1)

df_deg <- read_csv("data/degs/ko_vs_mieru.csv")

ko_symbols <- df_deg$ko_symbol %>% unique()
input_ko_symbol <- ko_symbols[1]

###########################################################
# Plot Venn Diagram
###########################################################

g_venn <- list()
for (input_ko_symbol in ko_symbols) {
    df_deg_symbol <- df_deg %>%
        filter(ko_symbol == input_ko_symbol) %>%
        pull(target_symbol) %>%
        unique()
    df_spliced_symbol <- df_spliced %>%
        filter(ko_symbol == input_ko_symbol) %>%
        pull(target_symbol) %>%
        unique()
    gene_list <- list(DEG = df_deg_symbol, DSG = df_spliced_symbol)
    g <- ggVennDiagram(gene_list) + labs(title = input_ko_symbol) + scale_fill_gradient(low = "#EEE", high = "#FF604E") + theme(text = element_text(size = 20))
    g_venn[[input_ko_symbol]] <- g
}

g_venn <- wrap_plots(g_venn, nrow = 3)

###########################################################
# Save the plot
###########################################################

width <- 18
height <- 18

dir.create("reports/Fig4/", showWarnings = FALSE)
ggsave("reports/Fig4/025-venn_overlap_dsg_deg.pdf", g_venn, width = width, height = height, family = "Arial", device = cairo_pdf)
ggsave("reports/Fig4/025-venn_overlap_dsg_deg.jpg", g_venn, width = width, height = height)

###########################################################
# Extract intersect
###########################################################

genes_intersect <-
    map_dfr(ko_symbols, function(input_ko_symbol) {
        df_deg_ko_symbol <- df_deg %>%
            filter(ko_symbol == input_ko_symbol) %>%
            select(target_symbol)
        df_rmats_ko_symbol <- df_spliced %>%
            filter(ko_symbol == input_ko_symbol) %>%
            select(target_symbol)
        gene_intersect <- intersect(df_deg_ko_symbol$target_symbol, df_rmats_ko_symbol$target_symbol)
        tibble(ko_symbol = input_ko_symbol, se_deg_overlapped = gene_intersect)
    })

write_csv(genes_intersect, "reports/Fig4/025-overlap_dsg_deg.csv")
