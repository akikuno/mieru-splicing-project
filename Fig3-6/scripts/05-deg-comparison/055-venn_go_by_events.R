library(tidyverse)
library(extrafont)
library(patchwork)
library(ggVennDiagram)
library(enrichR)

df_rmats <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv") %>% filter(fdr < 0.05, abs(dpsi) > 0.1)

df_deg <- read_csv("data/degs/ko_vs_mieru.csv")


##############################################
# Setup Enrichr
##############################################

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr")
}

##############################################
# Search WikiPathway
##############################################

dbs <- listEnrichrDbs() %>% as_tibble()
dbs_library <- pull(dbs, libraryName)
dbs_library[str_detect(dbs_library, "^GO_")]

dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023")

###########################################################
# Extract intersect
###########################################################

ko_symbols <- df_deg$ko_symbol %>% unique()
events <- df_rmats$event %>% unique()

enrichr_pathways <- tibble()
for (input_ko_symbol in ko_symbols) {
    for (input_event in events) {
        deg <- df_deg %>% filter(ko_symbol == input_ko_symbol) %>% pull(target_symbol) %>% unique() %>% toupper()
        rmats <- df_rmats %>% filter(ko_symbol == input_ko_symbol, event == input_event) %>% pull(target_symbol) %>% unique() %>% toupper()

        gene_sets <- list(deg = deg, rmats = rmats)

        for (set_name in names(gene_sets)) {
            genes <- gene_sets[[set_name]]
            print(c(input_ko_symbol, input_event, set_name, length(genes)))
            enriched <- enrichr(genes, dbs)
            for (db in dbs) {
                enrichr_pathway <- enriched[[db]] %>%
                    as_tibble() %>%
                    select(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>%
                    filter(Adjusted.P.value < 0.05) %>%
                    arrange(Adjusted.P.value) %>%
                    mutate(
                        ko_symbol = input_ko_symbol,
                        event = input_event,
                        type = set_name,
                        gene_number = length(genes),
                        db = db
                    )
                cat(c(input_ko_symbol, input_event, set_name, length(genes), db, nrow(enrichr_pathway)), "\n")
                enrichr_pathways <- bind_rows(enrichr_pathways, enrichr_pathway)
            }
        }
    }
}

write_csv(enrichr_pathways, "reports/Fig4/055-go_rmats_deg.csv")

###########################################################
# Plot Venn Diagram
###########################################################

width <- 20
height <- 18

dir.create("reports/Fig4/venn-go-by-events", showWarnings = FALSE)

for (input_event in events) {
    venn_list <-
        map(ko_symbols, function(input_ko_symbol) {
            go_deg <- enrichr_pathways %>% filter(ko_symbol == input_ko_symbol, type == "deg", event == input_event) %>% pull(Term)
            go_rmats <- enrichr_pathways %>% filter(ko_symbol == input_ko_symbol, type == "rmats", event == input_event) %>% pull(Term)
            list_names <- c("DEG", input_event)
            go_list <- list(go_deg, go_rmats)
            go_list <- setNames(go_list, list_names)
            ggVennDiagram(go_list) + labs(title = input_ko_symbol) + scale_fill_gradient(low="#EEE",high = "#FF604E")
        })

    g_venn <- wrap_plots(venn_list, nrow = 3)

    # Save the plot
    ggsave(str_glue("reports/Fig4/venn-go-by-events/055-venn_go_rmats_deg_{input_event}.pdf"), g_venn, width = width, height = height, family = "Arial", device = cairo_pdf)
    ggsave(str_glue("reports/Fig4/venn-go-by-events/055-venn_go_rmats_deg_{input_event}.jpg"), g_venn, width = width, height = height)
}
