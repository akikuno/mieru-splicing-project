library(enrichR)
library(tidyverse)

df_rmats <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv") %>% filter(fdr < 0.05, abs(dpsi) > 0.1)

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
enrichr_pathways <- tibble()
ko_symbols <- unique(pull(df_rmats, ko_symbol))

for (symbol in ko_symbols) {
    genes <- df_rmats %>% filter(ko_symbol == symbol) %>% pull(target_symbol) %>% unique()
    print(c(symbol, length(genes)))

    enriched <- enrichr(genes, dbs)
    for (db in dbs) {
        enrichr_pathway <- enriched[[db]] %>%
            as_tibble() %>%
            select(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>%
            filter(Adjusted.P.value < 0.05) %>%
            arrange(Adjusted.P.value) %>%
            mutate(ko_symbol = symbol, gene_number = length(genes), db = db)
        enrichr_pathways <- bind_rows(enrichr_pathways, enrichr_pathway)
    }
}

dir.create("reports/Fig4/", showWarnings = FALSE)
write_csv(enrichr_pathways, "reports/Fig4/enrichr.csv")
