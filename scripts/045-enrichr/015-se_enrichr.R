library(enrichR)
library(tidyverse)

df_se <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv") %>% filter(fdr < 0.05, abs(dpsi) > 0.1) %>% filter(event == "SE")

df_se <- df_se %>% mutate(dpsi_sign = ifelse(dpsi > 0, "skipping", "including"))

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
dbs_library[str_detect(dbs_library, "BioPlanet")]
dbs_library[str_detect(dbs_library, "WikiPathway")]
dbs_library[str_detect(dbs_library, "GO")]

dbs <- c("WikiPathway_2023_Human", "Reactome_2022", "KEGG_2021_Human", "BioPlanet_2019", "GO_Molecular_Function_2023", "GO_Biological_Process_2023", "GO_Cellular_Component_2023")

if (!file.exists("reports/enrichr_go_pathways.csv")) {
    enrichr_pathways <- tibble()
    ko_symbols <- unique(pull(df_se, ko_symbol))

    for (symbol in ko_symbols) {
        for (sign in c("skipping", "including")) {
            genes <- df_se %>% filter(ko_symbol == symbol, dpsi_sign == sign) %>% pull(target_symbol) %>% unique()
            print(c(symbol, sign, length(genes)))

            enriched <- enrichr(genes, dbs)
            for (db in dbs) {
                enrichr_pathway <- enriched[[db]] %>%
                    as_tibble() %>%
                    select(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>%
                    filter(Adjusted.P.value < 0.05) %>%
                    arrange(Adjusted.P.value) %>%
                    mutate(ko_symbol = symbol, se_event = sign, gene_number = length(genes), db = db)
                enrichr_pathways <- bind_rows(enrichr_pathways, enrichr_pathway)
            }
        }
    }

    write_csv(enrichr_pathways, "reports/enrichr_go_pathways.csv")
    }
