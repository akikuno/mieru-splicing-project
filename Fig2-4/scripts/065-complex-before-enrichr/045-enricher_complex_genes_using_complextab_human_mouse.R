library(enrichR)
library(tidyverse)
library(janitor)
library(patchwork)

dir.create("reports/Fig6", showWarnings = FALSE)

# ヒトとマウスのホモロジー遺伝子を取得する
if (!file.exists("data/Fig5/mgi_protein_coding_symbols.txt")) {
    col_names <- c("mgi_id", "mouse_symbol", "marker_id", "hgnc", "human_symbol", "ncbi")
    url <- "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt"
    df_homology <- read_tsv(url, col_names = col_names) %>% clean_names()
    df_homology %>%
        mutate(mouse = toupper(mouse_symbol), human = human_symbol) %>%
        select(mouse, human) %>%
        write_tsv("data/Fig5/mgi_homology_symbols.txt")
}
df_homology <- read_tsv("data/Fig5/mgi_homology_symbols.txt")
df_mgi_genes <- read_tsv("data/Fig5/mgi_protein_coding_symbols.txt") %>% mutate(symbol = toupper(symbol))
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

# Complexのデータの、ヒトの遺伝子をマウスの遺伝子に変更する
df_complextab <- read_csv("data/Fig5/complextab_name_go_symbol_organism.csv")
df_complextab <- df_complextab %>%
    left_join(df_homology, by = c("symbol" = "human")) %>% # human symbolに対応するmouse symbolを結合
    mutate(symbol = ifelse(organism == "human" & !is.na(mouse), mouse, symbol)) %>% # humanの場合のみsymbolを変換
    select(-mouse)

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(target_symbol = toupper(target_symbol)) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()
events <- df_all$event %>% unique()


###############################################################################
# ヒトとマウスの複合体
###############################################################################

complex_genes <- df_complextab %>%
    pull(symbol) %>%
    unique()

# input_ko_symbol <- ko_symbols[1]
# input_event <- "SE"
results_genes <- tibble()
for (input_ko_symbol in ko_symbols) {
    for (input_event in events) {
        spliced_genes <- df_spliced_genes %>%
            filter(ko_symbol == input_ko_symbol, event == input_event) %>%
            pull(target_symbol)

        overlap_spliced_complex <- spliced_genes %in% complex_genes

        print(c(input_ko_symbol, input_event, sum(overlap_spliced_complex)))

        results_genes <- bind_rows(results_genes, tibble(
            ko_symbol = input_ko_symbol,
            event = input_event,
            genes = spliced_genes[overlap_spliced_complex]
        ))
    }
}

###############################################################################
# Enrichr
###############################################################################

##############################################
# Setup Enrichr
##############################################

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr")
}

##############################################
# Search GOs
##############################################

dbs <- listEnrichrDbs() %>% as_tibble()
dbs_library <- pull(dbs, libraryName)
# dbs_library[str_detect(dbs_library, "BioPlanet")]
# dbs_library[str_detect(dbs_library, "WikiPathway")]
dbs_library[str_detect(dbs_library, "^GO_")]

dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023", "GO_Cellular_Component_2023")

enrichr_pathways <- tibble()

for (symbol in ko_symbols) {
    genes <- results_genes %>%
        filter(ko_symbol == symbol, event == "SE") %>%
        pull(genes) %>%
        unique()
    print(c(symbol, length(genes)))

    enriched <- enrichr(genes, dbs)
    for (db in dbs) {
        enrichr_pathway <- enriched[[db]] %>%
            as_tibble() %>%
            select(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>%
            filter(Adjusted.P.value < 0.05) %>%
            arrange(Adjusted.P.value) %>%
            mutate(ko_symbol = symbol, event = "SE", database = db, gene_number = length(genes)) %>%
            clean_names() %>%
            distinct()
        enrichr_pathways <- bind_rows(enrichr_pathways, enrichr_pathway)
    }
}

genes <- complex_genes
enriched <- enrichr(genes, dbs)
for (db in dbs) {
    enrichr_pathway <- enriched[[db]] %>%
        as_tibble() %>%
        select(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>%
        filter(Adjusted.P.value < 0.05) %>%
        arrange(Adjusted.P.value) %>%
        mutate(ko_symbol = "complex_genes", event = "SE", database = db, gene_number = length(genes)) %>%
        clean_names() %>%
        distinct()
    enrichr_pathways <- bind_rows(enrichr_pathways, enrichr_pathway)
}


enrichr_pathways <- enrichr_pathways %>%
    select(ko_symbol, event, database, term, overlap, adjusted_p_value, odds_ratio, combined_score, genes, gene_number) %>%
    arrange(ko_symbol, event, adjusted_p_value)

write_csv(enrichr_pathways, "reports/Fig6/enrichr_se_complextab_human_mouse.csv")


#########################################
# Plot Adjusted P-value
#########################################

input_ko_symbol <- "Qk"
input_event <- "SE"
input_db <- "GO_Molecular_Function_2023"
g_list_enrichr <- list()


for (input_ko_symbol in c(ko_symbols)) {
    for (input_db in dbs) {
        complex_pathways <- enrichr_pathways %>% filter(ko_symbol == "complex_genes", event == "SE", database == input_db)

        g_plot <- enrichr_pathways %>%
            filter(ko_symbol == input_ko_symbol, event == input_event, database == input_db) %>%
            # Complexでも有意となったのパスウェイを除く
            anti_join(complex_pathways, by = "term") %>%
            mutate(p_value = -log10(adjusted_p_value)) %>%
            mutate(term = fct_reorder(term, desc(adjusted_p_value))) %>%
            # slice_max(p_value, n = 10) %>%
            ggplot(aes(x = p_value, y = term)) +
            geom_col(fill = "red", alpha = 0.5) +
            geom_text(aes(label = term, x = 0.2), # テキストをX=0から配置
                hjust = 0, # 左寄せ
                size = 4, # フォントサイズ
                vjust = 0.5
            ) +
            # 0.05の線を引く
            geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#333") +
            theme_bw() +
            labs(title = input_ko_symbol, subtitle = input_db, x = "-log10(Adjusted P-value)", y = NULL) +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()
            )
        list_name <- str_c(input_ko_symbol, input_db, sep = "-")
        g_list_enrichr[[list_name]] <- g_plot
    }
}

g_wrap_enrichr <- wrap_plots(g_list_enrichr, ncol = 3)
ggsave("reports/Fig6/enrichr_se_complextab_human_mouse.jpg", width = 30, height = 45)
ggsave("reports/Fig6/enrichr_se_complextab_human_mouse.pdf", width = 30, height = 50)
