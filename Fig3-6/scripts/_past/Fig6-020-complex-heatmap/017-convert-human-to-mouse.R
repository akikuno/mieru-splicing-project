library(tidyverse)
###########################################################
# Convert Human genes to Mouse genes
###########################################################

# ヒトとマウスのホモロジー遺伝子を取得する
if (!file.exists("data/Fig5/mgi_protein_coding_symbols.txt")) {
    col_names <- c("mgi_id", "mouse_symbol", "marker_id", "hgnc", "human_symbol", "ncbi")
    url <- "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt"
    df_homology <- read_tsv(url, col_names = col_names) %>% clean_names()
    df_homology %>%
        mutate(mouse = mouse_symbol, human = human_symbol) %>%
        select(mouse, human) %>%
        write_tsv("data/Fig5/mgi_homology_symbols.txt")
}
df_homology <- read_tsv("data/Fig5/mgi_homology_symbols.txt") %>% rename(symbol = human)

df_go <- read_csv("data/Fig6/go_human_mouse.csv")

df_converted <-
    full_join(df_go, df_homology, by = "symbol") %>%
    mutate(symbol = case_when(
        taxon_id == 9606 & !is.na(mouse) ~ mouse,
        taxon_id == 10090 ~ symbol
    )) %>%
    filter(!is.na(symbol)) %>%
    select(!c(mouse, taxon_id)) %>%
    distinct()

write_csv(df_converted, "data/Fig6/go_annotation_human_mouse.csv")
