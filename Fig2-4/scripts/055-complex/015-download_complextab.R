library(tidyverse)
library(janitor)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
dir.create("reports/Fig5", showWarnings = FALSE)

if (!file.exists("reports/Fig5/9606.rds")) {
    df_hs <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv") %>% clean_names()
    df_mm <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/10090.tsv") %>% clean_names()
    saveRDS(df_hs, "reports/Fig5/9606.rds")
    saveRDS(df_mm, "reports/Fig5/10090.rds")
} else {
    df_hs <- readRDS("reports/Fig5/9606.rds")
    df_mm <- readRDS("reports/Fig5/10090.rds")
}

df_combine <- bind_rows(df_hs, df_mm)
df_format <- df_combine %>%
    dplyr::select(go = go_annotations, id = identifiers_and_stoichiometry_of_molecules_in_complex) %>%
    # カッコの中身だけを削除する
    mutate(go = str_remove_all(go, "\\([^\\)]*\\)")) %>%
    mutate(id = str_remove_all(id, "\\([^\\)]*\\)")) %>%
    # | で区切られたところを展開する
    separate_longer_delim(go, delim = "|") %>%
    separate_longer_delim(id, delim = "|")

# UNIPROT IDをGene Symbolに変換する

keys <- df_format$id %>% unique()

map_hs <- select(x = org.Hs.eg.db,
    keys = keys,
    keytype = "UNIPROT",
    columns = c("UNIPROT", "SYMBOL"))

map_mm <- select(x = org.Mm.eg.db,
    keys = keys,
    keytype = "UNIPROT",
    columns = c("UNIPROT", "SYMBOL"))

map_all <- bind_rows(map_hs, map_mm) %>%
    as_tibble() %>%
    # Remove NA
    filter(!is.na(SYMBOL)) %>%
    dplyr::select(id = UNIPROT, symbol = SYMBOL) %>%
    mutate(symbol = toupper(symbol))

df_complex_go_symbol <- df_format %>%
    inner_join(map_all, by = c("id" = "id"), relationship = "many-to-many") %>%
    dplyr::select(go, symbol) %>%
    distinct() %>%
    arrange(go, symbol)


df_complex_go_symbol %>% write_csv("reports/Fig5/complextab_go_symbol.csv")
df_complex_go_symbol %>% filter(go == "GO:0003723")
