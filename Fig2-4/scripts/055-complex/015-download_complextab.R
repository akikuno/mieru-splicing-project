# Select human/mouse complex proteins with GO terms from ComplexTab

library(tidyverse)
library(janitor)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

dir.create("data/Fig5", showWarnings = FALSE)
dir.create("reports/Fig5", showWarnings = FALSE)

if (!file.exists("data/Fig5/10090.tsv")) {
    df_hs <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv") %>% clean_names()
    df_mm <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/10090.tsv") %>% clean_names()
    write_tsv(df_hs, "data/Fig5/9606.tsv")
    write_tsv(df_mm, "data/Fig5/10090.tsv")
}
df_hs <- read_tsv("data/Fig5/9606.tsv") %>% mutate(organism = "human")
df_mm <- read_tsv("data/Fig5/10090.tsv") %>% mutate(organism = "mouse")

# df_hs %>% filter(recommended_name == "Cyclin L1-CDK11A(p110) complex") %>% select(accession = number_complex_ac, id = identifiers_and_stoichiometry_of_molecules_in_complex)
df_hs %>%
    dplyr::select(id = identifiers_and_stoichiometry_of_molecules_in_complex) %>%
    mutate(id = str_remove_all(id, "\\([^\\)]*\\)")) %>%
    separate_longer_delim(id, delim = "|") %>%
    distinct()
df_mm %>%
    dplyr::select(id = identifiers_and_stoichiometry_of_molecules_in_complex) %>%
    mutate(id = str_remove_all(id, "\\([^\\)]*\\)")) %>%
    separate_longer_delim(id, delim = "|") %>%
    distinct()

# df_mm %>% select(recommended_name, identifiers_and_stoichiometry_of_molecules_in_complex)

df_combine <- bind_rows(df_hs, df_mm)
# df_combine <- df_mm
df_format <- df_combine %>%
    dplyr::select(name = recommended_name, go = go_annotations, id = identifiers_and_stoichiometry_of_molecules_in_complex, organism) %>%
    # カッコの中身だけを削除する
    mutate(go = str_remove_all(go, "\\([^\\)]*\\)")) %>%
    mutate(id = str_remove_all(id, "\\([^\\)]*\\)")) %>%
    # | で区切られたところを展開する
    separate_longer_delim(go, delim = "|") %>%
    separate_longer_delim(id, delim = "|")

# UNIPROT IDをGene Symbolに変換する

keys <- df_format$id %>% unique()

map_hs <- select(
    x = org.Hs.eg.db,
    keys = keys,
    keytype = "UNIPROT",
    columns = c("UNIPROT", "SYMBOL")
)

map_mm <- select(
    x = org.Mm.eg.db,
    keys = keys,
    keytype = "UNIPROT",
    columns = c("UNIPROT", "SYMBOL")
)

map_all <- bind_rows(map_hs, map_mm) %>%
    as_tibble() %>%
    # Remove NA
    filter(!is.na(SYMBOL)) %>%
    dplyr::select(id = UNIPROT, symbol = SYMBOL) %>%
    mutate(symbol = toupper(symbol))

df_complex_go_symbol <- df_format %>%
    inner_join(map_all, by = c("id" = "id"), relationship = "many-to-many") %>%
    dplyr::select(name, go, symbol, organism) %>%
    distinct() %>%
    arrange(symbol, organism)


df_complex_go_symbol %>% filter(go == "GO:0003723")
df_complex_go_symbol %>% write_csv("data/Fig5/complextab_name_go_symbol_organism.csv")
