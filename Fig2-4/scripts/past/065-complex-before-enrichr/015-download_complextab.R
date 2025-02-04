# Select human/mouse complex proteins with GO terms from ComplexTab

library(tidyverse)
library(janitor)

dir.create("data/Fig5", showWarnings = FALSE)
dir.create("reports/Fig5", showWarnings = FALSE)

if (!file.exists("data/Fig5/10090.tsv")) {
    df_hs <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv") %>% clean_names()
    df_mm <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/10090.tsv") %>% clean_names()
    write_tsv(df_hs, "data/Fig5/9606.tsv")
    write_tsv(df_mm, "data/Fig5/10090.tsv")
}
df_hs <- read_tsv("data/Fig5/9606.tsv") %>% select(taxon_id = taxonomy_identifier, uniprot = identifiers_and_stoichiometry_of_molecules_in_complex, ac = number_complex_ac, name = recommended_name)
df_mm <- read_tsv("data/Fig5/10090.tsv") %>% select(taxon_id = taxonomy_identifier, uniprot = identifiers_and_stoichiometry_of_molecules_in_complex, ac = number_complex_ac, name = recommended_name)


###########################################################
# Convert UNIPROT ID to Gene Symbol
###########################################################

df_uniprot_id <- read_csv("data/Fig5/uniprot_gene_symbol_mouse_human.csv")

df_combine <- bind_rows(df_hs, df_mm)
df_format <-
    df_combine %>%
    # カッコの中身だけを削除する
    mutate(uniprot = str_remove_all(uniprot, "\\([^\\)]*\\)")) %>%
    # "-"以下を削除する
    mutate(uniprot = str_remove(uniprot, "-.*")) %>%
    # | で区切られたところを展開する
    separate_longer_delim(uniprot, delim = "|")

# UNIPROT IDをGene Symbolに変換する

df_format <- inner_join(df_format, df_uniprot_id, by = c("taxon_id", "uniprot")) %>%
    select(!uniprot)


###########################################################
# Convert Human genes to Mouse genes
###########################################################

# TODO 面倒くさいので後回し

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

# symbol, nameをヒトに合わせる

df_format <- df_format %>% mutate(name_upper = toupper(name))

df_unique_name <-
    df_format %>%
    filter(taxon_id == 9606) %>%
    select(ac, name_upper) %>%
    distinct()

df_complex <- full_join(df_format, df_homology, by = "symbol") %>%
    mutate(symbol = case_when(
        taxon_id == 9606 & !is.na(mouse) ~ toupper(mouse),
        taxon_id == 10090 ~ toupper(symbol)
    )) %>%
    filter(!is.na(symbol)) %>%
    select(!c(mouse, taxon_id)) %>%
    distinct()


df_complex %>% filter(symbol == "Cdk11")
filter(!is.na(taxon_id)) %>%
    df_symbol() %>%
    filter(symbol == "Atf1")
# keys <- df_format$id %>% unique()

# map_hs <- select(
#     x = org.Hs.eg.db,
#     keys = keys,
#     keytype = "UNIPROT",
#     columns = c("UNIPROT", "SYMBOL")
# )

# map_mm <- select(
#     x = org.Mm.eg.db,
#     keys = keys,
#     keytype = "UNIPROT",
#     columns = c("UNIPROT", "SYMBOL")
# )

# map_all <- bind_rows(map_hs, map_mm) %>%
#     as_tibble() %>%
#     # Remove NA
#     filter(!is.na(SYMBOL)) %>%
#     dplyr::select(id = UNIPROT, symbol = SYMBOL) %>%
#     mutate(symbol = toupper(symbol))

# df_complex_go_symbol <- df_format %>%
#     inner_join(map_all, by = c("id" = "id"), relationship = "many-to-many") %>%
#     dplyr::select(name, go, symbol, organism) %>%
#     distinct() %>%
#     arrange(symbol, organism)


# df_complex_go_symbol %>% filter(go == "GO:0003723")
# df_complex_go_symbol %>% write_csv("data/Fig5/complextab_name_go_symbol_organism.csv")
