# Select human/mouse complex proteins with GO terms from ComplexTab

library(tidyverse)
library(janitor)

dir.create("reports/Fig5", showWarnings = FALSE)

if (!file.exists("data/Fig5/10090.tsv")) {
    df_hs <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv") %>% clean_names()
    df_mm <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/10090.tsv") %>% clean_names()
    write_tsv(df_hs, "data/Fig5/9606.tsv")
    write_tsv(df_mm, "data/Fig5/10090.tsv")
}
df_hs <-
    read_tsv("data/Fig5/9606.tsv") %>%
    select(taxon_id = taxonomy_identifier, uniprot = identifiers_and_stoichiometry_of_molecules_in_complex, go = go_annotations, name = recommended_name)

df_mm <-
    read_tsv("data/Fig5/10090.tsv") %>%
    select(taxon_id = taxonomy_identifier, uniprot = identifiers_and_stoichiometry_of_molecules_in_complex, go = go_annotations, name = recommended_name)


###########################################################
# Convert UNIPROT ID to Gene Symbol
###########################################################

df_combine <- bind_rows(df_hs, df_mm)

df_format <-
    df_combine %>%
    # カッコの中身だけを削除する
    mutate(uniprot = str_remove_all(uniprot, "\\([^\\)]*\\)")) %>%
    mutate(go = str_remove_all(go, "\\([^\\)]*\\)")) %>%
    # "-"以下を削除する
    mutate(uniprot = str_remove(uniprot, "-.*")) %>%
    mutate(go = str_remove(go, "-.*")) %>%
    # | で区切られたところを展開する
    separate_longer_delim(uniprot, delim = "|") %>%
    separate_longer_delim(go, delim = "|")


# UNIPROT IDをGene Symbolに変換する
df_uniprot_id <- read_csv("data/Fig5/uniprot_gene_symbol_mouse_human.csv")
df_uniprot_id %>% filter(symbol == "CCNL1")

df_format %>% filter(str_detect(uniprot, "Q9UK58"))
df_converted <-
    inner_join(df_format, df_uniprot_id, by = c("taxon_id", "uniprot")) %>%
    select(!uniprot)

df_hs %>% filter(str_detect(uniprot, "Q9UK58"))
df_converted %>% filter(symbol == "CCNL1")

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

df_complex <-
    full_join(df_converted, df_homology, by = "symbol") %>%
    mutate(symbol = case_when(
        taxon_id == 9606 & !is.na(mouse) ~ mouse,
        taxon_id == 10090 ~ symbol
    )) %>%
    filter(!is.na(symbol)) %>%
    select(!c(mouse, taxon_id)) %>%
    distinct()


write_csv(df_complex, "data/Fig5/complextab_human_mouse.csv")
