# Select human/mouse complex proteins with GO terms from ComplexTab

library(tidyverse)
library(janitor)

dir.create("data/Fig5", showWarnings = FALSE)
dir.create("reports/Fig5", showWarnings = FALSE)

if (!file.exists("data/Fig5/10090.tsv")) {
    # df_hs <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv") %>% clean_names()
    df_mm <- read_tsv("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/10090.tsv") %>% clean_names()
    # write_tsv(df_hs, "data/Fig5/9606.tsv")
    write_tsv(df_mm, "data/Fig5/10090.tsv")
}
# df_hs <- read_tsv("data/Fig5/9606.tsv") %>% select(taxon_id = taxonomy_identifier, uniprot = identifiers_and_stoichiometry_of_molecules_in_complex, ac = number_complex_ac, name = recommended_name)
df_mm <- read_tsv("data/Fig5/10090.tsv") %>% select(taxon_id = taxonomy_identifier, uniprot = identifiers_and_stoichiometry_of_molecules_in_complex, ac = number_complex_ac, name = recommended_name)


###########################################################
# Convert UNIPROT ID to Gene Symbol
###########################################################

df_uniprot_id <- read_csv("data/Fig5/uniprot_gene_symbol_mouse_human.csv")

df_format <-
    df_mm %>%
    # カッコの中身だけを削除する
    mutate(uniprot = str_remove_all(uniprot, "\\([^\\)]*\\)")) %>%
    # "-"以下を削除する
    mutate(uniprot = str_remove(uniprot, "-.*")) %>%
    # | で区切られたところを展開する
    separate_longer_delim(uniprot, delim = "|")

# UNIPROT IDをGene Symbolに変換する

inner_join(df_format, df_uniprot_id, by = c("taxon_id", "uniprot")) %>%
    select(symbol, name) %>%
    add_count(name) %>%
    write_tsv("data/Fig5/complextab_mouse.tsv")

