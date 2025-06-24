# Select mouse complex proteins with GO terms from CORUM
# https://mips.helmholtz-muenchen.de/corum/download
# Corum 5.0 release (2024-09-03)

library(tidyverse)
library(janitor)

dir.create("data/Fig5", showWarnings = FALSE)

df_complex <- read_tsv("data/Fig5/corum_allComplexes.txt") %>% clean_names()

# df_complex %>% count(organism)

df_complex_go_symbol <-
    df_complex %>%
    filter(organism == "Mouse" | organism == "Human") %>%
    select(subunits_gene_name, functions_go_id, organism) %>%
    rename(symbol = subunits_gene_name, go = functions_go_id) %>%
    separate_longer_delim(symbol, delim = ";") %>%
    separate_longer_delim(go, delim = ";") %>%
    mutate(symbol = toupper(symbol), organism = tolower(organism)) %>%
    distinct() %>%
    select(go, symbol, organism) %>%
    arrange(go, symbol, organism)

df_complex_go_symbol %>% filter(go == "GO:0003723") # RNA binding
df_complex_go_symbol %>% write_csv("data/Fig5/corum_go_symbol_organism.csv")


