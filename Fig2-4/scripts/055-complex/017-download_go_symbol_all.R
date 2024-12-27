library(tidyverse)

df_all_symbol_go <- read_tsv("data/Fig5/mgi_symbol_go.txt", col_names = c("symbol", "go")) %>%
    mutate(symbol = toupper(symbol))

df_complex_symbol_go <- read_csv("reports/Fig5/complex_go_symbol.csv")
