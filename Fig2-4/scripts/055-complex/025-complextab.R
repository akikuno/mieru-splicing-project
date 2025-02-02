library(tidyverse)

df_all_symbol_go <- read_tsv("data/Fig5/mgi_symbol_go.txt", col_names = c("symbol", "go")) %>%
    mutate(symbol = toupper(symbol))

df_complex_symbol_go <- read_csv("data/Fig5/complextab_name_go_symbol_organism.csv")

#########################

df_go <- read_csv("reports/Fig4/enrichr_go_pathways.csv")

data <- df_go %>%
    select(Term, from = ko_symbol) %>%
    mutate(to = str_remove(Term, " \\(GO:.*"))

top_n <- 10
data_filtered <-
    data %>%
    count(to) %>%
    arrange(desc(n)) %>%
    slice_head(n = top_n) %>%
    inner_join(data, by = "to") %>%
    mutate(value = 1) %>%
    distinct() %>%
    select(from, to, value)

data_go_symbol <-
    data_filtered %>%
    select(to) %>%
    inner_join(data, by = "to", relationship = "many-to-many") %>%
    select(Term) %>%
    distinct() %>%
    inner_join(df_go, by = "Term") %>%
    select(Term, Genes) %>%
    group_by(Term) %>%
    summarise(Genes = paste(Genes, collapse = ", ")) %>%
    separate_longer_delim(Genes, delim = ";") %>%
    separate_longer_delim(Genes, delim = ", ") %>%
    distinct() %>%
    mutate(go = stringr::str_extract(Term, "\\(GO:\\d+\\)") %>% stringr::str_remove_all("[()]")) %>%
    select(Term, go, symbol = Genes)

results_fisher <- tibble()
# input_go <- "GO:0003723" # RNA binding
go_list <- data_go_symbol %>% pull(go) %>% unique()
for (input_go in go_list) {
    term <- data_go_symbol %>% filter(go == input_go) %>% pull(Term) %>% unique()
    data_symbols <- data_go_symbol %>% filter(go == input_go) %>% pull(symbol)
    all_symbols <- df_all_symbol_go %>% filter(go == input_go) %>% pull(symbol)
    complex_symbols <- df_complex_symbol_go %>% filter(go == input_go) %>% pull(symbol) %>% unique()

    overlap_data_complex <- data_symbols %in% complex_symbols
    overlap_all_complex <- all_symbols %in% complex_symbols

    a <- sum(overlap_data_complex)
    b <- sum(!overlap_data_complex)
    c <-sum (overlap_all_complex)
    d <- sum(!overlap_all_complex)
    vx <- matrix(c(a,b,c,d),nrow=2,byrow=T)
    result <- fisher.test(vx)
    sig <- ifelse(result$p.value < 0.05, "YES", "NO")
    results_fisher <- bind_rows(results_fisher, tibble(term = term, significance = sig, p_value = result$p.value, odds_ratio = result$estimate, se_forming_complex = a, se_not_forming_complex = b, all_go_genes_forming_complex = c, all_go_genes_not_forming_complex = d))
}

write_csv(results_fisher, "reports/Fig5/fisher_complextab.csv")

