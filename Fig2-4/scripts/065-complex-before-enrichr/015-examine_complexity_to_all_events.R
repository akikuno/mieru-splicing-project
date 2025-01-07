library(tidyverse)
library(janitor)

df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

if (!file.exists("data/Fig5/mgi_all_symbols.txt")) {
    url <- "https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt"
    read_tsv(url) %>%
        clean_names() %>%
        select(x3_marker_symbol) %>%
        mutate(x3_marker_symbol = toupper(x3_marker_symbol)) %>%
        rename(symbol = x3_marker_symbol) %>%
        distinct() %>%
        write_csv("data/Fig5/mgi_all_symbols.txt")
}
df_mgi_genes <- read_tsv("data/Fig5/mgi_all_symbols.txt")

df_signif_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(target_symbol = toupper(target_symbol)) %>%
    distinct()

df_nonsignif_genes <-
    full_join(df_mgi_genes, df_signif_genes, by = c("symbol" = "target_symbol")) %>%
    filter(is.na(event)) %>%
    select(symbol) %>%
    distinct()

df_complex_genes <- read_csv("reports/Fig5/complextab_go_symbol.csv") %>% select(symbol) %>% distinct()

ko_symbols <- df$ko_symbol %>% unique()
events <- df$event %>% unique()

input_ko_symbol <- ko_symbols[1]
input_event <- events[1]

non_spliced_genes <- df_nonsignif_genes %>% pull(symbol)
complex_genes <- df_complex_genes %>% pull(symbol)

results_fisher <- tibble()
for (input_ko_symbol in ko_symbols){
    for (input_event in events) {
        signif_genes <- df_signif_genes %>% filter(ko_symbol == input_ko_symbol, event == input_event) %>% pull(target_symbol)
        print(c(input_ko_symbol, input_event, length(signif_genes), length(complex_genes), length(non_spliced_genes)))

        overlap_spliced_complex <- signif_genes %in% complex_genes
        overlap_non_spliced_complex <- non_spliced_genes %in% complex_genes

        a <- sum(overlap_spliced_complex)
        b <- sum(!overlap_spliced_complex)
        c <- sum (overlap_non_spliced_complex)
        d <- sum(!overlap_non_spliced_complex)
        vx <- matrix(c(a,b,c,d),nrow=2,byrow=T)
        result <- fisher.test(vx)
        sig <- ifelse(result$p.value < 0.05, "YES", "NO")
        results_fisher <- bind_rows(results_fisher, tibble(
            ko_symbol = input_ko_symbol,
            event = input_event,
            significance = sig,
            p_value = result$p.value,
            odds_ratio = result$estimate,
            spliced_genes_forming_complex = a,
            spliced_genes_not_forming_complex = b,
            non_spliced_genes_forming_complex = c,
            non_spliced_genes_not_forming_complex = d)
            )
    }
}

results_fisher %>% filter(significance == "YES") %>% as.data.frame()

