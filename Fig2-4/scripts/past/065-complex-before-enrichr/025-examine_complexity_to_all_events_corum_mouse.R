library(tidyverse)
library(janitor)

if (!file.exists("data/Fig5/mgi_protein_coding_symbols.txt")) {
    url <- "https://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt"
    df_mgi <- read_tsv(url) %>% clean_names()
    # タンパク質コード遺伝子を抽出する
    df_mgi %>%
        select(marker_symbol, ref_seq_protein_i_ds) %>%
        rename(symbol = marker_symbol, id = ref_seq_protein_i_ds) %>%
        separate_longer_delim(id, delim = "|") %>%
        filter(str_detect(id, "^NP_")) %>%
        select(symbol) %>%
        distinct() %>%
        write_csv("data/Fig5/mgi_protein_coding_symbols.txt")
}
df_mgi_genes <- read_tsv("data/Fig5/mgi_protein_coding_symbols.txt") %>% mutate(symbol = toupper(symbol))
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")
df_corum <- read_csv("data/Fig5/corum_go_symbol_organism.csv")


df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(target_symbol = toupper(target_symbol)) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()
events <- df_all$event %>% unique()

###############################################################################
# マウスのみの複合体
###############################################################################

complex_genes <- df_corum %>% filter(organism == "mouse") %>% pull(symbol) %>% unique()
# input_ko_symbol <- ko_symbols[1]
# input_event <- events[1]
results_fisher <- tibble()
for (input_ko_symbol in ko_symbols){
    for (input_event in events) {
        spliced_genes <- df_spliced_genes %>% filter(ko_symbol == input_ko_symbol, event == input_event) %>% pull(target_symbol)
        non_spliced_genes <- df_mgi_genes %>% filter(!symbol %in% spliced_genes) %>% pull(symbol)

        print(c(input_ko_symbol, input_event, length(spliced_genes), length(non_spliced_genes)))

        overlap_spliced_complex <- spliced_genes %in% complex_genes
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

results_fisher %>% write_csv("reports/Fig5/fisher_corum_mouse.csv")

results_fisher %>% count(significance)

results_fisher %>% count(event, significance)
