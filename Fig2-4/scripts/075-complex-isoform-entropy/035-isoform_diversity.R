library(tidyverse)
library(patchwork)
library(ggsignif)

df_isoforms <- read_tsv("data/Fig7/tpm_isoforms_all.tsv") %>%
    mutate(genes = toupper(gene_symbol)) %>%
    mutate(group = str_remove(sample, "_.*$"))

df_homology <- read_tsv("data/Fig5/mgi_homology_symbols.txt")
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

# Complexのデータの、ヒトの遺伝子をマウスの遺伝子に変更する
df_complextab <- read_csv("data/Fig5/complextab_name_go_symbol_organism.csv")
df_complextab <- df_complextab %>%
    left_join(df_homology, by = c("symbol" = "human")) %>% # human symbolに対応するmouse symbolを結合
    mutate(symbol = ifelse(organism == "human" & !is.na(mouse), mouse, symbol)) %>% # humanの場合のみsymbolを変換
    select(-mouse)

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(target_symbol = toupper(target_symbol)) %>%
    distinct()

ko_symbols <- df_isoforms %>%
    filter(group != "MIERU") %>%
    pull(group) %>%
    unique()
events <- df_all$event %>% unique()


###############################################################################
# ヒトとマウスの複合体
###############################################################################

complex_genes <- df_complextab %>%
    pull(symbol) %>%
    unique()

# input_ko_symbol <- ko_symbols[1]
# input_event <- "SE"
df_complex <- tibble()
for (input_ko_symbol in ko_symbols) {
    for (input_event in events) {
        spliced_genes <- df_spliced_genes %>%
            filter(ko_symbol == input_ko_symbol, event == input_event) %>%
            pull(target_symbol)

        overlap_spliced_complex <- spliced_genes %in% complex_genes
        df_complex <- bind_rows(df_complex, tibble(
            ko_symbol = input_ko_symbol,
            event = input_event,
            genes = spliced_genes[overlap_spliced_complex]
        ))
    }
}

###########################################################
# エントロピーによるisoformの多様性を検定
###########################################################

calculate_entropy <- function(values) {
    total_sum <- sum(values)
    proportions <- values / total_sum
    -sum(proportions * log(proportions), na.rm = TRUE)
}

g_list <- list()
df_t_test <- tibble()
input_ko_symbol <- ko_symbols[1]

for (input_ko_symbol in ko_symbols) {
    genes_se_complex <- df_complex %>%
        filter(ko_symbol == input_ko_symbol) %>%
        select(genes)
    df_ko_isoforms <- df_isoforms %>%
        filter(str_detect(sample, input_ko_symbol)) %>%
        inner_join(genes_se_complex, by = "genes", relationship = "many-to-many")
    df_mieru_isoforms <- df_isoforms %>%
        filter(str_detect(sample, "MIERU")) %>%
        inner_join(genes_se_complex, by = "genes", relationship = "many-to-many")

    df_entropy <- tibble()

    df_ko_entropy <- df_ko_isoforms %>%
        select(sample, genes, tpm, group) %>%
        group_by(sample, genes) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        select(sample, group, entropy) %>%
        distinct()


    df_mieru_entropy <- df_mieru_isoforms %>%
        select(sample, genes, tpm, group) %>%
        group_by(sample, genes) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        select(sample, group, entropy) %>%
        distinct()


    df_entropy <- bind_rows(df_ko_entropy, df_mieru_entropy)

    df_entropy <- df_entropy %>%
        mutate(group = factor(group, levels = c("MIERU", input_ko_symbol)))

    fill_values <- c("MIERU" = "white") %>%
        c(setNames("#AAA", input_ko_symbol))

    g_plot <- ggplot(df_entropy, aes(x = group, y = entropy, fill = group)) +
        geom_violin() +
        geom_boxplot(width = 0.1, fill = "white") +
        geom_signif(
            comparisons = list(c("MIERU", input_ko_symbol)),
            test = "t.test", na.rm = FALSE, map_signif_level = TRUE, col = "black", step_increase = 0.1
        ) +
        scale_fill_manual(values = fill_values) +
        labs(x = "Sample", y = "Isoform diversity (entropy)") +
        theme_bw()

    g_list[[input_ko_symbol]] <- g_plot
    df_t_test <- bind_rows(df_t_test, tibble(ko_symbol = input_ko_symbol, p_value = t.test(df_ko_entropy$entropy, df_mieru_entropy$entropy)$p.value))
}

g_wrap <- wrap_plots(g_list, ncol = 4)

ggsave("reports/Fig7/isoform_diversity.png", g_wrap, width = 15, height = 5)
# df_ko_isoforms <- df_isoforms %>%
#     filter(!str_detect(sample, "MIERU")) %>%
#     inner_join(genes_se_complex, by = "genes", relationship = "many-to-many")
