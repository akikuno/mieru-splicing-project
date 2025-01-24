###########################################################
# エントロピーとdPSIの相関
###########################################################

library(tidyverse)
library(patchwork)
library(ggsignif)

df_isoforms <- read_tsv("data/Fig7/tpm_isoforms_all.tsv.gz") %>%
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
    genes_complex <- df_complex %>%
        filter(ko_symbol == input_ko_symbol) %>%
        select(genes)
    df_ko_isoforms <- df_isoforms %>%
        filter(str_detect(sample, input_ko_symbol)) %>%
        inner_join(genes_complex, by = "genes", relationship = "many-to-many")
    df_mieru_isoforms <- df_isoforms %>%
        filter(str_detect(sample, "MIERU")) %>%
        inner_join(genes_complex, by = "genes", relationship = "many-to-many")

    df_entropy <- tibble()

    df_ko_entropy <- df_ko_isoforms %>%
        select(sample, genes, tpm, group) %>%
        # すべてのサンプルにおいて、IsoformのTPMの総和が10以上の遺伝子のみを抽出
        group_by(sample, genes) %>%
        filter(sum(tpm) >= 10) %>%
        ungroup() %>%
        # グループごとに、エントロピーを計算
        group_by(group, genes) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        select(group, genes, entropy) %>%
        distinct()

    df_mieru_entropy <- df_mieru_isoforms %>%
        select(sample, genes, tpm, group) %>%
        # すべてのサンプルにおいて、IsoformのTPMの総和が10以上の遺伝子のみを抽出
        group_by(sample, genes) %>%
        filter(sum(tpm) >= 10) %>%
        ungroup() %>%
        # グループごとに、エントロピーを計算
        group_by(group, genes) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        select(group, genes, entropy) %>%
        distinct()


    df_entropy <- inner_join(df_ko_entropy, df_mieru_entropy, by = "genes", suffix = c("_ko", "_mieru"))
    cat(input_ko_symbol)
    print(t.test(df_entropy$entropy_ko, df_entropy$entropy_mieru, paired = TRUE)$p.value)

    df_entropy_longer <-
        df_entropy %>%
        pivot_longer(
            cols = starts_with("group"),
            names_to = "key_group",
            values_to = "group"
        ) %>%
        pivot_longer(
            cols = starts_with("entropy"),
            names_to = "key_entropy",
            values_to = "entropy"
        ) %>%
        filter((key_group == "group_ko" & key_entropy == "entropy_ko") | (key_group == "group_mieru" & key_entropy == "entropy_mieru")) %>%
        select(-key_group, -key_entropy)

    df_entropy_longer <- df_entropy_longer %>%
        mutate(group = factor(group, levels = c("MIERU", input_ko_symbol)))

    fill_values <- c("MIERU" = "white") %>%
        c(setNames("#AAA", input_ko_symbol))

    g_plot <- ggplot(df_entropy_longer, aes(x = group, y = entropy, fill = group)) +
        geom_violin() +
        geom_boxplot(width = 0.1, fill = "white") +
        geom_signif(
            comparisons = list(c("MIERU", input_ko_symbol)),
            test = "t.test", test.args = list(paired = TRUE), na.rm = FALSE, map_signif_level = TRUE, col = "black", step_increase = 0.1
        ) +
        scale_fill_manual(values = fill_values) +
        labs(x = "", y = "Isoform diversity (entropy)") +
        theme_bw()

    g_list[[input_ko_symbol]] <- g_plot
    df_t_test <- bind_rows(df_t_test, tibble(
        ko_symbol = input_ko_symbol,
        p_value = t.test(df_ko_entropy$entropy, df_mieru_entropy$entropy)$p.value,
        num_complex = nrow(genes_complex),
        average_entropy_ko = mean(df_ko_entropy$entropy),
        average_entropy_mieru = mean(df_mieru_entropy$entropy)
    ))
}

g_wrap <- wrap_plots(g_list, ncol = 5)

ggsave("reports/Fig7/isoform_diversity.png", g_wrap, width = 20, height = 10)
df_t_test %>%
    write_csv("reports/Fig7/isoform_diversity_t_test.csv")
