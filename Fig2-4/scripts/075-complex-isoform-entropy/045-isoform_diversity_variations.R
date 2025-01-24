###########################################################
# Q 群間でのエントロピーのばらつきは大きいか？
# → デンドログラム、ヒートマップで示す
###########################################################

library(tidyverse)
library(patchwork)
library(ggsignif)

df_isoforms <- read_tsv("data/Fig7/tpm_isoforms_all.tsv.gz") %>%
    mutate(genes = toupper(gene_symbol)) %>%
    mutate(group = str_remove(sample, "_.*$")) %>%
    # すべてのサンプルにおいて、IsoformのTPMの総和が10以上の遺伝子のみを抽出
    group_by(sample, genes) %>%
    filter(sum(tpm) >= 10) %>%
    ungroup()

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

    df_ko_entropy <- df_ko_isoforms %>%
        select(sample, genes, tpm, group) %>%
        # サンプルごとに、エントロピーを計算
        group_by(sample, genes) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        select(!tpm) %>%
        distinct()

    df_mieru_entropy <- df_mieru_isoforms %>%
        select(sample, genes, tpm, group) %>%
        # サンプルごとに、エントロピーを計算
        group_by(sample, genes) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        select(!tpm) %>%
        distinct()

    df_entropy <- bind_rows(df_ko_entropy, df_mieru_entropy)

    pca_data <- df_entropy %>%
        select(sample, genes, entropy) %>%
        pivot_wider(names_from = genes, values_from = entropy, values_fill = 0) %>%
        column_to_rownames(var = "sample")

    pca_data <- pca_data[, apply(pca_data, 2, var) != 0]
    pca_result <- prcomp(pca_data, scale. = TRUE)

    # PCA 結果をデータフレーム化
    pca_df <- as_tibble(pca_result$x)
    pca_df$sample <- rownames(pca_data)

    # Group 情報を追加
    pca_df <- df_entropy %>%
        select(sample, group) %>%
        distinct() %>%
        left_join(pca_df, by = "sample")

    # プロット作成
    color_pallete <- c("MIERU" = "#555") %>%
        c(setNames("#E65C5C", input_ko_symbol))

    g_list[[input_ko_symbol]] <-
        ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 5) +
        scale_color_manual(values = color_pallete) +
        theme_bw() +
        labs(
            title = str_glue("{input_ko_symbol} PCA of Entropy"),
            x = "PC1",
            y = "PC2",
            color = "Group"
        )
    print(input_ko_symbol)
}

g_wrap <- wrap_plots(g_list, ncol = 5)

ggsave("reports/Fig7/isoform_diversity_pca.png", g_wrap, width = 22, height = 10)
