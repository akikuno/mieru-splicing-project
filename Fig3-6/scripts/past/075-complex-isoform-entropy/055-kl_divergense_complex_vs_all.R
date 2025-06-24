library(tidyverse)
library(patchwork)
library(ggsignif)

df_isoforms <- read_tsv("data/Fig7/tpm_isoforms_all.tsv.gz") %>%
    mutate(genes = toupper(gene_symbol)) %>%
    mutate(group = str_remove(sample, "_.*$")) %>%
    # すべてのサンプルにおいて、IsoformのTPMの総和が10以上の遺伝子のみを抽出
    group_by(sample, genes) %>%
    filter(sum(tpm) >= 10) %>%
    ungroup() %>%
    select(group, sample, genes, tpm)

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
# ASイベントを受けたタンパク質複合体と、すべての遺伝子群の多様性を、
# KLダイバージェンス（KO vs Control）を使って検定する
###########################################################

kl_divergence <- function(A, B) {
    epsilon <- 1e-10  # 小さな値を加える
    P <- (A + epsilon) / sum(A + epsilon)
    Q <- (B + epsilon) / sum(B + epsilon)
    sum(P * log(P / Q), na.rm = TRUE)
}

A = c(0, 0.2, 0.3, 0.4)
B = c(0.1, 0.2, 0.3, 0.4)
print(kl_divergence(A, B))
# ------------------------------------------------------------
# すべての遺伝子群におけるKLダイバージェンスを計算
# ------------------------------------------------------------

input_ko_symbol <- ko_symbols[1]

df_all_mieru_isoforms <- df_isoforms %>% filter(str_detect(sample, "MIERU"))
df_all_mieru_percentile <- df_all_mieru_isoforms %>%
        select(!sample) %>%
        group_by(group, genes) %>%
        # 0.25, 0.5, 0.75分位数を計算
        reframe(
            percentile = c(0.25, 0.5, 0.75),
            tpm = quantile(tpm, probs = c(0.25, 0.5, 0.75))
        ) %>% ungroup()

df_kl_divergence <- tibble()

for (input_ko_symbol in ko_symbols) {
    df_ko_isoforms <- df_isoforms %>% filter(str_detect(sample, input_ko_symbol))
    df_ko_percentile <- df_ko_isoforms %>%
        select(!sample) %>%
        group_by(group, genes) %>%
        # 0.25, 0.5, 0.75分位数を計算
        reframe(
            percentile = c(0.25, 0.5, 0.75),
            tpm = quantile(tpm, probs = c(0.25, 0.5, 0.75))
        ) %>% ungroup()

    df_entropy <- tibble()

    df_all_kl_divergence <- inner_join(df_ko_percentile, df_all_mieru_percentile, by = c("genes", "percentile"), suffix = c("_ko", "_mieru")) %>%
        select(genes, tpm_ko, tpm_mieru) %>%
        group_by(genes) %>%
        reframe(kl = kl_divergence(tpm_ko, tpm_mieru)) %>%
        ungroup() %>%
        mutate(
            ko_symbol = input_ko_symbol,
            filter = "all_genes"
        )


    df_kl_divergence <- bind_rows(df_kl_divergence, df_all_kl_divergence)
    print(input_ko_symbol)
    print(nrow(df_kl_divergence))
}


# ------------------------------------------------------------
# タンパク質複合体に絞ったKLダイバージェンスを計算
# ------------------------------------------------------------

input_ko_symbol <- ko_symbols[1]

for (input_ko_symbol in ko_symbols) {
    genes_complex <- df_complex %>%
        filter(ko_symbol == input_ko_symbol) %>%
        select(genes)

    df_ko_percentile <- df_isoforms %>%
        filter(str_detect(sample, input_ko_symbol)) %>%
        inner_join(genes_complex, by = "genes", relationship = "many-to-many") %>%
        select(!sample) %>%
        group_by(group, genes) %>%
        # 0.25, 0.5, 0.75分位数を計算
        reframe(
            percentile = c(0.25, 0.5, 0.75),
            tpm = quantile(tpm, probs = c(0.25, 0.5, 0.75))
        ) %>% ungroup()

    df_mieru_percentile <- df_isoforms %>%
        filter(str_detect(sample, "MIERU")) %>%
        inner_join(genes_complex, by = "genes", relationship = "many-to-many") %>%
        select(!sample) %>%
        group_by(group, genes) %>%
        # 0.25, 0.5, 0.75分位数を計算
        reframe(
            percentile = c(0.25, 0.5, 0.75),
            tpm = quantile(tpm, probs = c(0.25, 0.5, 0.75))
        ) %>% ungroup()

    df_complex_kl_divergence <- inner_join(df_ko_percentile, df_mieru_percentile, by = c("genes", "percentile"), suffix = c("_ko", "_mieru")) %>%
            select(genes, tpm_ko, tpm_mieru) %>%
            group_by(genes) %>%
            reframe(kl = kl_divergence(tpm_ko, tpm_mieru)) %>%
            ungroup() %>%
        mutate(
            ko_symbol = input_ko_symbol,
            filter = "complex_genes"
        )


    df_kl_divergence <- bind_rows(df_kl_divergence, df_complex_kl_divergence)
    print(input_ko_symbol)
    print(nrow(df_kl_divergence))

}


df_kl_divergence %>%
    ggplot(aes(x = filter, y = kl, fill = filter)) +
    # geom_violin() +
    # geom_boxplot(width = 0.1, fill = "white", outliers = FALSE) +
    geom_boxplot(width = 0.1, outliers = FALSE) +
    facet_wrap(~ko_symbol, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
        title = "KL divergence of all genes and complex genes",
        x = "KL divergence",
        y = "Density"
    ) +
    ggsave("data/Fig7/kl_divergence_all_genes_vs_complex_genes.png", width = 10, height = 5)
