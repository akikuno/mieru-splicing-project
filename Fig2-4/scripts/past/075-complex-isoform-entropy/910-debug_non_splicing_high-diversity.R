###########################################################
# スプライシングを受けていないのに、MIERUと比べてエントロピーが高い遺伝子は、なぜなのか？本当に多様性が増えているのか？
###########################################################

library(tidyverse)
library(patchwork)
library(ggsignif)

tmp_isoforms <- read_tsv("data/Fig7/tpm_isoforms_all.tsv.gz") %>%
    mutate(group = str_remove(sample, "_.*$"))

tmp_isoforms_counts <-
    tmp_isoforms %>%
    select(sample, group) %>%
    distinct() %>%
    group_by(group) %>%
    add_count(group, name = "sample_number") %>%
    ungroup() %>%
    select(group, sample_number) %>%
    inner_join(tmp_isoforms, by = "group", relationship = "many-to-many") %>%
    distinct()

tmp_isoforms_tpm <-
    tmp_isoforms_counts %>%
    # サンプルにおいて、IsoformのTPMの平均が10以上の遺伝子のみを抽出
    group_by(sample, gene_symbol) %>%
    filter(mean(tpm) >= 10) %>%
    ungroup()

df_isoforms <-
    # すべてのサンプルが条件をみたすもののみを抽出
    tmp_isoforms_tpm %>%
    group_by(group, gene_symbol) %>%
    mutate(n = n_distinct(sample)) %>%
    ungroup() %>%
    filter(n == sample_number) %>%
    select(-c(n, sample_number))

# df_isoforms %>% filter(gene_symbol == "Rrm2") %>% as.data.frame()
# df_isoforms %>% filter(gene_symbol == "Rrm2") %>% filter(str_detect(sample, "MIERU"))

# df_isoforms %>% filter(group == "Ybx1" | group == "MIERU") %>% filter(gene_symbol == "Rrm2") %>%
#     ggplot(aes(x = transcript_symbol, y = tpm, color = group)) +
#     geom_point() +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))

ko_symbols <- df_isoforms %>%
    filter(group != "MIERU") %>%
    pull(group) %>%
    unique()

###########################################################
# エントロピーによるisoformの多様性を検定 (非スプライシング遺伝子)
###########################################################

calculate_entropy <- function(values) {
    total_sum <- sum(values)
    proportions <- values / total_sum
    -sum(proportions * log(proportions), na.rm = TRUE)
}

# df_isoforms %>% filter(group == "Ybx1" | group == "MIERU") %>% filter(gene_symbol == "Rrm2") %>%
#     # サンプルごとに、エントロピーを計算
#     group_by(sample, gene_symbol) %>%
#     mutate(entropy = calculate_entropy(tpm)) %>%
#     ungroup() %>%
#     # グループごとに、エントロピーの平均を計算
#     group_by(group) %>%
#     mutate(entropy_mean = mean(entropy)) %>%
#     ungroup() %>%
#     select(group, gene_symbol, entropy_mean) %>%
#     distinct()

df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

g_list <- list()
df_t_test <- tibble()
input_ko_symbol <- ko_symbols[1]

df_mieru_isoforms <- df_isoforms %>% filter(str_detect(sample, "MIERU"))

df_mieru_entropy <- df_mieru_isoforms %>%
    # サンプルごとに、エントロピーを計算
    group_by(sample, gene_symbol) %>%
    mutate(entropy = calculate_entropy(tpm)) %>%
    ungroup() %>%
    # グループごとに、エントロピーの平均を計算
    group_by(group, gene_symbol) %>%
    mutate(entropy_mean = mean(entropy)) %>%
    ungroup() %>%
    select(group, gene_symbol, entropy_mean) %>%
    distinct()

# df_mieru_entropy %>% filter(gene_symbol == "Rrm2") %>% as.data.frame()

# df_mieru_entropy <- df_mieru_isoforms %>%
#     select(sample, gene_symbol, tpm, group) %>%
#     # グループごとに、エントロピーを計算
#     group_by(group, gene_symbol) %>%
#     mutate(entropy = calculate_entropy(tpm)) %>%
#     ungroup() %>%
#     select(group, gene_symbol, entropy) %>%
#     distinct()

df_mieru_entropy <- inner_join(df_mieru_entropy, df_non_spliced_genes, by = "gene_symbol")

df_entropy_all <- tibble()
for (input_ko_symbol in ko_symbols) {

    df_non_spliced_genes <-
        df_all %>%
        filter(ko_symbol == input_ko_symbol) %>%
        filter(fdr > 0.05, abs(dpsi) < 0.1) %>%
        select(event, ko_symbol, target_symbol) %>%
        rename(gene_symbol = target_symbol) %>%
        select(gene_symbol) %>%
        distinct()

    df_ko_isoforms <-
        df_isoforms %>%
        filter(str_detect(sample, input_ko_symbol)) %>%
        inner_join(df_non_spliced_genes, by = "gene_symbol")

    df_ko_entropy <-
        df_ko_isoforms %>%
        # サンプルごとに、エントロピーを計算
        group_by(sample, gene_symbol) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        # グループごとに、エントロピーの平均を計算
        group_by(group, gene_symbol) %>%
        mutate(entropy_mean = mean(entropy)) %>%
        ungroup() %>%
        select(group, gene_symbol, entropy_mean) %>%
        distinct()

    df_entropy <- inner_join(df_ko_entropy, df_mieru_entropy, by = "gene_symbol", suffix = c("_ko", "_mieru"))
    print(input_ko_symbol)
    print(nrow(df_non_spliced_genes))
    print(t.test(df_entropy$entropy_mean_ko, df_entropy$entropy_mean_mieru, paired = TRUE)$p.value)

    df_entropy_all <- bind_rows(df_entropy_all, df_entropy)
}
