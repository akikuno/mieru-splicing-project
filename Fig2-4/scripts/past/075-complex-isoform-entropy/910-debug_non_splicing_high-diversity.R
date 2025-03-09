###########################################################
# スプライシングを受けていないのに、MIERUと比べてエントロピーが高い遺伝子は、なぜなのか？本当に多様性が増えているのか？
###########################################################

library(tidyverse)
library(patchwork)
library(ggsignif)

tmp_isoforms <- read_tsv("data/Fig7/tpm_isoforms_all.tsv.gz")
tmp_isoforms <- tmp_isoforms %>%
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
    # サンプルにおいて、IsoformのTPMの平均が100以上の遺伝子のみを抽出
    group_by(sample, gene_symbol) %>%
    filter(mean(tpm) >= 100) %>%
    ungroup()

df_isoforms <-
    # すべてのサンプルが条件をみたすもののみを抽出
    tmp_isoforms_tpm %>%
    group_by(group, gene_symbol) %>%
    mutate(n = n_distinct(sample)) %>%
    ungroup() %>%
    filter(n == sample_number) %>%
    select(-c(n, sample_number))

df_isoforms %>% filter(gene_symbol == "Rrm2") %>% as.data.frame()
df_isoforms %>% filter(gene_symbol == "Rrm2") %>% filter(str_detect(sample, "MIERU"))


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

df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

g_list <- list()
df_t_test <- tibble()
input_ko_symbol <- ko_symbols[1]

df_mieru_isoforms <- df_isoforms %>% filter(str_detect(sample, "MIERU"))
df_mieru_entropy <- df_mieru_isoforms %>%
    select(sample, gene_symbol, tpm, group) %>%
    # グループごとに、エントロピーを計算
    group_by(group, gene_symbol) %>%
    mutate(entropy = calculate_entropy(tpm)) %>%
    ungroup() %>%
    select(group, gene_symbol, entropy) %>%
    distinct()

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
        # グループごとに、エントロピーを計算
        group_by(group, gene_symbol) %>%
        mutate(entropy = calculate_entropy(tpm)) %>%
        ungroup() %>%
        select(group, gene_symbol, entropy) %>%
        distinct()

    df_entropy <- inner_join(df_ko_entropy, df_mieru_entropy, by = "gene_symbol", suffix = c("_ko", "_mieru"))
    print(input_ko_symbol)
    print(nrow(df_non_spliced_genes))
    print(t.test(df_entropy$entropy_ko, df_entropy$entropy_mieru, paired = TRUE)$p.value)

    df_entropy_all <- bind_rows(df_entropy_all, df_entropy)
}

# DEBUG
df_entropy_all %>%
    mutate(entropy_diff = entropy_mieru - entropy_ko) %>%
    slice_max(entropy_diff, n = 10)

# Tra2b, Galt
df_all %>% filter(target_symbol == "Eif5a")
df_all %>% filter(ko_symbol == "Qk", target_symbol == "Eif5a")

df_ko_isoforms %>% filter(gene_symbol == "Eif5a") %>% as.data.frame()
df_mieru_isoforms %>% filter(gene_symbol == "Eif5a") %>% as.data.frame()


df_isoforms %>% filter(gene_symbol == "Rrm2") %>% filter(str_detect(sample, "Tra2b")) %>% as.data.frame()
df_isoforms %>% filter(gene_symbol == "Rrm2") %>% filter(str_detect(sample, "MIERU"))

library(dplyr)

# Shannon entropy 計算関数
shannon_entropy <- function(tpm) {
  p <- tpm / sum(tpm)  # 各 isoform の発現割合
  p <- p[p > 0]        # 0を除外（log(0)の回避）
  -sum(p * log2(p))    # Shannon entropy
}

# データ統合
df_combined <- bind_rows(df_ko_isoforms, df_mieru_isoforms)

# 各groupごとにエントロピー計算
df_entropy <- df_combined %>%
  group_by(group, sample) %>%
  summarise(entropy = shannon_entropy(tpm), .groups = "drop")

# groupごとの平均エントロピーを比較
df_entropy_summary <- df_entropy %>%
  group_by(group) %>%
  summarise(mean_entropy = mean(entropy), sd_entropy = sd(entropy), .groups = "drop")

# 結果の表示
print(df_entropy_summary)

t.test(entropy ~ group, data = df_entropy)


library(corrr)

# TPMの順位を計算
df_rank <- df_combined %>%
  group_by(group, sample) %>%
  mutate(rank = rank(-tpm)) %>%  # 発現量が高いものから順にランク付け
  ungroup() %>%
  select(group, rank) %>%
  pivot_wider(names_from = group, values_from = rank, values_fill = NA) %>%
  filter(!is.na(MIERU), !is.na(Ybx1))


# スピアマン順位相関
correlation_value <- cor(df_rank$Ybx1, df_rank$MIERU, method = "spearman", use = "complete.obs")

print(correlation_value)
