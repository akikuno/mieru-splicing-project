library(tidyverse)
library(patchwork)
library(ggsignif)

df_isoforms <- read_tsv("data/Fig7/tpm_isoforms_all.tsv.gz") %>%
    mutate(genes = toupper(gene_symbol)) %>%
    mutate(group = str_remove(sample, "_.*$")) %>%
    # グループにおいて、IsoformのTPMの平均が1以上の遺伝子のみを抽出
    group_by(group, genes) %>%
    filter(mean(tpm) >= 1) %>%
    ungroup()

ko_symbols <- df_isoforms %>%
    filter(group != "MIERU") %>%
    pull(group) %>%
    unique()

###########################################################
# エントロピーによるisoformの多様性を検定 (全遺伝子)
###########################################################

calculate_entropy <- function(values) {
    total_sum <- sum(values)
    proportions <- values / total_sum
    -sum(proportions * log(proportions), na.rm = TRUE)
}

g_list <- list()
df_t_test <- tibble()
input_ko_symbol <- ko_symbols[1]

df_mieru_isoforms <- df_isoforms %>% filter(str_detect(sample, "MIERU"))
df_mieru_entropy <- df_mieru_isoforms %>%
    select(sample, genes, tpm, group) %>%
    # グループごとに、エントロピーを計算
    group_by(group, genes) %>%
    mutate(entropy = calculate_entropy(tpm)) %>%
    ungroup() %>%
    select(group, genes, entropy) %>%
    distinct()

for (input_ko_symbol in ko_symbols) {
    df_ko_isoforms <- df_isoforms %>% filter(str_detect(sample, input_ko_symbol))

    df_entropy <- tibble()

    df_ko_entropy <- df_ko_isoforms %>%
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
        average_entropy_ko = mean(df_ko_entropy$entropy),
        average_entropy_mieru = mean(df_mieru_entropy$entropy)
    ))
}

g_wrap <- wrap_plots(g_list, ncol = 5)

ggsave("reports/Fig7/isoform_diversity_all.png", g_wrap, width = 20, height = 10)
df_t_test %>%
    write_csv("reports/Fig7/isoform_diversity_t_test.csv")


###########################################################
# エントロピーによるisoformの多様性を検定 (非スプライシング遺伝子)
###########################################################

df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")
df_non_spliced_genes <- df_all %>%
    filter(fdr > 0.05, abs(dpsi) < 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(genes = toupper(target_symbol)) %>%
    select(genes) %>%
    distinct()
events <- df_all$event %>% unique()


df_mieru_entropy <- inner_join(df_mieru_entropy, df_non_spliced_genes, by = "genes")

g_list_non_significants <- list()
for (input_ko_symbol in ko_symbols) {
    df_ko_isoforms <- df_isoforms %>%
        filter(str_detect(sample, input_ko_symbol)) %>%
        inner_join(df_non_spliced_genes, by = "genes")

    df_entropy <- tibble()

    df_ko_entropy <- df_ko_isoforms %>%
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

    g_list_non_significants[[input_ko_symbol]] <- g_plot
    df_t_test <- bind_rows(df_t_test, tibble(
        ko_symbol = input_ko_symbol,
        p_value = t.test(df_ko_entropy$entropy, df_mieru_entropy$entropy)$p.value,
        average_entropy_ko = mean(df_ko_entropy$entropy),
        average_entropy_mieru = mean(df_mieru_entropy$entropy)
    ))
}

g_wrap <- wrap_plots(g_list_non_significants, ncol = 5)

ggsave("reports/Fig7/isoform_diversity_non_significants.png", g_wrap, width = 20, height = 10)
df_t_test %>%
    write_csv("reports/Fig7/isoform_diversity_non_significants_t_test.csv")



###########################################################
# エントロピーによるisoformの多様性を検定 (スプライシング遺伝子)
###########################################################

df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")
df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(genes = toupper(target_symbol)) %>%
    select(genes) %>%
    distinct()
events <- df_all$event %>% unique()


df_mieru_entropy <- inner_join(df_mieru_entropy, df_spliced_genes, by = "genes")

g_list_significants <- list()
for (input_ko_symbol in ko_symbols) {
    df_ko_isoforms <- df_isoforms %>%
        filter(str_detect(sample, input_ko_symbol)) %>%
        inner_join(df_spliced_genes, by = "genes")

    df_entropy <- tibble()

    df_ko_entropy <- df_ko_isoforms %>%
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

    g_list_significants[[input_ko_symbol]] <- g_plot
    df_t_test <- bind_rows(df_t_test, tibble(
        ko_symbol = input_ko_symbol,
        p_value = t.test(df_ko_entropy$entropy, df_mieru_entropy$entropy)$p.value,
        average_entropy_ko = mean(df_ko_entropy$entropy),
        average_entropy_mieru = mean(df_mieru_entropy$entropy)
    ))
}

g_wrap <- wrap_plots(g_list_significants, ncol = 5)

ggsave("reports/Fig7/isoform_diversity_significants.png", g_wrap, width = 20, height = 10)
df_t_test %>%
    write_csv("reports/Fig7/isoform_diversity_significants_t_test.csv")




###########################################################
# エントロピーによるisoformの多様性を検定 (２つのFigをまとめる)
###########################################################

df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(genes = toupper(target_symbol)) %>%
    select(genes) %>%
    distinct()


df_non_spliced_genes <- df_all %>%
    filter(fdr > 0.05, abs(dpsi) < 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(genes = toupper(target_symbol)) %>%
    select(genes) %>%
    distinct()


df_ko_entropy_spliced <- df_isoforms %>%
    inner_join(df_spliced_genes, by = "genes") %>%
    # グループごとに、エントロピーを計算
    group_by(group, genes) %>%
    mutate(entropy = calculate_entropy(tpm)) %>%
    ungroup() %>%
    select(group, genes, entropy) %>%
    distinct() %>%
    mutate(category = "SAG")

df_ko_entropy_non_spliced <- df_isoforms %>%
    inner_join(df_non_spliced_genes, by = "genes") %>%
    # グループごとに、エントロピーを計算
    group_by(group, genes) %>%
    mutate(entropy = calculate_entropy(tpm)) %>%
    ungroup() %>%
    select(group, genes, entropy) %>%
    distinct() %>%
    mutate(category = "non-SAG")

df_ko_entropy <- bind_rows(df_ko_entropy_spliced, df_ko_entropy_non_spliced)


g_plot <- ggplot(df_ko_entropy, aes(x = category, y = entropy, fill = group)) +
        geom_violin() +
        geom_boxplot(width = 0.1, fill = "white") +
        geom_signif(
            comparisons = list(c("MIERU", input_ko_symbol)),
            test = "t.test", test.args = list(paired = TRUE), na.rm = FALSE, map_signif_level = TRUE, col = "black", step_increase = 0.1
        ) +
        scale_fill_manual(values = fill_values) +
        labs(x = "", y = "Isoform diversity (entropy)") +
        theme_bw() +
        facet_wrap(~group)
