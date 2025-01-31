###############################################################################
# Aim: To calculate the odds ratio of the enrichment of spliced genes forming complex and not forming complex
###############################################################################
library(tidyverse)
library(janitor)

# ヒトとマウスのホモロジー遺伝子を取得する
if (!file.exists("data/Fig5/mgi_protein_coding_symbols.txt")) {
    col_names <- c("mgi_id", "mouse_symbol", "marker_id", "hgnc", "human_symbol", "ncbi")
    url <- "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt"
    df_homology <- read_tsv(url, col_names = col_names) %>% clean_names()
    df_homology %>%
        mutate(mouse = toupper(mouse_symbol), human = human_symbol) %>%
        select(mouse, human) %>%
        write_tsv("data/Fig5/mgi_homology_symbols.txt")
}
df_homology <- read_tsv("data/Fig5/mgi_homology_symbols.txt")
df_mgi_genes <- read_tsv("data/Fig5/mgi_protein_coding_symbols.txt") %>% mutate(symbol = toupper(symbol))
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

# Complexのデータの、ヒトの遺伝子をマウスの遺伝子に変更する
df_complextab <- read_csv("data/Fig5/complextab_go_symbol_organism.csv")
df_complextab <- df_complextab %>%
    left_join(df_homology, by = c("symbol" = "human")) %>% # human symbolに対応するmouse symbolを結合
    mutate(symbol = ifelse(organism == "human" & !is.na(mouse), mouse, symbol)) %>% # humanの場合のみsymbolを変換
    select(-mouse)

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(target_symbol = toupper(target_symbol)) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()
events <- df_all$event %>% unique()


###############################################################################
# ヒトとマウスの複合体
###############################################################################

complex_genes <- df_complextab %>%
    pull(symbol) %>%
    unique()
# input_ko_symbol <- ko_symbols[1]
# input_event <- events[1]
results_fisher <- tibble()
for (input_ko_symbol in ko_symbols) {
    for (input_event in events) {
        spliced_genes <- df_spliced_genes %>%
            filter(ko_symbol == input_ko_symbol, event == input_event) %>%
            pull(target_symbol)
        non_spliced_genes <- df_mgi_genes %>%
            filter(!symbol %in% spliced_genes) %>%
            pull(symbol)

        print(c(input_ko_symbol, input_event, length(spliced_genes), length(non_spliced_genes)))

        overlap_spliced_complex <- spliced_genes %in% complex_genes
        overlap_non_spliced_complex <- non_spliced_genes %in% complex_genes

        a <- sum(overlap_spliced_complex)
        b <- sum(!overlap_spliced_complex)
        c <- sum(overlap_non_spliced_complex)
        d <- sum(!overlap_non_spliced_complex)
        vx <- matrix(c(a, b, c, d), nrow = 2, byrow = T)
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
            non_spliced_genes_not_forming_complex = d
        ))
    }
}

results_fisher %>%
    filter(significance == "YES") %>%
    as.data.frame()

results_fisher %>% write_csv("reports/Fig5/fisher_complextab_human_mouse.csv")

results_fisher %>% count(significance)

###############################################################################
# Plot Odds Ratio
###############################################################################

# p_valueに応じたアスタリスクの列を追加
results_fisher <- results_fisher %>%
    mutate(asterisk = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01 ~ "**",
        p_value <= 0.05 ~ "*",
        TRUE ~ "" # 条件を満たさない場合は空白
    ))

colors <- c("#44ED8B", "#FF2FC1", "#3FAFFF", "#FFE270", "#FF604E")
names(colors) <- c("A3SS", "A5SS", "MXE", "RI", "SE")

# ggplot2で棒グラフとアスタリスクを描画
g_barplot <- results_fisher %>%
    ggplot(aes(x = event, y = odds_ratio, fill = event)) +
    geom_col(position = position_dodge(width = 0.9), color = "#333") + # 棒グラフ
    geom_hline(yintercept = 1, linetype = "dashed", color = "#333") + # y=1に線を描画
    geom_text(
        aes(label = asterisk, y = odds_ratio + 0.1), # アスタリスクをodds_ratioの少し上に配置
        position = position_dodge(width = 0.9),
        vjust = 0
    ) +
    scale_fill_manual(name = "Event", values = colors) +
    theme_bw() +
    # X軸のラベルを45度回転
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Enrichment (odds ratio)", legend = "Event") +
    facet_wrap(~ko_symbol, scales = "fixed", nrow = 2)


ggsave("reports/Fig5/barplot_odds_complextab_human_mouse.jpg", g_barplot, width = 15, height = 8)
ggsave("reports/Fig5/barplot_odds_complextab_human_mouse.pdf", g_barplot, width = 15, height = 8)
