###############################################################################
# Aim: To calculate the odds ratio of the enrichment of spliced genes forming complex and not forming complex
###############################################################################
library(tidyverse)

df_mgi_genes <- read_tsv("data/Fig5/mgi_protein_coding_symbols.txt")
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

df_complextab <- read_csv("data/Fig5/complextab_human_mouse.csv")

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()
events <- df_all$event %>% unique()


###############################################################################
# ヒトとマウスの複合体
###############################################################################
results_fisher <- read_csv("reports/Fig5/fisher_complextab_human_mouse.csv")
results_fisher_by_events <- read_csv("reports/Fig5/fisher_complextab_human_mouse_by_events.csv")


###############################################################################
# Plot Odds Ratio
###############################################################################

# -----------------------------------------------------------------------------
# KOごと
# -----------------------------------------------------------------------------
# p_valueに応じたアスタリスクの列を追加
results_fisher <- results_fisher %>%
    mutate(asterisk = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01 ~ "**",
        p_value <= 0.05 ~ "*",
        TRUE ~ "" # 条件を満たさない場合は空白
    ))

# ggplot2で棒グラフとアスタリスクを描画
g_barplot <- results_fisher %>%
    ggplot(aes(x = ko_symbol, y = odds_ratio)) +
    geom_col(position = position_dodge(width = 0.9), color = "#333", fill = "#AAA") + # 棒グラフ
    geom_hline(yintercept = 1, linetype = "dashed", color = "#333") + # y=1に線を描画
    geom_text(
        aes(label = asterisk, y = odds_ratio + 0.05), # アスタリスクをodds_ratioの少し上に配置
        position = position_dodge(width = 0.9),
        vjust = 0,
        size = 6
    ) +
    theme_bw() +
    # X軸のラベルを45度回転
    theme(
        text = element_text(size = 28), # 全体のフォントサイズを大きく
        axis.text.x = element_text(angle = 45, hjust = 1), # X軸のラベルを45度回転してフォントを大きく
    ) +
    labs(x = "", y = "Enrichment (odds ratio)", fill = "SF-KO")


ggsave("reports/Fig5/barplot_odds_complextab_human_mouse.jpg", g_barplot, width = 15, height = 8)
ggsave("reports/Fig5/barplot_odds_complextab_human_mouse.pdf", g_barplot, width = 15, height = 8)

# -----------------------------------------------------------------------------
# KOごと + Eventごと
# -----------------------------------------------------------------------------

# p_valueに応じたアスタリスクの列を追加
results_fisher_by_events <- results_fisher_by_events %>%
    mutate(asterisk = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01 ~ "**",
        p_value <= 0.05 ~ "*",
        TRUE ~ "" # 条件を満たさない場合は空白
    ))

colors <- c("#44ED8B", "#FF2FC1", "#3FAFFF", "#FFE270", "#FF604E")
names(colors) <- c("A3SS", "A5SS", "MXE", "RI", "SE")

# ggplot2で棒グラフとアスタリスクを描画
g_barplot <- results_fisher_by_events %>%
    ggplot(aes(x = event, y = odds_ratio, fill = event)) +
    geom_col(position = position_dodge(width = 0.9), color = "#333") + # 棒グラフ
    geom_hline(yintercept = 1, linetype = "dashed", color = "#333") + # y=1に線を描画
    geom_text(
        aes(label = asterisk, y = odds_ratio + 0.05), # アスタリスクをodds_ratioの少し上に配置
        position = position_dodge(width = 0.9),
        vjust = 0,
        size = 6
    ) +
    scale_fill_manual(name = "Event", values = colors) +
    theme_bw() +
    # X軸のラベルを45度回転
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Enrichment (odds ratio)", legend = "Event") +
    facet_wrap(~ko_symbol, scales = "fixed", nrow = 2)


ggsave("reports/Fig5/barplot_odds_complextab_human_mouse_by_events.jpg", g_barplot, width = 15, height = 8)
ggsave("reports/Fig5/barplot_odds_complextab_human_mouse_by_events.pdf", g_barplot, width = 15, height = 8)
