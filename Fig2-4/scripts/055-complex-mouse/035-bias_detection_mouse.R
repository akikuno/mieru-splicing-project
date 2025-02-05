# 複合体を形成する遺伝子にスプライシングが受けやすい遺伝子が濃縮されていたら、DAGが濃縮されるのは当然なので、そのバイアスの有無を検証する
library(tidyverse)
library(janitor)
library(ggsignif)


col_names <- c("symbol", "refseq")
df_refflat_mm <- read_tsv("https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz", col_names = col_names) %>% select(all_of(col_names))

df_num_isoform_mm <- df_refflat_mm %>%
    filter(str_detect(refseq, "^NM_")) %>%
    count(symbol)
df_complextab <- read_tsv("data/Fig5/complextab_mouse.tsv")

df_results <- tibble()

df_complex_genes <- df_complextab %>%
    select(symbol) %>%
    distinct()
num_isoform_with_complex <- df_num_isoform_mm %>%
    filter(symbol %in% df_complex_genes$symbol) %>%
    pull(n)
num_isoform_without_complex <- df_num_isoform_mm %>%
    filter(!symbol %in% df_complex_genes$symbol) %>%
    pull(n)
p_value <- t.test(num_isoform_with_complex, num_isoform_without_complex)$p.value

mean_isoform_with_complex <- mean(num_isoform_with_complex)
mean_isoform_without_complex <- mean(num_isoform_without_complex)
print(c(p_value, mean_isoform_with_complex, mean_isoform_without_complex))

df_results <- tibble(
    type = c(rep("complex", length(num_isoform_with_complex)), rep("non-complex", length(num_isoform_without_complex))),
    number_of_isoform = c(num_isoform_with_complex, num_isoform_without_complex),
    mean_number_of_isoform = c(rep(mean_isoform_with_complex, length(num_isoform_with_complex)), rep(mean_isoform_without_complex, length(num_isoform_without_complex))),
    p_value = p_value
)

# プロット作成
g_plot <-
    df_results %>%
    ggplot(aes(x = type, y = log2(number_of_isoform), fill = type)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    annotate(
        "text",
        x = 1.5, # 2群の中央
        y = max(log2(df_results$number_of_isoform)) + 0.1, # 最大値より少し上
        label = "N.S.",
        size = 5,
        vjust = 0
    ) +
    labs(x = "", y = "log2(Number of isoform)") +
    theme_bw()


df_report <- df_results %>%
    select(-number_of_isoform) %>%
    distinct() %>%
    mutate(significance = ifelse(p_value < 0.05, "YES", "NO"))

ggsave("reports/Fig5/complex_bias_number_of_isoform.jpg", g_plot, width = 10, height = 5)
ggsave("reports/Fig5/complex_bias_number_of_isoform.pdf", g_plot, width = 10, height = 5)
write_csv(df_report, "reports/Fig5/bias_number_of_isoform.csv")
