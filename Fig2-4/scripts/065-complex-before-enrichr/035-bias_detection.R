# 複合体を形成する遺伝子にスプライシングが受けやすい遺伝子が濃縮されていたら、DAGが濃縮されるのは当然なので、そのバイアスの有無を検証する
library(tidyverse)
library(janitor)

df_complextab <- read_csv("data/Fig5/complextab_go_symbol_organism.csv")

col_names <- c("symbol", "refseq")
df_refflat_hs <- read_tsv("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz", col_names = col_names) %>% select(col_names)
df_refflat_mm <- read_tsv("https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz", col_names = col_names) %>% select(col_names)

df_num_isoform_hs <- df_refflat_hs %>% filter(str_detect(refseq, "^NM_")) %>% count(symbol)
df_num_isoform_mm <- df_refflat_mm %>% filter(str_detect(refseq, "^NM_")) %>% count(symbol) %>% mutate(symbol = toupper(symbol))

df_complextab <- read_csv("data/Fig5/complextab_go_symbol_organism.csv")
df_corum <- read_csv("data/Fig5/corum_go_symbol_organism.csv")

db <- "ComplexTab"
og <- "human"
df_results <- tibble()
for (db in c("ComplexTab", "CORUM")) {
    for (og in c("human", "mouse")) {
        if (db == "ComplexTab") {
            df_complex <- df_complextab
        } else {
            df_complex <- df_corum
        }
        if (og == "human") {
            df_num_isoform <- df_num_isoform_hs
        } else {
            df_num_isoform <- df_num_isoform_mm
        }

        df_complex_genes <- df_complex %>% filter(organism == og) %>% select(symbol) %>% distinct()
        num_isoform_with_complex <- df_num_isoform %>% filter(symbol %in% df_complex_genes$symbol) %>% pull(n)
        num_isoform_without_complex <- df_num_isoform %>% filter(!symbol %in% df_complex_genes$symbol) %>% pull(n)
        p_value <- t.test(num_isoform_with_complex, num_isoform_without_complex)$p.value

        mean_isoform_with_complex <- mean(num_isoform_with_complex)
        mean_isoform_without_complex <- mean(num_isoform_without_complex)
        print(c(db, og, p_value, mean_isoform_with_complex, mean_isoform_without_complex))

        df_results <- bind_rows(df_results, tibble(
            database = db,
            organism = og,
            type = c(rep("complex", length(num_isoform_with_complex)), rep("non-complex", length(num_isoform_without_complex))),
            number_of_isoform = c(num_isoform_with_complex, num_isoform_without_complex),
            mean_number_of_isoform = c(rep(mean_isoform_with_complex, length(num_isoform_with_complex)), rep(mean_isoform_without_complex, length(num_isoform_without_complex))),
            p_value = p_value
        ))
    }
}

g_plot <-
    df_results %>%
    ggplot(aes(x = type, y = log2(number_of_isoform), fill = type)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill="white") +
    theme_bw() +
    facet_grid(database ~ organism)


df_report <- df_results %>% select(-number_of_isoform) %>% distinct() %>%
    mutate(significance = ifelse(p_value < 0.05, "YES", "NO"))

ggsave("reports/Fig6/complex_bias_number_of_isoform.jpg", g_plot, width = 10, height = 5)
ggsave("reports/Fig6/complex_bias_number_of_isoform.pdf", g_plot, width = 10, height = 5)
write_csv(df_report, "reports/Fig6/complex_bias_number_of_isoform.csv")
