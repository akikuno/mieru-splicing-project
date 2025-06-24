library(DESeq2)
library(tidyverse)

dir.create("data/degs", showWarnings = FALSE)

df <- read_tsv("data/counts/featurecounts_gene_name.tsv.gz")

count <- as.matrix(select(df, -c(Geneid, Length)))
rownames(count) <- df$Geneid

dim(count)

###########################################################
# Pre-filtering
###########################################################

keep <- rowSums(count) >= 10
count_filtered <- count[keep, ]

dim(count_filtered)


###########################################################
# MIERU vs each KO
###########################################################

ko_symbols <- colnames(count_filtered) %>% str_remove("_KO_[0-9]+$") %>% unique()
ko_symbols <- ko_symbols[!grepl("^MIERU_CH_", ko_symbols)]

results_ko_degs <- map_dfr(ko_symbols, function(input_ko_symbol) {

    print(input_ko_symbol)

    df_ko_symbol <- df %>% select(matches(c("Geneid", input_ko_symbol, "MIERU_CH")))

    count_ko_symbol <- as_tibble(count_filtered) %>% select(matches(c(input_ko_symbol, "MIERU_CH"))) %>% as.matrix()
    rownames(count_ko_symbol) <- rownames(count_filtered)

    group <-
        colnames(count_ko_symbol) %>%
        str_remove("_[0-9]+$") %>%
        as_tibble() %>%
        pull(value)

    group <- data.frame(con = factor(group))

    dds <- DESeqDataSetFromMatrix(
        countData = count_ko_symbol,
        colData = group,
        design = ~con
    )

    degs <- DESeq(dds) %>% DESeq2::results()

    results <-
        degs %>%
        as.data.frame() %>%
        rownames_to_column("Geneid") %>%
        filter(padj < 0.05) %>%
        as_tibble() %>%
        inner_join(df_ko_symbol, by = "Geneid")

    print(nrow(results))

    results %>%
        rowwise() %>%
        mutate(
            ko_ave = mean(c_across(matches(input_ko_symbol)), na.rm = TRUE),
            control_ave = mean(c_across(starts_with("MIERU_CH")), na.rm = TRUE)
        ) %>%
        ungroup() %>%
        mutate(up_down = ifelse(ko_ave > control_ave, "ko_up", "ko_down")) %>%
        mutate(log2FoldChange = ifelse(ko_ave > control_ave, abs(log2FoldChange), -abs(log2FoldChange))) %>%
        mutate(ko_symbol = input_ko_symbol) %>%
        rename(Geneid = "target_symbol") %>%
        select(ko_symbol, target_symbol, baseMean, log2FoldChange, padj, up_down)
})

write_csv(results_ko_degs, "data/degs/ko_vs_mieru.csv")
