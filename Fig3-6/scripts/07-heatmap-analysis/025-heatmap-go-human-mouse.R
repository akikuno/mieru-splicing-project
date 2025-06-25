###############################################################################
# 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合をヒートマップで可視化
# X軸：各遺伝子の各イベント
# Y軸：顕著な複合体名
# 値： 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合 （例：複合体形成遺伝子が1つで、その複合体には３つの遺伝子が係る場合→ 1/3）
###############################################################################

library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(pheatmap)

df_mgi_genes <- read_tsv("data/Fig5/mgi_protein_coding_symbols.txt")
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

df_complextab <- read_csv("data/Fig5/complextab_human_mouse.csv")
df_go_annotation <- read_csv("data/Fig6/go_annotation_human_mouse.csv")

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()
events <- df_all$event %>% unique()

input_ko_symbol <- ko_symbols[1]
input_event <- events[1]

if(!file.exists("data/Fig6/go_enrichments.csv")) {
    go_enrichments <- tibble()

    genes_complex <- df_complextab$symbol %>% unique()
    for (input_ko_symbol in ko_symbols) {
        spliced_genes <- df_spliced_genes %>%
            filter(ko_symbol == input_ko_symbol) %>%
            select(symbol = target_symbol) %>%
            distinct()
        spliced_genes_with_complex <- spliced_genes %>%
            filter(symbol %in% genes_complex) %>%
            pull(symbol) %>%
            unique()

        go_enrichment <- enrichGO(
            gene = spliced_genes_with_complex,
            universe = genes_complex,
            OrgDb = org.Mm.eg.db,
            keyType = "SYMBOL",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.01,
            qvalueCutoff = 0.05
        )
        go_enrichments <- as_tibble(go_enrichment) %>%
            mutate(ko_symbol = input_ko_symbol) %>%
            bind_rows(go_enrichments)
    }

    go_enrichments <- go_enrichments %>% janitor::clean_names()

    write_csv(go_enrichments, "data/Fig6/go_enrichments.csv")
}

go_enrichments <- read_csv("data/Fig6/go_enrichments.csv")

go_overlap <-
    go_enrichments %>%
    add_count(description, name = "go_count") %>%
    filter(go_count > 9) %>%
    select(ko_symbol, description, fold_enrichment) %>%
    complete(ko_symbol, description, fill = list(fold_enrichment = 0)) %>%
    distinct()

mat_overlap <-
    go_overlap %>%
    pivot_wider(names_from = description, values_from = fold_enrichment) %>%
    column_to_rownames("ko_symbol") %>%
    as.matrix() %>%
    t()

dir.create("reports/Fig6", showWarnings = FALSE)
pdf("reports/Fig6/heatmap_go_human_mouse.pdf", width = 15, height = 10)
pheatmap(mat_overlap, scale = "none", color = colorRampPalette(c("white", "#fd7e00"))(10), fontsize = 20)
dev.off()

jpeg("reports/Fig6/heatmap_go_human_mouse.jpg", width = 15, height = 10, units = "in", res = 600)
pheatmap(mat_overlap, scale = "none", color = colorRampPalette(c("white", "#fd7e00"))(10), fontsize = 20)
dev.off()
