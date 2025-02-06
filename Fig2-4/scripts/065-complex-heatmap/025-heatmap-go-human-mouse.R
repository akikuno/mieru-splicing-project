###############################################################################
# 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合をヒートマップで可視化
# X軸：各遺伝子の各イベント
# Y軸：顕著な複合体名
# 値： 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合 （例：複合体形成遺伝子が1つで、その複合体には３つの遺伝子が係る場合→ 1/3）
###############################################################################

library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)

df_mgi_genes <- read_tsv("data/Fig5/mgi_protein_coding_symbols.txt")
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

df_complextab <- read_csv("data/Fig5/complextab_human_mouse.csv")
df_go_annotation <- read_csv("data/Fig6/go_annotation_human_mouse.csv")

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.01) %>%
    select(event, ko_symbol, target_symbol) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()
events <- df_all$event %>% unique()

input_ko_symbol <- ko_symbols[1]
input_event <- events[1]

go_complex_genes <- df_complextab %>% count(go)

genes_complex <- df_complextab$symbol %>% unique()

go_enrichments <- tibble()

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


# go_enrichments %>% filter(ko_symbol == "Strap") %>%
#     select(ko_symbol, ID, Description, FoldEnrichment, qvalue, geneID, GeneRatio)

go_enrichments <- read_csv("data/Fig6/go_enrichments.csv")

go_overlap <-
    go_enrichments %>%
    add_count(description, name = "go_count") %>%
    filter(go_count > 8) %>%
    select(ko_symbol, description, fold_enrichment) %>%
    distinct()

go_overlap <- go_overlap %>% filter(description == "regulation of cell cycle")
# gene_idの要素が90%以上マッチしているものは、同一のGOとして、qvalueがもっとも低いものを選択する


# gene_id をリストに変換
go_overlap <- go_overlap %>%
    mutate(gene_list = strsplit(gene_id, "/"))

# 類似度判定関数（90%以上の遺伝子が一致するか）
similarity_check <- function(x, y) {
    common_genes <- length(intersect(x, y))
    min_length <- min(length(x), length(y))
    return(common_genes / min_length)
}

# 全ペアの組み合わせを取得
pairs <- expand.grid(idx1 = 1:nrow(go_overlap), idx2 = 1:nrow(go_overlap)) %>%
    filter(idx1 < idx2)

# 類似するペアをフィルタリング
edges <-
    pairs %>%
    mutate(similar = pmap_dbl(list(idx1, idx2), function(i, j) {
        similarity_check(go_overlap$gene_list[[i]], go_overlap$gene_list[[j]])
    })) %>%
    filter(similar > 0.75) %>%
    select(idx1, idx2)
