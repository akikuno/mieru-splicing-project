library(tidyverse)
library(pheatmap)

df_rmats <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")
df_spliceosome_family <- read_csv("data/Fig7/spliceosome_family_mouse.csv")

df_spliced_genes <- df_rmats %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(ko_symbol, target_symbol, dpsi) %>%
    distinct()


ko_symbols <- df_rmats$ko_symbol %>% unique()
input_ko_symbol <- ko_symbols[1]

df_spliceosome_family_disturbed <- tibble()

for (input_ko_symbol in ko_symbols) {
    df_spliceosome_family_disturbed <-
        df_spliced_genes %>%
        filter(ko_symbol == input_ko_symbol) %>%
        select(target_symbol, dpsi) %>%
        distinct() %>%
        inner_join(df_spliceosome_family, by = c("target_symbol" = "mouse"), relationship = "many-to-many") %>%
        mutate(ko_symbol = input_ko_symbol) %>%
        bind_rows(df_spliceosome_family_disturbed)
}

df_spliceosome_family_disturbed %>% select(target_symbol) %>% distinct() %>% arrange(target_symbol)
df_spliceosome_family %>% select(mouse) %>% distinct() %>% arrange(mouse)

order_target_symbol <- df_spliceosome_family %>% select(mouse) %>% distinct() %>% pull(mouse)

mat_dpsi <-
    df_spliceosome_family_disturbed %>%
    right_join(select(df_spliceosome_family, mouse), by = c("target_symbol" = "mouse"), relationship = "many-to-many") %>%
    mutate(dpsi = ifelse(is.na(dpsi), 0, dpsi)) %>%
    arrange(desc(dpsi)) %>%
    select(ko_symbol, target_symbol, dpsi) %>%
    # abs(dpsi)が最大のものを選ぶ
    group_by(ko_symbol, target_symbol) %>%
    slice_max(abs(dpsi), n=1) %>%
    ungroup() %>%
    distinct() %>%
    # matrix変換
    pivot_wider(
        names_from = ko_symbol,
        values_from = dpsi
    ) %>%
    # NAを0に置換
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    arrange(match(target_symbol, order_target_symbol)) %>%
    column_to_rownames("target_symbol") %>%
    as.matrix()


dir.create("reports/Fig7", showWarnings = FALSE)
pdf("reports/Fig7/heatmap_rbp.pdf", width = 10, height = 20)
pheatmap(mat_dpsi, scale = "none", cluster_rows = FALSE, color = colorRampPalette(c("#223a70", "#FFF", "#ec6d51"))(100), fontsize = 16)
dev.off()

svg("reports/Fig7/heatmap_rbp.svg", width = 10, height = 20)
pheatmap(mat_dpsi, scale = "none", cluster_rows = FALSE, color = colorRampPalette(c("#223a70", "#FFF", "#ec6d51"))(100), fontsize = 16)
dev.off()
