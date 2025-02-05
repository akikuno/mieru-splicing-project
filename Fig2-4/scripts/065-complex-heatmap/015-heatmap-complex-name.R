###############################################################################
# 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合をヒートマップで可視化
# X軸：各遺伝子の各イベント
# Y軸：顕著な複合体名
# 値： 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合 （例：複合体形成遺伝子が1つで、その複合体には３つの遺伝子が係る場合→ 1/3）
###############################################################################

library(tidyverse)


df_complextab <- read_tsv("data/Fig5/complextab_mouse.tsv")

df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")
df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.01) %>%
    select(event, ko_symbol, target_symbol) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()
events <- df_all$event %>% unique()

input_ko_symbol <- ko_symbols[1]
input_event <- events[1]

for (input_ko_symbol in ko_symbols) {
df_spliced_genes %>%
    filter(ko_symbol == input_ko_symbol) %>%
    select(target_symbol) %>%
    distinct() %>%
    inner_join(df_complextab, by = c("target_symbol" = "symbol"), relationship = "many-to-many") %>%
    group_by(target_symbol, name) %>%
    add_count(name = "count") %>%
    mutate(percentage = count / n * 100) %>%
    arrange(desc(percentage)) %>%
    filter(n > 1) %>%
    print(.)
}

# count > 1のものがないので、ヒートマップは使えなさそう
# Complex PortalにアノテーションされているGOで、共通性の高いものを探してみる？
