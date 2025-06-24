###############################################################################
# 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合をヒートマップで可視化
# X軸：各遺伝子の各イベント
# Y軸：顕著な複合体名
# 値： 複合体形成遺伝子とその複合体に関わる遺伝子数の重複割合 （例：複合体形成遺伝子が1つで、その複合体には３つの遺伝子が係る場合→ 1/3）
###############################################################################

library(tidyverse)
library(readxl)
library(janitor)

df_spliceosome_process <- read_excel("data/Fig7/spliceosome_gene_list.xlsx", sheet = 1) %>% clean_names()
df_spliceosome_family <- read_excel("data/Fig7/spliceosome_gene_list.xlsx", sheet = 2) %>% clean_names()
df_spliceosome_family <- read_csv("data/Fig7/spliceosome_gene_list.csv") %>% clean_names()

df_spliceosome_process_longer <-
    df_spliceosome_process %>%
    select(-mouse) %>%
    pivot_longer(cols = -c(human, spliceosome), names_to = "process", values_to = "bool") %>%
    filter(!is.na(bool)) %>%
    select(-bool)
    # mutate(bool = ifelse(is.na(bool), FALSE, TRUE))

df_spliceosome_family_longer <-
    df_spliceosome_family %>%
    pivot_longer(cols = -c(human, spliceosome), names_to = "family", values_to = "bool") %>%
    filter(!is.na(bool)) %>%
    select(-bool)
    # mutate(bool = ifelse(is.na(bool), FALSE, TRUE))

# Human-Mouse orthologs
df_homology <- read_tsv("data/Fig5/mgi_homology_symbols.txt")

df_spliceosome_process_longer_mouse <-
    df_spliceosome_process_longer %>%
    inner_join(df_homology) %>%
    select(-human)


df_spliceosome_family_longer_mouse <-
    df_spliceosome_family_longer %>%
    inner_join(df_homology) %>%
    select(-human)

df_spliceosome_process_longer_mouse %>% count(process)
df_spliceosome_family_longer_mouse %>% count(family)


write_csv(df_spliceosome_process_longer_mouse, "data/Fig7/spliceosome_process_mouse.csv")
write_csv(df_spliceosome_family_longer_mouse, "data/Fig7/spliceosome_family_mouse.csv")


df_spliceosome_family_longer_mouse <-
    df_spliceosome_family_longer %>%
    left_join(df_homology) %>%
    select(human, mouse) %>%
    distinct() %>%
    arrange(mouse) %>%
    write_csv("data/Fig7/spliceosome_human_mouse_converter.csv")
