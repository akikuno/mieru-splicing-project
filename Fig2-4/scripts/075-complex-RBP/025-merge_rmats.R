library(tidyverse)

df_rmats <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")
df_spliceosome_process <- read_csv("data/Fig7/spliceosome_process_mouse.csv")
df_spliceosome_family <- read_csv("data/Fig7/spliceosome_family_mouse.csv")

# df_spliceosome_family_rmats <-
#     inner_join(df_rmats, df_spliceosome_family, by = c("target_symbol" = "mouse"), relationship = "many-to-many")

# df_spliceosome_family %>% select(mouse) %>% distinct()
# df_spliceosome_family_rmats %>% select(target_symbol) %>% distinct()

# df_spliceosome_family_rmats <-
#     df_spliceosome_family_rmats %>%
#     mutate(is_significant = ifelse(fdr < 0.05 & abs(dpsi) > 0.1, TRUE, FALSE)) %>%
#     select(-fdr)


df_spliced_genes <- df_rmats %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    distinct()


ko_symbols <- df_rmats$ko_symbol %>% unique()
events <- df_rmats$event %>% unique()

input_ko_symbol <- ko_symbols[1]
input_event <- events[1]

df_spliceosome_family_disturbed <- tibble()

for (input_ko_symbol in ko_symbols) {
    df_spliceosome_family_disturbed <-
        df_spliced_genes %>%
        filter(ko_symbol == input_ko_symbol) %>%
        select(target_symbol) %>%
        distinct() %>%
        inner_join(df_spliceosome_family, by = c("target_symbol" = "mouse"), relationship = "many-to-many") %>%
        mutate(ko_symbol = input_ko_symbol) %>%
        bind_rows(df_spliceosome_family_disturbed)
}

df_spliceosome_family_disturbed %>% select(target_symbol) %>% distinct() %>% arrange(target_symbol)
df_spliceosome_family %>% select(mouse) %>% distinct() %>% arrange(mouse)

df_spliceosome_family_counts <-
    df_spliceosome_family_disturbed %>%
    select(ko_symbol, target_symbol, spliceosome) %>%
    distinct() %>%
    count(target_symbol) %>%
    right_join(df_spliceosome_family, by = c("target_symbol" = "mouse"), relationship = "many-to-many") %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    arrange(desc(n))


df_spliceosome_family_counts %>%
    select(target_symbol, spliceosome, n) %>%
    distinct() %>%
    ggplot(aes(x = n, y = spliceosome, color= spliceosome)) +
    geom_jitter(size=3, width = 0) +
    theme_bw()


df_spliceosome_family_counts %>%
    select(target_symbol, family, n) %>%
    distinct() %>%
    ggplot(aes(x = n, y = family, color= family)) +
    geom_jitter(size=3, width = 0) +
    theme_bw()
