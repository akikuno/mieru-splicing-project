library(tidyverse)

df_se <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv") %>% filter(fdr < 0.05, abs(dpsi) > 0.1) %>% filter(event == "SE")

df_deg <- read_csv("data/degs/ko_vs_mieru.csv")

df_se_deg <- df_se %>% left_join(df_deg, by = c("ko_symbol", "target_symbol"))


