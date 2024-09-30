library(tidyverse)
library(patchwork)

######################################
# Exon Length and %GC
######################################

df <- read_csv("data/exon_conservations/exon_length_gc.csv")

summary(df$exon_length)
summary(df$gc)

p_length <- ggplot(df, aes(x=ko_symbol, y=log2(exon_length))) +
    geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    labs(x="RBP KOs", y="log2(Exon length)") +
    theme_bw()

p_gc <- ggplot(df, aes(x=ko_symbol, y=gc)) +
    geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    labs(x="RBP KOs", y="% GC") + 
    theme_bw()

######################################
# Exon Conservation
######################################

df_cons <- read_csv("data/exon_conservations/exon_conservation.csv.gz")

df_cons_mean <-
    df_cons %>%
    group_by(ko_symbol, target_symbol, fdr, dpsi) %>%
    summarise(mean_phylop = mean(phylop)) %>%
    ungroup()

p_cons <- ggplot(df_cons_mean, aes(x=ko_symbol, y=mean_phylop)) +
    geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    labs(x="RBP KOs", y="phyloP Conservation Score") + 
    theme_bw()


######################################
# Save plot
######################################

p <- p_length / p_gc / p_cons + plot_annotation(title = "All SE events", tag_levels = 'A')

dir.create("reports/rmats/exon_length_gc", showWarnings = FALSE)
ggsave("reports/rmats/exon_length_gc/violin_exon_length_gc_conservation.jpg", p, width=10, height=10)
ggsave("reports/rmats/exon_length_gc/violin_exon_length_gc_conservation.pdf", p, width=10, height=10)
