library(tidyverse)
library(patchwork)

######################################
# Exon Length and %GC
######################################

df <- read_tsv("data/exon_characteristics/exon_length_gc.tsv")

summary(df$exon_length)
summary(df$gc)

p_length <- ggplot(df, aes(x=ko_symbol, y=log2(exon_length))) +
    geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    labs(x="SF KOs", y="log2(Exon length)") +
    theme_bw()

p_gc <- ggplot(df, aes(x=ko_symbol, y=gc)) +
    geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    labs(x="SF KOs", y="% GC") +
    theme_bw()

######################################
# Exon Conservation
######################################

df_cons <- read_tsv("data/exon_characteristics/exon_conservation.tsv")

df_cons_mean <-
    df_cons %>%
    group_by(exon_id) %>%
    summarise(mean_phylop = mean(phylop)) %>%
    ungroup()

p_cons <- ggplot(df_cons_mean, aes(x=ko_symbol, y=mean_phylop)) +
    geom_violin() +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    labs(x="SF KOs", y="phyloP Conservation Score") +
    theme_bw()


######################################
# Save plot
######################################

p <- p_length / p_gc / p_cons + plot_annotation(title = "All SE events", tag_levels = 'A')

ggsave("reports/figure/027-violin_exon_length_gc_conservation.jpg", p, width=10, height=10)
ggsave("reports/figure/027-violin_exon_length_gc_conservation.pdf", p, width=10, height=10)
