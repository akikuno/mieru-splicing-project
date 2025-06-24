library(tidyverse)
library(circlize)
library(patchwork)

df_enrichr <- read_csv("reports/Fig5/enrichr.csv")

data <- df_enrichr %>%
    select(Term, from = ko_symbol) %>%
    mutate(to = str_remove(Term, " \\(GO:.*"))

top_n <- 10
df_chord <-
    data %>%
    count(to) %>%
    arrange(desc(n)) %>%
    slice_head(n = top_n) %>%
    inner_join(data, by = "to") %>%
    mutate(value = 1) %>%
    distinct() %>%
    # to の文字列に、何個の遺伝子の重複があるのかをアノテーションする
    mutate(to = paste(to, "\n(", n, ")", sep = "")) %>%
    select(from, to, value)

unique(df_chord$to)
nrows <- length(unique(df_chord$from))
ncols <- length(unique(df_chord$to))

jpeg(file = "reports/Fig5/023-circos.jpg", width = 1500, height = 1500, units = "px", res = 300)
par(cex = 0.4, family = "Arial")
circos.par(gap.after = c(rep(3, nrows - 1), 10, rep(3, ncols - 1), 10))
chordDiagram(df_chord,
    transparency = 0.5,
    annotationTrack = c("name", "grid"),
    directional = -1,
    direction.type = c("diffHeight")
)
circos.clear()
dev.off()

pdf(file = "reports/Fig5/023-circos.pdf")
circos.clear()
par(cex = 0.5)
circos.par(gap.after = c(rep(3, nrows - 1), 10, rep(3, ncols - 1), 10))
chordDiagram(df_chord,
    transparency = 0.5,
    annotationTrack = c("name", "grid"),
    directional = -1,
    direction.type = c("diffHeight")
)
dev.off()
