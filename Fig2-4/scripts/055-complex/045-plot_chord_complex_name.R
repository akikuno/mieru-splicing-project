########################################
# 具体的にどの様なタンパク質複合体が影響を受けやすいのか？
########################################
library(tidyverse)
library(circlize)
library(patchwork)

library(tidyverse)
library(janitor)

# ヒトとマウスのホモロジー遺伝子を取得する
if (!file.exists("data/Fig5/mgi_protein_coding_symbols.txt")) {
    col_names <- c("mgi_id", "mouse_symbol", "marker_id", "hgnc", "human_symbol", "ncbi")
    url <- "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt"
    df_homology <- read_tsv(url, col_names = col_names) %>% clean_names()
    df_homology %>%
        mutate(mouse = toupper(mouse_symbol), human = human_symbol) %>%
        select(mouse, human) %>%
        write_tsv("data/Fig5/mgi_homology_symbols.txt")
}
df_homology <- read_tsv("data/Fig5/mgi_homology_symbols.txt")
df_mgi_genes <- read_tsv("data/Fig5/mgi_protein_coding_symbols.txt") %>% mutate(symbol = toupper(symbol))
df_all <- read_csv("data/rmats/all_events_ko_target_fdr_dpsi.csv")

# Complexのデータの、ヒトの遺伝子をマウスの遺伝子に変更する
df_complextab <- read_csv("data/Fig5/complextab_name_go_symbol_organism.csv")
df_complextab <- df_complextab %>%
    left_join(df_homology, by = c("symbol" = "human")) %>% # human symbolに対応するmouse symbolを結合
    mutate(symbol = ifelse(organism == "human" & !is.na(mouse), mouse, symbol)) %>% # humanの場合のみsymbolを変換
    select(-mouse, -go, -organism) %>%
    distinct() %>%
    add_count(name)

df_spliced_genes <- df_all %>%
    filter(fdr < 0.05, abs(dpsi) > 0.1) %>%
    select(event, ko_symbol, target_symbol) %>%
    mutate(target_symbol = toupper(target_symbol)) %>%
    distinct()

ko_symbols <- df_all$ko_symbol %>% unique()


for (input_ko_symbol in ko_symbols) {
    df_spliced_genes %>%
        filter(ko_symbol == input_ko_symbol) %>%
        select(target_symbol) %>%
        distinct() %>%
        inner_join(df_complextab, by = c("target_symbol" = "symbol"), relationship = "many-to-many") %>%
        group_by(target_symbol, name) %>%
        add_count() %>%
        mutate(percentage = nn / n * 100) %>%
        arrange(desc(percentage)) %>%
        select(name) %>%
        distinct() %>%
        write_tsv(paste0("data/Fig5/", input_ko_symbol, ".txt"))
}

df_enrichr <- read_csv("reports/Fig4/enrichr.csv")

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

jpeg(file = "reports/Fig4/023-circos.jpg", width = 1500, height = 1500, units = "px", res = 300)
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

pdf(file = "reports/Fig4/023-circos.pdf")
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
