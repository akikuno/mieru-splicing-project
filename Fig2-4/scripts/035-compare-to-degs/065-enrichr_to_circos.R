library(tidyverse)
library(circlize)
library(patchwork)

df_go <- read_csv("reports/Fig3/055-go_rmats_deg.csv")

events <- df_go$event %>% unique()
input_event <- "RI"

plot_list <- list()

dir.create("reports/Fig3/065-enrichr_to_circos", showWarnings = FALSE)
for (input_event in events) {
    data <- df_go %>%
        filter(event == input_event) %>%
        select(Term, from = ko_symbol) %>%
        mutate(to = str_remove(Term, " \\(GO:.*")) %>%
        select(-Term)

    top_n <- 10
    data_filtered <-
        data %>%
        count(to) %>%
        arrange(desc(n)) %>%
        slice_head(n = top_n) %>%
        inner_join(data, by = "to") %>%
        mutate(value = 1) %>%
        distinct() %>%
        select(from, to, value)

    jpeg(file = paste0("reports/Fig3/065-enrichr_to_circos/", input_event, ".jpg"), width = 1500, height = 1500, units = "px", res = 300)
    par(cex = 0.5)
    chordDiagram(data_filtered, transparency = 0.5, big.gap = 30)
    title(input_event)
    dev.off()

    # pdf(file = paste0("reports/Fig3/065-enrichr_to_circos/", input_event, ".pdf"), width = 1500, height = 1500)
    # par(cex = 0.5)
    # chordDiagram(data_filtered, transparency = 0.5, big.gap = 30)
    # title(input_event)
    # dev.off()

    data_filtered %>%
        select(top_10 = to) %>%
        distinct() %>%
        write_csv(paste0("reports/Fig3/065-enrichr_to_circos/", input_event, ".csv"))
}

# p <- wrap_elements(full = ~ chordDiagram(data_filtered))
# wrap_plots(plot_list)

# 全ての組み合わせを作成
all_combinations <- expand_grid(
    Var1 = unique(data_filtered$Var1),
    feature = unique(data_filtered$feature)
)

# 組み合わせにFreq列を追加
result <- data_filtered %>%
    full_join(all_combinations, by = c("Var1", "feature"), keep = TRUE) %>%
    mutate(Freq = if_else(is.na(Var1.x), 0, 1)) %>%
    select(feature = feature.y, Freq) %>%
    mutate(Var1 = "KO")

# 結果を表示
print(result)

chordDiagram(result, transparency = 0.5)

# Rコードでデータフレームを作成
data <- data.frame(
  gene = c("HSD17B12", "C19ORF40", "TRAF1", "ERCC1", "ITGA3", "CYCS", 
           "PRKAR2A", "XIAP", "GLUL", "POLM", "XRCC2", "XRCC5", "POLH", "MRE11A"),
  GOterm1 = c(0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0),
  GOterm2 = c(1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
  GOterm3 = c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0),
  GOterm4 = c(0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1),
  GOterm5 = c(0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0)
)

# 作成したデータフレームを表示
print(data)

chordDiagram(data, transparency = 0.5)

# patients <- c(rep("patient1",20), rep("patient2",10))
# cell.types <- c(rep("cell1",12), rep("cell2",8),rep("cell1",6), rep("cell2",4))
# features <- c(paste("feature",1:12,sep="_"), paste("feature",9:16,sep="_"), paste("feature",c(1,2,9,10,17,18),sep="_"), paste("feature",c(1,18,19,20),sep="_"))
# dat <- data.frame(patient=patients, cell.type=cell.types, feature=features)
# dat
# dat <- with(dat, table(paste(patient,cell.type,sep='|'), feature))
# dat

# as.data.frame(dat)
# chordDiagram(as.data.frame(dat), transparency = 0.5)



# 入力データ
data <- tibble(
  Var1 = c("A", "A", "A", "B", "B", "C"),
  feature = c("foo", "bar", "hoge", "foo", "bar", "foo")
)

# 全ての組み合わせを作成
all_combinations <- expand_grid(
  Var1 = unique(data$Var1),
  feature = unique(data$feature)
)

# 組み合わせにFreq列を追加
result <- data %>%
  full_join(all_combinations, by = c("Var1", "feature"), keep = TRUE) %>%
  mutate(Freq = if_else(is.na(Var1.x), 0, 1)) %>%
  select(feature = feature.y, Freq) %>%
  mutate(Var1 = "KO")

# 結果を表示
print(result)

chordDiagram(result, transparency = 0.5)
