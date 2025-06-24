library(tidyverse)

df_all_symbol_go <- read_tsv("data/Fig5/mgi_symbol_go.txt", col_names = c("symbol", "go")) %>%
    mutate(symbol = toupper(symbol))

df_complex_symbol_go <- read_csv("reports/Fig5/corum_go_symbol.csv")

#########################

df_go <- read_csv("reports/Fig4/enrichr_go_pathways.csv")

data <- df_go %>%
    select(Term, from = ko_symbol) %>%
    mutate(to = str_remove(Term, " \\(GO:.*"))

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

data_go_symbol <-
    data_filtered %>%
    select(to) %>%
    inner_join(data, by = "to", relationship = "many-to-many") %>%
    select(Term) %>%
    distinct() %>%
    inner_join(df_go, by = "Term") %>%
    select(Term, Genes) %>%
    group_by(Term) %>%
    summarise(Genes = paste(Genes, collapse = ", ")) %>%
    separate_longer_delim(Genes, delim = ";") %>%
    separate_longer_delim(Genes, delim = ", ") %>%
    distinct() %>%
    mutate(go = stringr::str_extract(Term, "\\(GO:\\d+\\)") %>% stringr::str_remove_all("[()]")) %>%
    select(Term, go, symbol = Genes)


results_fisher <- tibble()
# input_go <- "GO:0003723" # RNA binding
go_list <- data_go_symbol %>% pull(go) %>% unique()
for (input_go in go_list) {
    term <- data_go_symbol %>% filter(go == input_go) %>% pull(Term) %>% unique()
    data_symbols <- data_go_symbol %>% filter(go == input_go) %>% pull(symbol)
    all_symbols <- df_all_symbol_go %>% filter(go == input_go) %>% pull(symbol)
    complex_symbols <- df_complex_symbol_go %>% filter(go == input_go) %>% pull(symbol) %>% unique()

    overlap_data_complex <- data_symbols %in% complex_symbols
    overlap_all_complex <- all_symbols %in% complex_symbols

    a <- sum(overlap_data_complex)
    b <- sum(!overlap_data_complex)
    c <-sum (overlap_all_complex)
    d <- sum(!overlap_all_complex)
    vx <- matrix(c(a,b,c,d),nrow=2,byrow=T)
    result <- fisher.test(vx)
    sig <- ifelse(result$p.value < 0.05, "YES", "NO")
    results_fisher <- bind_rows(results_fisher, tibble(term = term, significance = sig, p_value = result$p.value, odds_ratio = result$estimate, se_forming_complex = a, se_not_forming_complex = b, all_go_genes_forming_complex = c, all_go_genes_not_forming_complex = d))
}

write_csv(results_fisher, "reports/Fig5/fisher_corum.csv")
# p <- wrap_elements(full = ~ chordDiagram(data_filtered))
# wrap_plots(plot_list)

# # 全ての組み合わせを作成
# all_combinations <- expand_grid(
#     Var1 = unique(data_filtered$Var1),
#     feature = unique(data_filtered$feature)
# )

# # 組み合わせにFreq列を追加
# result <- data_filtered %>%
#     full_join(all_combinations, by = c("Var1", "feature"), keep = TRUE) %>%
#     mutate(Freq = if_else(is.na(Var1.x), 0, 1)) %>%
#     select(feature = feature.y, Freq) %>%
#     mutate(Var1 = "KO")

# # 結果を表示
# print(result)

# chordDiagram(result, transparency = 0.5)

# # Rコードでデータフレームを作成
# data <- data.frame(
#   gene = c("HSD17B12", "C19ORF40", "TRAF1", "ERCC1", "ITGA3", "CYCS", 
#            "PRKAR2A", "XIAP", "GLUL", "POLM", "XRCC2", "XRCC5", "POLH", "MRE11A"),
#   GOterm1 = c(0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0),
#   GOterm2 = c(1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
#   GOterm3 = c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0),
#   GOterm4 = c(0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1),
#   GOterm5 = c(0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0)
# )

# # 作成したデータフレームを表示
# print(data)

# chordDiagram(data, transparency = 0.5)

# patients <- c(rep("patient1",20), rep("patient2",10))
# cell.types <- c(rep("cell1",12), rep("cell2",8),rep("cell1",6), rep("cell2",4))
# features <- c(paste("feature",1:12,sep="_"), paste("feature",9:16,sep="_"), paste("feature",c(1,2,9,10,17,18),sep="_"), paste("feature",c(1,18,19,20),sep="_"))
# dat <- data.frame(patient=patients, cell.type=cell.types, feature=features)
# dat
# dat <- with(dat, table(paste(patient,cell.type,sep='|'), feature))
# dat

# as.data.frame(dat)
# chordDiagram(as.data.frame(dat), transparency = 0.5)



# # 入力データ
# data <- tibble(
#   Var1 = c("A", "A", "A", "B", "B", "C"),
#   feature = c("foo", "bar", "hoge", "foo", "bar", "foo")
# )

# # 全ての組み合わせを作成
# all_combinations <- expand_grid(
#   Var1 = unique(data$Var1),
#   feature = unique(data$feature)
# )

# # 組み合わせにFreq列を追加
# result <- data %>%
#   full_join(all_combinations, by = c("Var1", "feature"), keep = TRUE) %>%
#   mutate(Freq = if_else(is.na(Var1.x), 0, 1)) %>%
#   select(feature = feature.y, Freq) %>%
#   mutate(Var1 = "KO")

# # 結果を表示
# print(result)

# chordDiagram(result, transparency = 0.5)
