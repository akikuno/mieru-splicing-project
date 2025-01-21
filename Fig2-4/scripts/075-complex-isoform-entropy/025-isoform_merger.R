library(tidyverse)
library(janitor)

directory_path <- "data/rsem/bam"

df_file_path <- list.files(directory_path, full.names = TRUE) %>%
    as.data.frame() %>%
    set_names("file_path") %>%
    filter(str_detect(file_path, "isoforms.results")) %>%
    filter(str_detect(file_path, "Cd2bp2|MIERU|Qk|Rbm24|Spen|Trim71|Ubr5|Wt1|Ybx1"))

df_isoform <- tibble()
for (file_path in df_file_path$file_path) {
    sample_name <- basename(file_path) %>% str_remove("\\..*$")

    file_data <- read_tsv(file_path, show_col_types = FALSE) %>%
        clean_names() %>%
        mutate(sample = sample_name)

    df_isoform <- bind_rows(df_isoform, file_data)
}

# print(df_isoform)

dir.create("data/Fig7", showWarnings = FALSE)
df_isoform %>%
    mutate(gene_symbol = str_remove(gene_id, "^.*_")) %>%
    mutate(transcript_symbol = str_remove(transcript_id, "^.*_")) %>%
    select(sample, gene_symbol, transcript_symbol, tpm) %>%
    write_tsv("data/Fig7/tpm_isoforms_all.tsv.gz")
