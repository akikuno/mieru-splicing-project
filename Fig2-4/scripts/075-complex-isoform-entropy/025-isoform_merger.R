library(tidyverse)
library(janitor)

directory_path <- "data/rsem/bam"

df <- list.files(directory_path, full.names = TRUE) %>%
    as_tibble() %>%
    filter(str_detect(value, "isoforms.results")) %>%
    map_dfr(., function(file_path) {
        sample_name <- basename(file_path) %>% str_remove("\\..*$")
        print(sample_name)
        read_csv(file_path) %>%
            clean_names() %>%
            mutate(sample = sample_name)
    })


library(tidyverse)
library(janitor)

directory_path <- "data/rsem/bam"

df_file_path <- list.files(directory_path, full.names = TRUE) %>%
    as_tibble() %>%
    set_names("file_path") %>%
    filter(str_detect(file_path, "isoforms.results")) 
    
df_isoform <- tibble()
for (file_path in df_file_path$file_path) {
    sample_name <- basename(file_path) %>% str_remove("\\..*$")
    
    file_data <- read_tsv(file_path, show_col_types = FALSE) %>%
        clean_names() %>%
        mutate(sample = sample_name)
    
    df_isoform <- bind_rows(df_isoform, file_data)
}

print(df_isoform)

dir.create("data/Fig7",)
df_isoform %>%
    mutate(gene_symbol = str_remove(gene_id, "^.*_")) %>%
    mutate(transcript_symbol = str_remove(transcript_id, "^.*_")) %>%
    select(sample, gene_symbol, transcript_symbol, tpm) %>%
    write_tsv("data/Fig7/tpm_isoforms_all.tsv")

