###############################################################################
# Install Packages
###############################################################################

library(tidyverse)
library(ggrepel)

###############################################################################
# Load Data
###############################################################################

df_counts <- read_tsv("data/counts/featurecounts_gene_name.tsv.gz")
marker_genes <- read_tsv("data/markers.tsv")

#------------------------------------------------------------------------------
# TPM Calculation
#------------------------------------------------------------------------------

# Function to calculate TPM for a single column (sample)
calculate_tpm <- function(counts, lengths) {
    rpk <- counts / (lengths / 1000)  # Reads per kilobase (RPK)
    per_million_scaling_factor <- sum(rpk) / 1e6  # Scaling factor
    tpm <- rpk / per_million_scaling_factor  # TPM
    return(tpm)
}

# Apply the TPM calculation for each sample (i.e., each column of counts)
df_tpm <- df_counts %>%
    mutate(across(starts_with("MIERU_CH"), ~ calculate_tpm(.x, Length), .names = "{col}")) %>%
    select(!Length)

df_expression <-
    df_tpm %>%
    pivot_longer(cols = -Geneid, names_to = "sample", values_to = "expression") %>%
    as_tibble()

df_mean_expression <-
    df_expression %>%
    mutate(expression = expression + 1) %>%
    group_by(sample, Geneid) %>%
    summarize(logmean = log2(mean(expression, na.rm = TRUE))) %>%
    ungroup()

df_markers <-
    marker_genes %>%
    pivot_longer(everything(), names_to = "tissue", values_to = "Geneid")

df_plot <-
    left_join(df_mean_expression, df_markers, by = "Geneid") %>%
    mutate(tissue = if_else(!is.na(tissue), tissue, "other")) %>%
    mutate(label = if_else(tissue == "other", as.character(NA), Geneid)) %>%
    mutate(colors = case_when(
        tissue == "other" ~ "#DDDDDD",
        tissue == "exoderm" ~ "#ff4b00",
        tissue == "mesoderm" ~ "#005aff",
        tissue == "endoderm" ~ "#03af7a",
        tissue == "visceral endoderm" ~ "#8B8000"
    ))


df_plot$sample <-
    df_plot$sample %>%
    str_replace("MIERU_CH_1050$", "TC (\u00D715)") %>%
    str_replace("MIERU_CH_", "TC (\u00D71) ") %>%
    str_replace("MIERU_Homo_", "NM (\u00D71) ")

levels <- c("TC (\u00D71) 1", "TC (\u00D71) 2", "TC (\u00D71) 3", "TC (\u00D71) 4",
    "NM (\u00D71) 1", "NM (\u00D71) 2", "NM (\u00D71) 3", "NM (\u00D71) 4",
    "TC (\u00D715)")

df_plot$sample <- factor(df_plot$sample, levels = levels)

###############################################################################
# Violin plot
###############################################################################

options(ggrepel.max.overlaps = Inf)

g_violin <-
    ggplot(df_plot, aes(x = "", y = logmean, label = label, fill = colors)) +
    geom_violin(fill = "white") +
    geom_label_repel(
        color = "white",
        segment.colour = "black",
        point.padding = NA,
        box.padding = 0.5
    ) +
    scale_fill_identity(colors) +
    labs(x = "", y = "log2(TPM + 1)") +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    facet_wrap(~sample, scales = "free_y", ncol = 4)

ggsave("reports/violinplot_markers.jpg", g_violin, width = 15, height = 12, dpi = 600)
ggsave("reports/violinplot_markers.pdf", g_violin, width = 15, height = 12)
