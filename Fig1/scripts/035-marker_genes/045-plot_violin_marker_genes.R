###############################################################################
# Install Packages
###############################################################################

options(repos = "http://cran.us.r-project.org")
options(readr.show_col_types = FALSE)
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, tidyr, tidyfast, dplyr, stringr, ggplot2, ggrepel)

###############################################################################
# Download and format the data
###############################################################################

tpm_matrix <- read_tsv("data/TPM_Matrix_RIAS.tsv")
marker_genes <- read_tsv("data/markers.tsv")

tbl_expression <-
    tpm_matrix %>%
    dt_pivot_longer(cols = -gene_sym, names_to = "sample", values_to = "expression") %>%
    as_tibble()

tbl_mean_expression <-
    tbl_expression %>%
    group_by(sample, gene_sym) %>%
    summarize(logmean = log2(mean(expression, na.rm = TRUE))) %>%
    # if mean value == 0, log transformation makes it as -Inf,
    # so match the minimum value for each group.
    mutate(logmean = if_else(logmean == -Inf, Inf, logmean)) %>%
    mutate(logmean = if_else(logmean == Inf, min(logmean), logmean)) %>%
    ungroup()

tbl_markers <-
    marker_genes %>%
    pivot_longer(everything(), names_to = "tissue", values_to = "gene_sym")

df_plot <-
    left_join(tbl_mean_expression, tbl_markers, by = "gene_sym") %>%
    mutate(tissue = if_else(!is.na(tissue), tissue, "other")) %>%
    mutate(label = if_else(tissue == "other", as.character(NA), gene_sym)) %>%
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
# df_plot$sample <- factor(df_plot$sample, levels = c(c("MIERU_CH_1", "MIERU_CH_2", "MIERU_CH_3", "MIERU_CH_4", "MIERU_Homo_1", "MIERU_Homo_2", "MIERU_Homo_3", "MIERU_Homo_4", "MIERU_CH_1050")))

# df_plot %>% filter(gene_sym == "Sox17")
# df_plot %>% filter(gene_sym == "Otx2")
# df_plot %>% filter(gene_sym == "T")


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
    labs(x = "", y = "log2 TPM") +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    facet_wrap(~sample, scales = "free_y", ncol = 4)

ggsave("reports/plot_violin_markers.png", g_violin, width = 15, height = 12, dpi = 600)
ggsave("reports/plot_violin_markers.pdf", g_violin, width = 15, height = 12)


# g_violin <- ggplot(df_plot) +
#     aes(
#         x = tissue, y = logmean, label = label,
#         fill = factor(colors)
#     ) +
#     geom_violin(position = position_dodge(6)) +
#     # geom_boxplot(
#     #     outlier.shape = NA,
#     #     position = position_dodge(1),
#     #     color = "#333333",
#     #     width = 0.2,
#     #     alpha = 0.2
#     # ) +
#     scale_fill_identity(colors) +
#     labs(x = "", y = "log2 TPM") +
#     theme_bw() +
#     facet_wrap(~sample, scales = "free_y", ncol = 4)

# ggsave("reports/plot_violin_tissues.png", g_violin, width = 15, height = 5, dpi = 600)
