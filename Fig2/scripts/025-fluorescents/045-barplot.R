library(tidyverse)

df <- read_csv("data/counts/read_counts_fluorescents.csv",
    col_names = c("sample", "total_reads", "fluorescentes", "reads", "% of reads")
    )

df <- df %>%
    mutate(colors = case_when(
        fluorescentes == "EGFP" ~ "#03af7a",
        fluorescentes == "tdTomato" ~ "#ff4b00",
        fluorescentes == "TagBFP" ~ "#005aff"
    ))

df_plot <- df %>%
    group_by(sample) %>%
    mutate(total_reads_fluorescents = sum(reads)) %>%
    mutate(`% of reads in fluorescents` = reads / total_reads_fluorescents * 100) %>%
    ungroup()

df_plot$sample <-
    df_plot$sample %>%
    str_replace("MIERU_CH_1050$", "TC (\u00D715)") %>%
    str_replace("MIERU_CH_", "TC (\u00D71) ") %>%
    str_replace("MIERU_Homo_", "NM (\u00D71) ")

levels <- c("TC (\u00D71) 1", "TC (\u00D71) 2", "TC (\u00D71) 3", "TC (\u00D71) 4",
    "NM (\u00D71) 1", "NM (\u00D71) 2", "NM (\u00D71) 3", "NM (\u00D71) 4",
    "TC (\u00D715)")

df_plot$sample <- factor(df_plot$sample, levels = levels)

g_reads_number <- ggplot(df_plot, aes(x = fluorescentes, y = reads, fill = colors)) +
    geom_col() +
    theme_bw() +
    labs(x = "", y = "Number of reads") +
    scale_fill_identity(colors) +
    facet_wrap(~sample, scales = "free_y", ncol = 4)

ggsave("reports/fluorescents_read_counts.jpg", g_reads_number, width = 10, height = 5, units = "in", dpi = 600)
ggsave("reports/fluorescents_read_counts.pdf", g_reads_number, width = 10, height = 5, units = "in")

g_reads_percentage <- ggplot(df_plot, aes(x = fluorescentes, y = `% of reads in fluorescents`, fill = colors)) +
    geom_col() +
    theme_bw() +
    labs(x = "", y = "Percentage of reads on fluorescents") +
    scale_fill_identity(colors) +
    ylim(0, 100) +
    facet_wrap(~sample, scales = "free_y", ncol = 4)

ggsave("reports/fluorescents_read_percentage.jpg", g_reads_percentage, width = 10, height = 5, units = "in", dpi = 600)
ggsave("reports/fluorescents_read_percentage.pdf", g_reads_percentage, width = 10, height = 5, units = "in")


g_percentage_total_reads <- ggplot(df, aes(x = fluorescentes, y = `% of reads`, fill = colors)) +
    geom_col() +
    theme_bw() +
    labs(x = "", y = "Number of reads on fluorescents") +
    scale_fill_identity(colors) +
    facet_wrap(~sample, scales = "free_y", ncol = 4)

ggsave("reports/fluorescents_percentage_total_reads.jpg", g_percentage_total_reads, width = 10, height = 5, units = "in", dpi = 600)
ggsave("reports/fluorescents_percentage_total_reads.pdf", g_percentage_total_reads, width = 10, height = 5, units = "in")
