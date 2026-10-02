# AMAS summary outlier check

# install.packages(c("ggplot2", "ggrepel", "tidyr"))
library(ggplot2)
library(ggrepel)
library(tidyr)

# ! Mandatorial
summary_file <- "concatenated.out-seq-summary.txt"
# change with your own input file
# this file is produce by AMAS.py (Borowiec, 2016, doi: 10.7717/peerj.1660).
# run: AMAS.py summary -f fasta -d dna -i concatenated.fasta --by-taxon

# ! Mandatorial
columns_to_plot <- c("Missing_percent")
# which numeric columns from summary_file to check (names must match exactly)

# ! Optional
groups_file <- "groups.csv"
# groups.csv: a CSV file containing the group information for your species. The species names in the 'id' column must match the Taxon_name column in summary_file. The 'group' column indicates the group to which each species belongs ('color' is unused here, kept only so the same file also works with molPCA.R).
# example of groups.csv:
# group,color,id
# Group1,red,Species1
# Group1,red,Species2

# ! Optional
groups_to_keep <- c("Agamidae", "Serpentes")
# only used if groups_file is provided: groups not listed here collapse into "Other"

# prettier facet titles for columns_to_plot; anything not listed here keeps its raw column name
column_labels <- c(
  Missing_percent = "Missing (%)",
  Undetermined_characters = "Undetermined characters (count)"
)

# separate filenames
output_file_groups <- "AMAS_summary_groups.png"
output_file_taxa <- "AMAS_summary_taxa.png"

# read summary file and pivot to long format
summary_df <- read.delim(summary_file)
long_df <- pivot_longer(summary_df, cols = all_of(columns_to_plot), names_to = "Statistic", values_to = "Value")
long_df$Statistic <- ifelse(long_df$Statistic %in% names(column_labels), column_labels[long_df$Statistic], long_df$Statistic)

# groups_file is optional: with it, compare groups via boxplots instead of flagging individual outliers
has_groups <- file.exists(groups_file)

if (has_groups) {
  groups_df <- read.csv(groups_file, strip.white = TRUE)
  species_to_group <- setNames(groups_df$group, groups_df$id)
  long_df$Group <- ifelse(long_df$Taxon_name %in% names(species_to_group), species_to_group[long_df$Taxon_name], "Unassigned")
  if (length(groups_to_keep) > 0) long_df$Group <- ifelse(long_df$Group %in% groups_to_keep, long_df$Group, "Other")

  # flag standard boxplot outliers (beyond 1.5*IQR) within each Group x Statistic
  long_df$Outlier <- as.logical(ave(long_df$Value, long_df$Group, long_df$Statistic, FUN = function(x) {
    q <- quantile(x, c(0.25, 0.75))
    iqr <- diff(q)
    x < q[1] - 1.5 * iqr | x > q[2] + 1.5 * iqr
  }))

  ggplot(long_df, aes(x = Group, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = subset(long_df, Outlier), color = "red", size = 2) +
    geom_text_repel(data = subset(long_df, Outlier), aes(label = Taxon_name), color = "red", size = 3) +
    facet_wrap(~ Statistic, scales = "free", strip.position = "bottom") +
    coord_flip() +
    theme_minimal() +
    theme(
    text = element_text(family = "Helvetica Neue"), # font family
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 10),  # group names and values
    plot.title = element_text(size = 10),  # title
    strip.text = element_text(size = 10)   # panel labels, e.g. Missing_percent
    ) +
    theme(panel.grid = element_blank(), axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"), strip.placement = "outside") +
    labs(x = NULL, y = NULL, title = "AMAS summary - by group")
  ggsave(output_file_groups, width = 4, height = 5, dpi = 300)
  #ggsave(sub("png$", "svg", output_file_groups), width = 8, height = 4)

} else {
  # outlier = more than 2 standard deviations above the mean, per column
  long_df$Outlier <- as.logical(ave(long_df$Value, long_df$Statistic, FUN = function(x) x > mean(x) + 2 * sd(x)))
  thresholds <- aggregate(Value ~ Statistic, long_df, function(x) mean(x) + 2 * sd(x))
  names(thresholds)[2] <- "cutoff"

  ggplot(long_df, aes(x = reorder(Taxon_name, Value), y = Value, color = Outlier)) +
    geom_hline(data = thresholds, aes(yintercept = cutoff), linetype = "dashed", color = "grey50") +
    geom_point(size = 2) +
    geom_text_repel(data = subset(long_df, Outlier), aes(label = Taxon_name), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
    facet_wrap(~ Statistic, scales = "free", strip.position = "bottom") +
    coord_flip() +
    theme_minimal() +
    theme(
    text = element_text(family = "Helvetica Neue"),
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 10),  # taxon names and values
    plot.title = element_text(size = 10),
    strip.text = element_text(size = 10)
    ) +
    theme(panel.grid = element_blank(), axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"), strip.placement = "outside", legend.position = "none") +
    labs(x = NULL, y = NULL, title = "AMAS summary - outliers")
  ggsave(output_file_taxa, width = 4, height = 5, dpi = 300)
  #ggsave(sub("png$", "svg", output_file_taxa), width = 8, height = 8)
}

