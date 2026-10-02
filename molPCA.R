# molecular PCA

# install.packages(c("adegenet", "ggplot2", "ggrepel"))
library(adegenet)
library(ggplot2)
library(ggrepel)

# ! Mandatorial
fasta_file <- "genes.fasta"
# change with your own input file
# genes.fasta: a fasta file containing the aligned sequences of your species

# ! Optional
groups_file <- "groups.csv"
# groups.csv: a CSV file containing the group information for your species. The species names in the 'id' column must match the species names in the fasta file. The 'group' column indicates the group to which each species belongs, and the 'color' column specifies the color to be used for that group in the PCA plot.
# example of groups.csv:
# group,color,id
# Group1,red,Species1
# Group1,red,Species2   

molecular_data <- fasta2genlight(fasta_file, snpOnly = TRUE)
print(molecular_data)
# 'nf = 3' retains the top 3 Principal Components for downstream analysis
pca_results <- glPca(molecular_data, nf = 3)
# Convert the PCA scores matrix into a readable data frame
pca_df <- as.data.frame(pca_results$scores)
# Extract species names from the row names and save them into a dedicated column
pca_df$Species <- rownames(pca_df)
# Calculate the exact percentage of variance explained by each PC
variance_explained <- (pca_results$eig / sum(pca_results$eig)) * 100

# groups_file columns: group, color, id
# groups_file is optional: no file, no coloring, just one black color
has_groups <- file.exists(groups_file)
if (has_groups) {
  groups_df <- read.csv(groups_file, strip.white = TRUE)
  species_to_group <- setNames(groups_df$group, groups_df$id)
  group_to_color <- setNames(groups_df$color, groups_df$group)[!duplicated(groups_df$group)]
  pca_df$Group <- ifelse(pca_df$Species %in% names(species_to_group), species_to_group[pca_df$Species], "Unassigned")
  group_to_color["Unassigned"] <- "black"
} else {
  pca_df$Group <- "All"
  group_to_color <- c(All = "black")
}

# Generates a clear scatter plot mapping the genetic distance between your species
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Species)) +
  geom_point(size = 2) +
  scale_color_manual(values = group_to_color, breaks = names(group_to_color)[names(group_to_color) != "Unassigned"]) +
  geom_text_repel(max.overlaps = Inf, show.legend = FALSE,
    aes(label = Species),
    box.padding = 0.2, #space around the text label
    point.padding = 0.3, #space around the point
    segment.color = "grey50",
    segment.size = 0.3) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black")) +
  labs(title = "", x = paste0("PC1 (", round(variance_explained[1], 1), "%)"), y = paste0("PC2 (", round(variance_explained[2], 1), "%)"))

