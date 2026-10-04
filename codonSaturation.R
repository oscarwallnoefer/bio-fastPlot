# third-codon-position saturation plot

# install.packages(c("ape", "ggplot2"))
library(ape)
library(ggplot2)

# ! Mandatorial
fasta_files <- c("COX1.fasta", "ATP6.fasta")
# one or more aligned, in-frame, codon-aligned fasta files; more than one makes a multipanel plot, one gene per panel

# ! Optional
substitution_model <- "T92"
# Tamura 1992 correction, same as EMBOSS distmat's "Tamura" method; see ?dist.dna for alternatives

# ! Optional
positions_to_plot <- c("1st+2nd", "3rd")
# which codon position(s) to compute and overlay for comparison; remove entries you don't want

extract_positions <- function(alignment, position) {
  n <- ncol(alignment)
  switch(position,
    "All"     = alignment,
    "1st"     = alignment[, seq(1, n, by = 3)],
    "2nd"     = alignment[, seq(2, n, by = 3)],
    "3rd"     = alignment[, seq(3, n, by = 3)],
    "1st+2nd" = alignment[, sort(c(seq(1, n, by = 3), seq(2, n, by = 3)))]
  )
}

compute_saturation <- function(fasta_file) {
  alignment <- read.dna(fasta_file, format = "fasta")
  gene <- tools::file_path_sans_ext(basename(fasta_file))
  do.call(rbind, lapply(positions_to_plot, function(pos) {
    sub_aln <- extract_positions(alignment, pos)
    raw_dist <- dist.dna(sub_aln, model = "raw", pairwise.deletion = TRUE)
    corrected_dist <- dist.dna(sub_aln, model = substitution_model, pairwise.deletion = TRUE)
    df <- data.frame(Corrected = as.vector(corrected_dist), Observed = as.vector(raw_dist))
    df[!is.finite(df$Corrected) | !is.finite(df$Observed), ] <- NA
    df <- na.omit(df)
    df$Gene <- gene
    df$Position <- pos
    df
  }))
}

sat_df <- do.call(rbind, lapply(fasta_files, compute_saturation))

# same max on both axes, per gene, so the diagonal sits at a true 45 degrees
sat_df$Limit <- ave(pmax(sat_df$Corrected, sat_df$Observed), sat_df$Gene, FUN = max)

ggplot(sat_df, aes(x = Corrected, y = Observed)) +
  geom_blank(aes(x = Limit, y = Limit)) +
  geom_blank(aes(x = 0, y = 0)) +
  geom_abline(slope = 1, intercept = 0, color = "lightgrey", linewidth = 0.6) +
  geom_point(aes(color = Position), size = 2, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free") +
  theme_minimal() +
  theme(
  text = element_text(family = "Helvetica Neue"),
  axis.title = element_text(size = 10),
  axis.text  = element_text(size = 10),
  plot.title = element_text(size = 10),
  strip.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10)
  ) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"), aspect.ratio = 1) +
  labs(title = "Substitution saturation by codon position", x = paste0("Corrected (p) distance"), y = "Uncorrected (p) distance", color = "Position")

n_col <- ceiling(sqrt(length(unique(sat_df$Gene))))
n_row <- ceiling(length(unique(sat_df$Gene)) / n_col)
ggsave("saturation_3rd_pos.png", width = n_col * 3, height = n_row * 3, dpi = 300)
#ggsave("saturation_3rd_pos.svg", width = n_col * 4, height = n_row * 4)
