# correlation network

# install.packages("igraph")
library(igraph)

# ! Mandatorial
infile <- "toy_correlations_big.tsv"
# change with your own input file
# a TSV file with a header and (at least) 4 columns, in this order:
# gene1, gene2, Pearson correlation coefficient (r), FDR (already corrected p-value)
# column names do not matter: the first 4 columns are used, any other column is ignored
# pairs can appear in both directions (A-B and B-A): duplicates are removed
# example:
# geneA   geneB   r_pearson   FDR
# mod1_g01   mod1_g02   0.8156   3.165e-04
# mod1_g02   mod1_g01   0.8156   3.165e-04

# ! Optional
r_min   <- 0.5      # keep edges with |r| >= r_min (between 0 and 1)
fdr_max <- 0.001    # keep edges with FDR < fdr_max (> 0 and <= 0.05)
sign    <- "both"   # "both" (|r|), "positive" (r only) or "negative" (-r only)
labels  <- FALSE    # show gene names? (FALSE with many nodes is a good idea)

# Read the TSV, rename the first 4 columns and print the mapping (original >> new)
# so the user can check which columns are being used
load_edges <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  new <- c("gene1", "gene2", "r", "fdr")
  cat(sprintf("%s >> %s\n", names(df)[1:4], new), sep = "")
  if (ncol(df) > 4) cat(sprintf("ignored columns: %d\n", ncol(df) - 4))
  df <- df[, 1:4]
  names(df) <- new
  # remove self loops (geneA - geneA)
  df <- df[df$gene1 != df$gene2, ]
  # write every pair in the same (alphabetical) order, so A-B and B-A become identical
  a <- df$gene1
  b <- df$gene2
  m <- a < b
  df$gene1 <- ifelse(m, a, b)
  df$gene2 <- ifelse(m, b, a)
  # remove the duplicated pairs
  df[!duplicated(df[, c("gene1", "gene2")]), ]
}

# Keep only the edges above the r threshold and below the FDR threshold
filter_edges <- function(df, r_min, fdr_max, sign = c("both", "positive", "negative")) {
  sign <- match.arg(sign)
  # stop with a clear message if the thresholds are out of range
  if (r_min < 0 || r_min > 1) stop("r_min must be between 0 and 1")
  if (fdr_max <= 0 || fdr_max > 0.05) stop("fdr_max must be > 0 and <= 0.05")
  keep_r <- switch(sign,
                   both     = abs(df$r) >= r_min,
                   positive = df$r >= r_min,
                   negative = df$r <= -r_min)
  # which() drops rows with NA values instead of returning empty NA rows
  df[which(keep_r & df$fdr < fdr_max), ]
}

# Build an undirected graph: the first 2 columns are the nodes,
# the other columns (r, fdr) become edge attributes (E(g)$r, E(g)$fdr)
build_graph <- function(edges) {
  graph_from_data_frame(edges, directed = FALSE)
}

# Draw the graph in the plot window
# layout: Kamada-Kawai, nodes closer = fewer steps between them (r is not used)
# separate networks are arranged in a circle, each one inside a shaded area
plot_graph <- function(g, labels = TRUE) {
  edge_w    <- 1.5      # edge width
  node_size <- 2        # node size (keep it about proportional to edge_w)
  # one group of nodes per connected component (= per separate network)
  groups <- split(seq_len(vcount(g)), components(g)$membership)
  par(mar = c(1, 1, 1, 1), family = "Helvetica Neue")
  plot(g,
       layout              = layout_with_kk(g),
       mark.groups         = groups,
       mark.col            = adjustcolor("lightblue", alpha.f = 0.1),   # semi-transparent area
       mark.border         = NA,
       mark.expand         = 8,                                         # margin around the nodes
       vertex.size         = node_size,
       vertex.color        = "white",
       vertex.frame.color  = "grey30",
       vertex.label        = if (labels) V(g)$name else NA,
       vertex.label.family = "Helvetica Neue",                          # font family
       vertex.label.cex    = 1,                                         # font size (1 = 12 pt)
       vertex.label.color  = "black",
       vertex.label.dist   = 1.1,                                       # label distance from the node
       vertex.label.degree = -pi/2,                                     # label above the node
       edge.color          = ifelse(E(g)$r > 0, "steelblue", "red"),    # positive = blue, negative = red
       edge.width          = edge_w)
}

# Run: load -> filter -> graph -> plot
df    <- load_edges(infile)
edges <- filter_edges(df, r_min, fdr_max, sign)
g     <- build_graph(edges)
plot_graph(g, labels)
