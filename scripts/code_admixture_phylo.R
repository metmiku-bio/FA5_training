library(ggtree)

## visualisation of phylogenetic tree using ggtree
data_dir = '/Users/pmonsieurs/programming/ext_FA5_ethiopia_2026/data/'
tree <- read.tree(paste0(source_dir, "combined.filtered.min4.fasta.treefile"))
tree <- read.tree(paste0(source_dir, "combined.filtered.min4.no_outgroup.fasta.treefile"))
ggtree(tree) +
  geom_tiplab() +
  theme_tree2()