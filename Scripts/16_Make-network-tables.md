Make network tables for Cytoscape
================
Nguyen Huynh
01/03/2023

``` r
# Making of the network 
# Node from filtered with 3 TF counters 
node <- read.table("output/regulons_filt_shared_targets_node.txt", header = TRUE, stringsAsFactors = FALSE)
# Filtered targets 
filt_targ <- node$node[node$Type == "Target"]
# Edge 
edge_ipa <- read.delim("output/networkEdges_de.txt")
edge_ipa <- edge_ipa[edge_ipa$Target %in% filt_targ, ]
edge <- read.table("output/regulons_filt_shared_targets.txt", header = TRUE, stringsAsFactors = FALSE)

# IPA edges should contain all the nodes (filtered once more for deg only), plus IPA
a <- c(edge$TF, edge$Target) %>% unique()
b <- c(edge_ipa$Source, edge_ipa$Target) %>% unique()
b[!b %in% a]

# Use edge_ipa, make new node file by adding in the IPAs 
node_ipa <- node[node$node %in% b, ]
tmp <- data.frame(node = b[!b %in% a], 
                  ME = NA, 
                  age_strength = NA, 
                  Type = "IPA")
node_ipa <- rbind(node_ipa, tmp)
# write.table(node_ipa, "output/networkNodes_de_filt_NH.txt", quote = FALSE, sep = "\t", row.names = FALSE)
# write.table(edge_ipa, "output/networkEdges_de_filt_NH.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```

``` r
sessionInfo()
```
