Network Construction from SCENIC, WGCNA, and IPA
================
John Mariani
12/6/2022

``` r
library(data.table)
library(dplyr)
library(tidyr)
```

# Load Network edges from SCENIC and WGCNA analysis

``` r
networkEdges <- read.delim("output/regulons_filt_module.txt")
networkNodes <- data.frame(Node = unique(c(networkEdges$TF, networkEdges$Target)))
networkNodes$Node <- as.character(networkNodes$Node)
```

# Filter for genes DE in any comparison in either model and export for IPA

``` r
allDE <- read.delim("DE/allDE.txt")

networkEdges_de <- networkEdges[networkEdges$Target %in% allDE$Gene,]
networkNodes_de <- data.frame(Node = unique(c(networkEdges_de$TF, networkEdges_de$Target)))

write.table(networkNodes_de, "data_for_import/network_de_for_ipa.txt", sep = "\t", quote = F, row.names = F)
```

# Import IPA information and export network edges

``` r
files <- "data_for_import/ipa_blue_regulon_de.txt"
compNames <- c("TFregulons")
upstreamIPA <- fread(files, skip = "Upstream Regulators")
```

    ## Warning in fread(files, skip = "Upstream Regulators"): Stopped early on line
    ## 1619. Expected 13 fields but found 0. Consider fill=TRUE and comment.char=.
    ## First discarded non-empty line: <<Causal Networks for My Projects->HD_scRescue-
    ## >network_de_for_ipa - 2023-01-03 08:29 PM>>

``` r
upstreamIPA <- upstreamIPA[,c(3,10,11)]
names(upstreamIPA) <- c("Pathway", "pVal", "Genes")
upstreamIPA$Pathway <- paste0(upstreamIPA$Pathway, " Signaling")
canonicalIPA <- fread(files, skip = "Canonical Pathways")
```

    ## Warning in fread(files, skip = "Canonical Pathways"): Detected 5 column names
    ## but the data has 6 columns (i.e. invalid file). Added 1 extra default column
    ## name for the first column which is guessed to be row names or an index. Use
    ## setnames() afterwards if this guess is not correct, or fix the file write
    ## command that created the file to create a valid file.

    ## Warning in fread(files, skip = "Canonical Pathways"): Stopped early on line
    ## 488. Expected 6 fields but found 0. Consider fill=TRUE and comment.char=. First
    ## discarded non-empty line: <<Upstream Regulators for My Projects->HD_scRescue-
    ## >network_de_for_ipa - 2023-01-03 08:29 PM>>

``` r
canonicalIPA <- canonicalIPA[,c(1,2,5)]
names(canonicalIPA)<- c("Pathway", "pVal", "Genes")
functionalIPA <- fread(files, skip = "Diseases and Bio Functions")
```

    ## Warning in fread(files, skip = "Diseases and Bio Functions"): Stopped early on
    ## line 3215. Expected 11 fields but found 0. Consider fill=TRUE and comment.char=.
    ## First discarded non-empty line: <<Tox Functions for My Projects->HD_scRescue-
    ## >network_de_for_ipa - 2023-01-03 08:29 PM>>

``` r
functionalIPA <- functionalIPA[,c(3,4,9)]
names(functionalIPA) <- c("Pathway", "pVal", "Genes")
IPA <- rbind(upstreamIPA, canonicalIPA)
IPA <- rbind(IPA, functionalIPA)

edges <- IPA %>% 
  mutate(Gene = strsplit(as.character(Genes), ",")) %>% 
  unnest(Gene) %>% .[,-3]

edges <- edges[edges$Gene %in% networkNodes$Node,]

goCats <- c("EIF2 Signaling", "LARP1 Signaling", "MLXIPL Signaling", "mTOR Signaling", "MYC Signaling", "MYCN Signaling", "Oxidative Phosphorylation", "RICTOR Signaling", "Synthesis of protein", "TP53 Signaling", "YAP1 Signaling", "CEBPZ Signaling",
            "CTCF Signaling", "E2F1 Signaling", "ETV4 Signaling", "NFYB Signaling")

edges <- edges[edges$Pathway %in% goCats,]

names(edges)
```

    ## [1] "Pathway" "pVal"    "Gene"

``` r
networkEdges_de$Type <- "SCENIC"
names(networkEdges_de)[1] <- "Source"

edges$Pathway <- gsub(x = edges$Pathway, pattern = " Signaling", replacement = "")
edges$pVal <- NULL
names(edges) <- c("Source", "Target")
edges$Type <- "IPA"

scenicEdges <- networkEdges_de
scenicEdges$Weight <- NULL


networkEdges_Final <- rbind(scenicEdges, edges)


write.table(networkEdges_Final, "data_for_import/networkEdges_de.txt", quote = F, row.names = F, sep = "\t")
```

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] tidyr_1.1.3       dplyr_1.0.7       data.table_1.14.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.13   knitr_1.33        magrittr_2.0.1    tidyselect_1.1.1 
    ##  [5] R6_2.5.0          rlang_1.0.6       fansi_0.5.0       stringr_1.4.0    
    ##  [9] tools_4.1.0       xfun_0.24         utf8_1.2.2        DBI_1.1.1        
    ## [13] cli_3.4.1         htmltools_0.5.1.1 ellipsis_0.3.2    assertthat_0.2.1 
    ## [17] yaml_2.2.1        digest_0.6.27     tibble_3.1.3      lifecycle_1.0.3  
    ## [21] crayon_1.4.1      purrr_0.3.4       vctrs_0.5.1       glue_1.4.2       
    ## [25] evaluate_0.14     rmarkdown_2.9     stringi_1.7.3     compiler_4.1.0   
    ## [29] pillar_1.6.2      generics_0.1.0    pkgconfig_2.0.3
