Export data to use as input for pySCENIC
================
Nguyen Huynh
12/14/2022

``` r
library(Seurat)
object <- readRDS("RDS/mergedUpdated.rds")
scenic <- object[["RNA"]]@counts
dim(scenic)
# Filter to retain only genes that are expressed more than minimal counts per gene threshold
minCountsPerGene <- 3*0.01*ncol(scenic)
sum <- apply(scenic, 1, sum)
keep <- names(sum)[which(sum > minCountsPerGene)]
# Filter to retain only genes that are expressed in more than the minimal number of cells
minSamples <- ncol(scenic)*0.01
zeros <- apply(scenic, 1, function(x) sum(x>0))
keep2 <- names(zeros)[which(zeros > minSamples)]
keep <- keep[keep %in% keep2]
scenic <- as.matrix(scenic[keep,])
# write.csv(scenic, "output/Competition_forScenic.csv", quote = FALSE)
```

``` r
## keep genes with total counts > 1
total.counts <- apply(object[["RNA"]]@counts, 1, sum)
keep <- names(total.counts)[which(total.counts >= 1)]
## Write out for WGCNA 
# write.csv(object[["RNA"]]@counts[keep, ], "output/Competition_forWGCNA.csv", quote=FALSE)
```

``` r
## Genes and cell info for DCA
cells <- object@meta.data
cells$X <- row.names(cells)
# write.csv(cells, "output/DCA_cellInfo.csv", quote = FALSE)
genes <- data.frame(row.names = keep, Gene = keep)
# write.csv(genes, "output/DCA_geneInfo.csv", quote = FALSE)
```

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)
    ## 
    ## Matrix products: default
    ## BLAS:   /gpfs/fs1/sfw2/r/4.1.1/b1/lib64/R/lib/libRblas.so
    ## LAPACK: /gpfs/fs1/sfw2/r/4.1.1/b1/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.1.1  magrittr_2.0.2  fastmap_1.1.0   cli_3.3.0      
    ##  [5] tools_4.1.1     htmltools_0.5.2 rstudioapi_0.13 yaml_2.3.4     
    ##  [9] stringi_1.7.6   rmarkdown_2.11  knitr_1.37      stringr_1.4.0  
    ## [13] xfun_0.29       digest_0.6.29   rlang_1.0.1     evaluate_0.15
