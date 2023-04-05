WGCNA on the GPC subset
================
Nguyen Huynh
12/20/2022

``` r
library(WGCNA)
options(stringsAsFactors = FALSE)
```

``` r
dat <- read.csv("output/Denoised_Competition_matrix.csv", row.names = 1)
# Filter for GPC 
meta <- read.csv("output/DCA_cellInfo.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
all.equal(colnames(dat) ,row.names(meta))
dat <- dat[, colnames(dat) %in% row.names(meta)[which(meta$CellType == "GPC")]]
# Use all data
filt <- log2(as.matrix(dat)+1e-6)
filt <- t(filt)
```

``` r
# Pick soft thresholding power
powers <- c(seq(from=1, to=41, by = 2))
sft <- pickSoftThreshold(filt, dataIsExpr = TRUE, powerVector = powers, corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.1), networkType = "signed", verbose = TRUE)
```

``` r
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.9,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

``` r
softPower = 19
#turn expression into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(filt,networkType = "signed", TOMType = "signed", power = softPower, corType="bicor", maxPOutliers = 0.1)

# save.image("output/WGCNA.RData")
colnames(TOM) <- colnames(filt)
row.names(TOM) <- colnames(filt)
dissTOM <- 1- TOM
geneTree <- flashClust::flashClust(as.dist(dissTOM), method = "average")
```

``` r
minModuleSize <- 20
dynamicMod <- cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid", minClusterSize = minModuleSize, deepSplit = 4, pamRespectsDendro = FALSE, cutHeight = 0.99, respectSmallClusters = TRUE)
table(dynamicMod)
dynamicColor <- labels2colors(dynamicMod)

MEs <- moduleEigengenes(filt, dynamicColor)
MEG <- MEs$eigengenes

names(dynamicColor) <- colnames(filt)
# save(filt, dynamicColor, MEs, geneTree, file="output/WGCNA_output.RData")
# saveRDS(MEG, "output/MEG.RDS")
```

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /gpfs/fs2/scratch/sgoldman_lab/.conda/envs/nh.v8.4/lib/R/lib/libRblas.so
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
    ##  [1] compiler_3.5.1  magrittr_1.5    fastmap_1.1.0   tools_3.5.1    
    ##  [5] htmltools_0.5.2 yaml_2.2.0      stringi_1.2.4   rmarkdown_2.1  
    ##  [9] knitr_1.27      stringr_1.3.1   xfun_0.12       digest_0.6.25  
    ## [13] rlang_0.4.11    evaluate_0.14
