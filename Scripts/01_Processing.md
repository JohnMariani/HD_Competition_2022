Processing of scRNA-Seq Data
================
John Mariani
12/10/22

\#\#Load Libraries

``` r
library(Seurat)
```

    ## Attaching SeuratObject

    ## Attaching sp

``` r
library(textTinyR)
```

    ## Loading required package: Matrix

``` r
library(SeuratObject)


options(future.globals.maxSize = 16000 * 1024^2)

`%not in%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
```

\#\#Read in counts and create list of Seurat Objects

``` r
sampleList <- list.files("Matrices/")
sampleList
```

    ##  [1] "G19_RescueCtrl_eGFP"                    
    ##  [2] "G19_RescueCtrl_eGFP_alone"              
    ##  [3] "G19_RescueCtrl_mCherry"                 
    ##  [4] "G19_RescueCtrl_mCherry_alone"           
    ##  [5] "Mixed_Rescue_mCherry_G19_1"             
    ##  [6] "Rescue_eGFP_G20_2"                      
    ##  [7] "Rescue_mCherry_G19_2"                   
    ##  [8] "Sample_rag1_Rescue_Ctrl_G19_mCherry_Neg"
    ##  [9] "Sample_rag1_Rescue_Ctrl_G19_mCherry_Pos"
    ## [10] "Sample_rag1_Rescue_Ctrl_G20_eGFP_Neg"   
    ## [11] "Sample_rag1_Rescue_Ctrl_G20_eGFP_Pos"   
    ## [12] "Sample_rag1_rescue_GFP_g20_1"           
    ## [13] "Sample_rag1_rescue_mCherry_g19_1"

``` r
# Read in raw matrices and throw out droplets with fewer than 250 umi's
raw <- sapply(sampleList, function(x) {print(x) ; Read10X(paste0("Matrices/",x,"/star2/hs/raw"))})
```

    ## [1] "G19_RescueCtrl_eGFP"
    ## [1] "G19_RescueCtrl_eGFP_alone"
    ## [1] "G19_RescueCtrl_mCherry"
    ## [1] "G19_RescueCtrl_mCherry_alone"
    ## [1] "Mixed_Rescue_mCherry_G19_1"
    ## [1] "Rescue_eGFP_G20_2"
    ## [1] "Rescue_mCherry_G19_2"
    ## [1] "Sample_rag1_Rescue_Ctrl_G19_mCherry_Neg"
    ## [1] "Sample_rag1_Rescue_Ctrl_G19_mCherry_Pos"
    ## [1] "Sample_rag1_Rescue_Ctrl_G20_eGFP_Neg"
    ## [1] "Sample_rag1_Rescue_Ctrl_G20_eGFP_Pos"
    ## [1] "Sample_rag1_rescue_GFP_g20_1"
    ## [1] "Sample_rag1_rescue_mCherry_g19_1"

``` r
raw <- sapply(raw, function(x) x[,sparse_Sums(x, rowSums = F) > 250])
sets <- length(raw)


#  Convert to list of seurat objects and filter for quality
ObjectsH <- sapply(c(1:sets), function(x) CreateSeuratObject(raw[[x]], project = sampleList[x]))
```

\#\#Filter for quality and merge

``` r
for (i in 1:sets) {
  ObjectsH[[i]] <- PercentageFeatureSet(ObjectsH[[i]], pattern = "^MT-", col.name = "percent.mt")
  ObjectsH[[i]] <- subset(x = ObjectsH[[i]], subset = nFeature_RNA > 250 & percent.mt < 15)
}



merged <- merge(ObjectsH[[1]], y = ObjectsH[2:length(ObjectsH)])
```

    ## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
    ## duplicated across objects provided. Renaming to enforce unique cell names.

``` r
dim(merged)
```

    ## [1] 38743 11636

``` r
metaMerged <- merged@meta.data


sampleTable <- read.csv("data_for_import/HDsampleTable.csv")

metaMerged$cellName <- row.names(metaMerged)
metaMerged <- merge(metaMerged, sampleTable, by.x = "orig.ident", by.y = "sample")
```

## Remove ambiguous cells

``` r
DefaultAssay(merged) <- "RNA"


metaMerged <- merge(metaMerged, FetchData(merged, c("EGFP", "mCherry"), slot = "counts"), by.x = "cellName", by.y = 0)


doubleCells <- metaMerged[metaMerged$EGFP>0 & metaMerged$mCherry >0,]


metaMerged <- metaMerged[metaMerged$cellName %not in% c(doubleCells$cellName),]


# Assign cells that were sorted on the appropriate fluorophore and have appropriate fluorophore expression to the correct sample.
# Keep cells in the negative captures for those that were singly engrafted

metaMerged$SampleName <- ifelse(metaMerged$EGFP > 0, as.character(metaMerged$EGFP_assign), ifelse(metaMerged$mCherry > 0, as.character(metaMerged$mCherry_assign), as.character(metaMerged$zero_assign)))


metaMerged <- metaMerged[metaMerged$SampleName != "Filter",]
table(metaMerged$initialGroup)
```

    ## 
    ##        G19 EGFP Adult Alone Late        G19 EGFP Adult CtrlRescue 
    ##                             1910                             1871 
    ##     G19 mCherry Adult Alone Late         G19 mCherry Adult Rescue 
    ##                              494                              548 
    ##  G19 mCherry Neonatal Alone Late  G19 mCherry Neonatal CtrlRescue 
    ##                             1661                             3246 
    ##     G20 EGFP Neonatal Alone Late         G20 EGFP Neonatal Rescue 
    ##                               74                              206 
    ## Neg G19 mCherry Adult Alone Late Neg G20 EGFP Neonatal Alone Late 
    ##                              240                              160

``` r
dim(metaMerged)
```

    ## [1] 10410    18

``` r
sampleKey <- read.csv("data_for_import/HDsampleKey.csv")

metaMerged <- merge(metaMerged, sampleKey, by.x = "SampleName", by.y = "SampleName")

table(metaMerged$Group)
```

    ## 
    ##       G19 EGFP Adult Alone Late       G19 EGFP Adult CtrlRescue 
    ##                            1910                            1871 
    ##    G19 mCherry Adult Alone Late        G19 mCherry Adult Rescue 
    ##                             734                             548 
    ## G19 mCherry Neonatal Alone Late G19 mCherry Neonatal CtrlRescue 
    ##                            1661                            3246 
    ##    G20 EGFP Neonatal Alone Late        G20 EGFP Neonatal Rescue 
    ##                             234                             206

\#\#Filter Seurat object and export to H5AD for SCVI integration

``` r
merged <- subset(merged, cells = metaMerged$cellName)

dim(merged)
```

    ## [1] 38743 10410

``` r
saveRDS(merged, "RDS/HD_Competition_Preprocessed.rds")


library(SeuratDisk)
```

    ## Registered S3 method overwritten by 'SeuratDisk':
    ##   method            from  
    ##   as.sparse.H5Group Seurat

``` r
SaveH5Seurat(merged, filename = "H5AD/HD_Competition_Preprocessed2.h5Seurat", overwrite = T)
```

    ## Warning: Overwriting previous file H5AD/HD_Competition_Preprocessed2.h5Seurat

    ## Creating h5Seurat file for version 3.1.5.9900

    ## Adding counts for RNA

    ## Adding data for RNA

    ## No variable features found for RNA

    ## No feature-level metadata found for RNA

``` r
Convert("H5AD/HD_Competition_Preprocessed2.h5Seurat", dest = "h5ad", overwrite = T)
```

    ## Validating h5Seurat file

    ## Adding data from RNA as X

    ## Adding counts from RNA as raw

    ## Transfering meta.data to obs

\#\#Session Info

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
    ## other attached packages:
    ## [1] SeuratDisk_0.0.0.9020 textTinyR_1.1.7       Matrix_1.4-0         
    ## [4] sp_1.4-6              SeuratObject_4.1.0    Seurat_4.1.1         
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16            colorspace_2.0-2      deldir_1.0-6         
    ##   [4] ellipsis_0.3.2        ggridges_0.5.3        rgdal_1.5-28         
    ##   [7] rstudioapi_0.13       spatstat.data_3.0-0   leiden_0.3.9         
    ##  [10] listenv_0.8.0         bit64_4.0.5           ggrepel_0.9.1        
    ##  [13] fansi_1.0.2           codetools_0.2-18      splines_4.1.1        
    ##  [16] knitr_1.37            polyclip_1.10-0       jsonlite_1.7.3       
    ##  [19] ica_1.0-2             cluster_2.1.2         png_0.1-7            
    ##  [22] rgeos_0.5-9           uwot_0.1.11           shiny_1.7.1          
    ##  [25] sctransform_0.3.3     spatstat.sparse_3.0-0 compiler_4.1.1       
    ##  [28] httr_1.4.2            assertthat_0.2.1      fastmap_1.1.0        
    ##  [31] lazyeval_0.2.2        cli_3.3.0             later_1.3.0          
    ##  [34] htmltools_0.5.2       tools_4.1.1           igraph_1.3.0         
    ##  [37] gtable_0.3.0          glue_1.6.2            RANN_2.6.1           
    ##  [40] reshape2_1.4.4        dplyr_1.0.8           Rcpp_1.0.8           
    ##  [43] scattermore_0.8       vctrs_0.3.8           nlme_3.1-155         
    ##  [46] progressr_0.10.0      lmtest_0.9-39         spatstat.random_3.0-1
    ##  [49] xfun_0.29             stringr_1.4.0         globals_0.14.0       
    ##  [52] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1      
    ##  [55] irlba_2.3.5           goftest_1.2-3         future_1.24.0        
    ##  [58] MASS_7.3-55           zoo_1.8-11            scales_1.1.1         
    ##  [61] spatstat.core_2.4-2   promises_1.2.0.1      spatstat.utils_3.0-1 
    ##  [64] parallel_4.1.1        RColorBrewer_1.1-2    yaml_2.3.4           
    ##  [67] reticulate_1.24       pbapply_1.5-0         gridExtra_2.3        
    ##  [70] ggplot2_3.3.5         rpart_4.1.16          stringi_1.7.6        
    ##  [73] rlang_1.0.1           pkgconfig_2.0.3       matrixStats_0.61.0   
    ##  [76] evaluate_0.15         lattice_0.20-45       ROCR_1.0-11          
    ##  [79] purrr_0.3.4           tensor_1.5            patchwork_1.1.1      
    ##  [82] htmlwidgets_1.5.4     bit_4.0.4             cowplot_1.1.1        
    ##  [85] tidyselect_1.1.1      parallelly_1.30.0     RcppAnnoy_0.0.19     
    ##  [88] plyr_1.8.6            magrittr_2.0.2        R6_2.5.1             
    ##  [91] generics_0.1.2        DBI_1.1.2             withr_2.4.3          
    ##  [94] mgcv_1.8-38           pillar_1.7.0          fitdistrplus_1.1-6   
    ##  [97] survival_3.2-13       abind_1.4-5           tibble_3.1.6         
    ## [100] future.apply_1.8.1    hdf5r_1.3.5           crayon_1.5.1         
    ## [103] KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_3.0-3  
    ## [106] plotly_4.10.0         rmarkdown_2.11        grid_4.1.1           
    ## [109] data.table_1.14.2     digest_0.6.29         xtable_1.8-4         
    ## [112] tidyr_1.2.0           httpuv_1.6.5          munsell_0.5.0        
    ## [115] viridisLite_0.4.0
