Processing Ingenuity Pathway Analysis output of HD and Ctr competition
striata
================
John Mariani
12/6/2022

``` r
library(data.table)
library(tidyverse)
library(xlsx)
```

# Load IPA comparison output data (Canonical, Functional and Upstream Signaling) pValues and Z-Scores

``` r
deRubric <- read.csv("data_for_import/deRubric.csv")

goCats <- c("EIF2 Signaling", "LARP1 Signaling", "MLXIPL Signaling", "mTOR Signaling", "MYC Signaling", "MYCN Signaling", "Oxidative Phosphorylation", "RICTOR Signaling", "Synthesis of protein", "TP53 Signaling", "YAP1 Signaling")


CANpval <- fread("data_for_import/CanonicalP.txt", skip = 2)
CANzscore <- fread("data_for_import/CanonicalZ.txt", skip = 2)

FUNpval <- fread("data_for_import/FunctionalP.txt", skip = 2)
FUNzscore <- fread("data_for_import/FunctionalZ.txt", skip = 2)

UPpval <- fread("data_for_import/UpstreamP.txt", skip = 2)
UPzscore <- fread("data_for_import/UpstreamZ.txt", skip = 2)
```

# For IPA data

``` r
UPpval$`Upstream Regulators` <- paste0(UPpval$`Upstream Regulators`, " Signaling")
UPzscore$`Upstream Regulators` <- paste0(UPzscore$`Upstream Regulators`, " Signaling")


names(FUNpval)
```

    ## [1] "Diseases and Bio Functions"                                             
    ## [2] "WT.Adult.CtrlRescue.vs.WT.Adult.Alone.sig2 - 2023-02-07 03:42 PM"       
    ## [3] "HD.Neonatal.Rescue.vs.HD.Neonatal.Alone.sig - 2022-12-10 05:35 PM"      
    ## [4] "WT.Adult.Alone.vs.HD.Neonatal.Alone.sig - 2022-12-10 05:36 PM"          
    ## [5] "WT.Adult.Alone.vs.WT.Neonatal.Alone.sig - 2022-12-10 05:36 PM"          
    ## [6] "WT.Neonatal.CtrlRescue.vs.WT.Neonatal.Alone.sig - 2022-12-10 05:36 PM"  
    ## [7] "WT.Adult.Rescue.vs.WT.Adult.Alone.sig - 2022-12-10 05:35 PM"            
    ## [8] "WT.Adult.CtrlRescue.vs.WT.Neonatal.CtrlRescue.sig - 2022-12-10 05:35 PM"
    ## [9] "WT.Adult.Rescue.vs.HD.Neonatal.Rescue.sig - 2022-12-10 05:35 PM"

``` r
formatIPAp <- function(ipa){
  temp <- ipa
  names(temp) <- gsub("(.*).sig.*","\\1",names(temp))
  names(temp) <- gsub("(.*).Sig.*","\\1",names(temp))
  temp[temp == "N/A"] <- 0
  names(temp)[1] <- "Signaling"
  temp <-temp  %>% tidyr::pivot_longer(cols = names(temp)[2:9], names_to = "Comparison", values_to = "Pval")
  return(temp)
}

formatIPAz <- function(ipa){
  temp <- ipa
  names(temp) <- gsub("(.*).sig.*","\\1",names(temp))
  names(temp) <- gsub("(.*).Sig.*","\\1",names(temp))
  temp[temp == "N/A"] <- 0
  names(temp)[1] <- "Signaling"
  temp <-temp  %>% tidyr::pivot_longer(cols = names(temp)[2:9], names_to = "Comparison", values_to = "Zscore")
  return(temp)
}

CANpval <- formatIPAp(CANpval)
CANzscore <- formatIPAz(CANzscore)
FUNpval <- formatIPAp(FUNpval)
FUNzscore <- formatIPAz(FUNzscore)
UPpval <- formatIPAp(UPpval)
UPzscore <- formatIPAz(UPzscore)




UPgo <- merge(UPpval, UPzscore)
FUNgo <- merge(FUNpval, FUNzscore)
CANgo <- merge(CANpval, CANzscore)

UPgo$category <- "Upstream"
FUNgo$category <- "Functional"
CANgo$category <- "Canonical"

GO <- rbind(UPgo, FUNgo)
GO <- rbind(GO, CANgo)
GO <- merge(GO, deRubric, by.x = "Comparison", by.y = "Comparison")

GO$Zscore <- as.numeric(GO$Zscore)
```

# Write out extended table

``` r
formatGO <- function(GO, paradigm){
  temp <- GO[GO$Paradigm == paradigm,]
  temp <- temp[,c(2,3,4,5,7)]
  names(temp)[[5]] <- "Comparison"
  temp <- temp[order(temp$Pval, decreasing = T),]
  return(temp)
}

temp <- formatGO(GO, "Allograft")


# write.xlsx(formatGO(GO, "Allograft"), "ExtendedDataTables/Extended Data Table 2-WT vs HD hGPC in competition DE GOs and pathways.xlsx", 
#            sheetName= "All IPA Terms",
#            col.names=T, row.names=F, append=T)
# 
# write.xlsx(formatGO(GO, "Isograft"), "ExtendedDataTables/Extended Data Table 4-Isogenic Young vs Older-DE GOs and pathways.xlsx", 
#            sheetName= "All IPA Terms",
#            col.names=T, row.names=F, append=T)
```

# Write out table for figure graphs

``` r
GO <- GO[GO$Signaling %in% goCats,]


GOlevels <- c("WT Ctr vs HD Ctr", "WT Allograft vs HD Allograft", "WT Allograft vs WT Ctr", "HD Allograft vs HD Ctr",
              "Young WT Ctr vs Aged WT Ctr", "Young WT Isograft vs Aged WT Isograft", "Young WT Isograft vs Young WT Ctr", "Aged WT Isograft vs Aged WT Ctr")




GO[GO$Pval < -log10(0.05),]$Pval <- -1



GO$Signaling <- factor(GO$Signaling, levels = goCats)
GO$Label <- factor(GO$Label, levels = GOlevels)

write.csv(GO, "data_for_import/GO.csv")
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
    ##  [1] xlsx_0.6.5        forcats_0.5.1     stringr_1.4.0     dplyr_1.0.7      
    ##  [5] purrr_0.3.4       readr_2.0.0       tidyr_1.1.3       tibble_3.1.3     
    ##  [9] ggplot2_3.4.0     tidyverse_1.3.1   data.table_1.14.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.1  xfun_0.24         rJava_1.0-6       haven_2.4.3      
    ##  [5] colorspace_2.0-2  vctrs_0.5.1       generics_0.1.0    htmltools_0.5.1.1
    ##  [9] yaml_2.2.1        utf8_1.2.2        rlang_1.0.6       pillar_1.6.2     
    ## [13] glue_1.4.2        withr_2.5.0       DBI_1.1.1         dbplyr_2.1.1     
    ## [17] modelr_0.1.8      readxl_1.3.1      lifecycle_1.0.3   munsell_0.5.0    
    ## [21] gtable_0.3.0      cellranger_1.1.0  rvest_1.0.1       evaluate_0.14    
    ## [25] knitr_1.33        tzdb_0.1.2        xlsxjars_0.6.1    fansi_0.5.0      
    ## [29] broom_0.7.9       Rcpp_1.0.7        backports_1.2.1   scales_1.2.1     
    ## [33] jsonlite_1.7.2    fs_1.5.0          hms_1.1.0         digest_0.6.27    
    ## [37] stringi_1.7.3     grid_4.1.0        cli_3.4.1         tools_4.1.0      
    ## [41] magrittr_2.0.1    crayon_1.4.1      pkgconfig_2.0.3   ellipsis_0.3.2   
    ## [45] xml2_1.3.2        reprex_2.0.1      lubridate_1.7.10  assertthat_0.2.1 
    ## [49] rmarkdown_2.9     httr_1.4.2        rstudioapi_0.13   R6_2.5.0         
    ## [53] compiler_4.1.0
