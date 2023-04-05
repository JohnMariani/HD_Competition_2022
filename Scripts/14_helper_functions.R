# Function to extract p value from linear model object
# Source: http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html
# STAR
lmp <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
} # end lmp()
####################################################################################
####################################################################################
# correlate_to_variables()
# Function to calculate correlation coefficients between each module eigengene/ or principal components and each
# known covariate. The results are written into a table. Function also returns a data frame 
# of module eigengene/principal component v. covariate pvalues and r squared. 
# Parameters:
#  RESPONSE (d.f):      module eigengene or principal components 
#  sample_data (d.f):   sample information with appropriate covariates to test
#  cov_type (vector):   vector indicating covariate types (1=categorical or 2=continuous)
#  SIGNIFLEVEL (num):   how significant do we want to cut BH padj at 
#  suffix (string):     suffix that is appended to file names to keep them unique
#
correlate_to_variables <- function(RESPONSE, sample_data, cov_type, SIGNIFLEVEL, suffix) {
  # Big data frame containing all MEs/PCs (col 1) and covariate (col 2), and their corresonding lm results 
  lmResults.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
  colnames(lmResults.df) <- c("response","covariate","type","adj_r_squared","pvalue", "BH", "slope")
 
  for (p in 1:ncol(sample_data)) {
    tryCatch({
    if (length(unique(sample_data[,p])) > 1 ) { # if statement to disregard covariate where there is only one factor
    ### predictor specification
    if (cov_type[p] == 1) {
      type <- "categorical"
      # Only goes through this if there are more than 1 level 
      predictor <- sample_data[,p]
      lvl <- levels(as.factor(predictor))
      if (length(lvl) > 1) {
        predictor <- predictor
      } else {
        stop(paste0("There is only one level for covariate: ", colnames(sample_data)[p]))
      }
      
    } else if (cov_type[p] == 2) {
      type <- "continuous"
      predictor <- sample_data[,p]
    } else {
      cat("check typos - covariate is neither categorical or continuous")
    }# end else if
    
    predictor_name <- colnames(sample_data)[p]
     # Small data frame containing per ME/PC and all covariate. 
    lmResults <- lmResults.df[0,]
    for (r in 1:ncol(RESPONSE)) {
      
      ### module eigengene or principal component is the response
      response <- RESPONSE[,r]
      response_name <- colnames(RESPONSE)[r]
      
      ### linear model test
      test <- lm(response ~ predictor)
      
      ### extract corelation coefficient and P Value
      coef <- summary(test)$adj.r.squared
      pval <- lmp(test)
      slope <- test$coefficients[[2]]
      
      tmp <- data.frame(response = response_name, 
                        covariate = colnames(sample_data)[p], 
                        type = type, 
                        adj_r_squared = coef, 
                        pvalue = pval, 
                        BH = NA, 
                        slope = slope)

      lmResults <- rbind(lmResults, tmp)
    } # end for r
    ### Perform BH correction and write out linear module results 
    #lmResults$BH    <- p.adjust(as.vector(lmResults$pvalue), method = "BH")
    #lmResults <- lmResults[, order(lmResults$BH)]
    ### Pipe to bigger output 
    lmResults.df <- rbind(lmResults.df, lmResults)
    
    } else {
    cat(paste0("Covariate ", colnames(sample_data)[p], " only has 1 level of factors \n"))
  } # end if
  }, error=function(e){}) # end tryCatch
  } # end for p
  # Correct for multiple test per covariate 
  lmResults.df <- as.data.frame(lmResults.df %>% group_by(covariate) %>% mutate(BH = p.adjust(pvalue, method = "BH")))
  lmResults.df$signif <- ifelse(lmResults.df$BH <= SIGNIFLEVEL, "Yes", NA)
  if (!is.na(suffix)) {
    write.table(lmResults.df, paste0("./tables/lmResults_", suffix, ".txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
  return(lmResults.df)
}


####################################################################################
####################################################################################
# plotlmVar()
# Function to plot each ME against known covariate. 
# Heatmap is highlighted by adjusted BH value 

# Parameters:
#  lmRes (d.f):         Linear regression result - from correlate_to_variables
#  suffx (string):      What is the title of the plot? 
#  signif (num):        What is the cut-off for BH padj to highlight 

plotlmVar <- function(lmRes, suffix, signif) {
  # PCov.df = result from run_pca
  df <- lmRes 
  df$log <- -log10(df$BH)
  df$signif <- ifelse(df$BH >= signif, NA, round(df$BH, 3))
  p <- ggplot(df, aes(covariate, response, fill=log)) + geom_tile()
  p <- p + scale_fill_gradient2(low="white", mid="lightyellow", high="firebrick", name="-log10(p-value)")
  p <- p + geom_text(aes(covariate, response, label=signif), color="white", size=3)
  p <- p + ggtitle(suffix) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  return(p)
}

####################################################################################
####################################################################################
# calculate_relative_importance()
# Function to calculate relative importance of known covariates in the module
# eigengene (ME) variation. The importance is measured by the relative contribution
# of each covariate towards the explanatory value of the additive R^2. The results
# are saved in output files. See documentation of relaimpo package.
# Function returns a data frame of each covariate by ME - and their confidence intervals 
#
# Parameters:
#  suffix (string):     suffix saved results 
#  fit (lm):            linear model object containing the additive model with ME as
#                       the response variable and known covariates as predictors
#  MEcolor (character): color label identifying the ME in fit
#  cov_keep (vector):   a vector of covariate considered for this linear regression. hard-coded.
# STAR
calculate_relative_importance <- function(fit, MEcolor, suffix) {
  
  bootfit <- NULL
  # may come back as NULL when dependent covariates produce a singlularity error when
  # calculating the bootstrap intervals
  while (is.null(bootfit)) {
    bootfit <- boot.relimp(fit, b = 200, type = c("lmg", "first", "last"), rank = TRUE, diff = TRUE, rela = TRUE, fixed = TRUE)
  } # end while
  
  # plot the bootfit results
  pdf(paste0("./plots/Relative Covariate Importance/", suffix, "_Booteval Plot.pdf"),
      height = 6, width = 8, useDingbats = FALSE)
  plot(booteval.relimp(bootfit))
  graphics.off()
  
  bootout <- booteval.relimp(bootfit)
  
  res <- rbind(t(rbind(bootout@lmg, bootout@lmg.lower, bootout@lmg.upper)),
               t(rbind(bootout@first, bootout@first.lower, bootout@first.upper)), 
               t(rbind(bootout@last, bootout@last.lower, bootout@last.upper)))
  res <- as.data.frame(res)
  colnames(res)  <- c("confidence", "lower", "upper") 
  res$covariate  <- paste(rep(cov_keep, 3), rep(c("lmg", "first", "last"), each = length(cov_keep)), sep = "_")
  res$ME         <- rep(suffix, nrow(res))
  
  # save output to text files
  write.table(res,
              paste0("./tables/Relative Covariate Importance/", suffix, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # extract data frame to plot covariate contribution to ME
  res1           <- res[1:length(cov_keep),] # only intersted in the lmg method
  res1$covariate <- gsub("_lmg", "", res1$covariate)
  return(res1)
} # end calculate_relative_importance()

######################################################################################
######################################################################################
# plot_cov_relaimpo()
# Function to plot module to phenotype association by -log(BH padj) and relative contribution 
# of each covariate to investigated eigengenes
# Parameters: 
# lmRes (d.f):        data frame storing the information of each covariate and their BH adjusted p values to ME. Results from correlate_to_variables()
# rsq (d.f):          data frame of total R squared value for each ME achieved with the additive multi factor linear model taking all covariates together 
# proR (d.f):         data frame of each covariate and how they explain the ME (as a proportion of total R2). Results from calculate_relative_importance()
# signifME(v):        vector of significant MEs to highlight
# color_pal(v):       vector of colors to use - in order to keep covariate color consistent across different plots 
plot_cov_relaimpo <- function(lmRes, rsq, color_pal) {
  #
  lmRes$logBH <- -log10(lmRes$BH)
  p1 <- ggplot(lmRes, aes(x=response, y=logBH, fill=covariate)) + geom_col(position = position_dodge2(preserve = "single")) + scale_fill_manual(values = color_pal) + theme_bw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "top") + geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "firebrick") + ylab("-Log10(P)")
  #
  rsq$group <- rep(1, nrow(rsq))
  p2 <- ggplot(rsq, aes(x=ME, y=R2, col=group)) + geom_line(aes(group=group), color = "black") + theme_bw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + ylab(bquote(R^2))
  #
  con <- ifelse(unique(proR$ME) %in% signifME, 'firebrick', 'grey60')
  p3 <- ggplot(proR, aes(x=ME, y=confidence, fill = covariate)) + geom_col(position = position_dodge2(preserve = "single")) + scale_fill_manual(values = color_pal) + theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank(), legend.position = "none") + ylab(bquote('Proportion'~R^2))
  # Stitching 
  cowplot::plot_grid(p1, p2, p3, ncol=1, align = "v", labels = c("i", "ii", "iii"), label_size = 8, rel_heights = c(1,0.7,1))
} # end plot_cov_relaimpo

#######################################################################################
#######################################################################################
# compare()
# Function to loop through a list and run pairwise gene-set overlap test
# Function returns a data frame of each pairwise comparison and their corresponding p values 
# Parameters: 
# m (list):        a list of gene sets to compare 
# gs (numeric):    a number indicating the total number of genes compared (genome size)
compare <- function(m, gs) {
  combs <- t(combn(x= names(m), m = 2))
  
  comp_value <- apply(combs, 1, function(x) {
    go.obj <- newGeneOverlap(listA = m[[x[1]]],
                             listB = m[[x[2]]],
                             genome.size = gs)
    go.obj <- testGeneOverlap(go.obj)
    f_res <- go.obj@pval
    return(f_res)
  })
  df_out <- data.frame(combs, comp_value)
  names(df_out) <- c("Module1", "Module2", "Fisher_pVal")
  return(df_out)
}
########################################################################################
########################################################################################
# my_enrich.f()
# Function to loop through a list of enrichR results and return a bound data frame 
# Parameters: 
# enrichR_res (list):   list of data frame results of enrichR. When dbs is composed of multiple sources, this will be a list of data frame.
# pval_co (nr):         p value cut-off
# prefix (string):      name of the written table, if write-out is prefered
my_enrichr.f <- function(enrichR_res, pval_co, prefix){
  X <- plyr::ldply(enrichR_res, rbind)
  
  X <- X[order(X$Adjusted.P.value), ]
  
  signif_X <- X[which(X$Adjusted.P.value < pval_co),]
  
  X$Overlap <- gsub("\\/", "_", X$Overlap)
  if (!is.na(prefix)) {
  write.table(file=prefix, X, quote=FALSE, sep="\t")
  }
  return(signif_X)
}

########################################################################################
########################################################################################
# merge_GO.f()
# Function to loop through the previously resulted data frames of significant GO term and merge the similar terms by gene set overlaps 
# Function returns a data frame of significantly overlapped terms. Will need to be annotated manually 
# Parameters:
# signif_X (df):   data frame results of my_enriched.f()
# cs (numeric):    a numeric value indication the total size of genes within each cluster being investigated
merge_GO.f <- function(signif_X, cs) {
  ### Group together similar pathways that were returned more than once - based on the number of gene overlaps between those pathways
  gene_set <- strsplit(signif_X$Genes, ";")
  names(gene_set) <- signif_X$Term
  gene_set_merge <- compare(m = gene_set, gs = cs) # Requires the above compare(m) function
  gene_set_merge <- gene_set_merge[which(gene_set_merge$Fisher_pVal < 0.001),]

  return(gene_set_merge)
}

########################################################################################
########################################################################################
# annot.f()
# Function to rearrange the data frame of significantly enriched terms into genes as row names and with GO term it belongs to. 
# This will give the annotation data for heatmaps 
# Parameters:
# pway (df):       data frame results of my_enriched.f()
# allGO (str):     if a gene appears in more than 1 GO terms, do we report all terms? Yes=Report all terms or No=Only report the most signif term. 
annot.f <- function(pway, allGO) {
  my_GO <- pway
  v <- my_GO$Genes
  v <- strsplit(v, ";")
  y <- unique(unlist(v))
  
  ## What GO term does a gene belong to?
  group <- data.frame(Gene =y, Term=rep(0, times=length(y)))
  
  for (x in y) {
    term <- my_GO$Term[grep(x, my_GO$Genes)]
    beautiful <- strsplit(term," \\(")
    beautiful <- unlist(lapply(beautiful, function(z) z[1]))# This deletes the (GO: cat #)
    # If returns more than 1 terms, then get the first term (since my_GO is already sorted by p-value. The first time is likely the most significant term that contains the gene)
    if (allGO == "No") {
      beautiful <- rep(beautiful[1], length(beautiful))
    }
    tmp <- data.frame(Gene = rep(x, length(beautiful)), Term = beautiful)
    
    group <- rbind(group, tmp)
  }
  ## Genes that are not part of significant GO term are assigned "undertermined"
  group$Term <- gsub("0", "undetermined", group$Term)
  results <- group
  
  return(results)
}

#####################################################################################
#####################################################################################
# pullCigar.f()
# Function to return a data frame of genes of interest
# and see if they are up- or down- regulated in our list of DE comparions
# Parameters:
# gene_vector (chr):    a vector of gene names of interest
# dres_list (list):     a list of deseq2 results. Note that column names are hard coded 
# padj (nr):            padj cut-off to use to define whether a gene is signif or not
# This function returns a data frame of GOI as row names and comparison as colnames
#----------------------- baseMean log2FoldChange     lfcSE     stat       pvalue         padj geneName
#----- ENSG00000164600   153.9489       5.943271 0.8905399 6.673784 2.492910e-11 4.812584e-09  NEUROD6
require(dplyr)
pullCigar.f <- function(gene_vector, dres_list, padj) {
  df <- data.frame(row.names = gene_vector)
  sn <- lapply(dres_list, function(x) x[which(x$padj < padj),])
  sn_gene <- lapply(sn, function(x) x[x$geneName%in% gene_vector, c("geneName", "log2FoldChange")])
  df <- sn_gene %>% purrr::reduce(full_join, by = "geneName")
  colnames(df)[-1] <- names(sn_gene)
  return(df)
}

#####################################################################################
#####################################################################################
# pullUnique()
# Function to return a vector of genes that are unique to each list item 
# (in comparison to all ther items within the same list) 
# Parameters: 
# gene_list (list):      a list of gene IDs or gene names 
pullUnique <- function(gene_list) {
  return_vector <- vector()
  for (i in 1:length(gene_list)) {
    list_item  <- gene_list[[i]]
    list_other <- gene_list[!1:length(gene_list) %in% i]
    item  <- unlist(list_item)
    other <- unique(unlist(list_other))
    return_vector <- c(return_vector, item[! item %in% other])
  }
  return(return_vector)
}

#####################################################################################
#####################################################################################
# draw_heatmap()
# Function to return a heatmap based on the matrix provided 
# Function relates to the egos gene ontology term for top_node
# Parameters: 
# df (df):           a df with gene name, comparisons, and log2FC as values. Colnames is hard-coded
# annot_row (df):    a data frame to annotate the genes above 
# suffix (chr):      what to name the heatmap
draw_heatmap <- function(df, annot_row, suffix) {
  my_color <- c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(9, "Set3"))
  # values of colors so that it stays consistent among plots 
  my_color.sub <- c("ghostwhite", my_color[1:ncol(annot_row)-1])
  names(my_color.sub) <- c("",colnames(annot_row)[2:ncol(annot_row)][order(colnames(annot_row)[2:ncol(annot_row)])])
  # main title
  suffix <- paste0(suffix, " Targets")
  # Order genes
  data <- tidyr::spread(df[, 1:3], comparison, log2FC)
  row.names(data) <- data$geneName
  data <- data[,-1]
  if (nrow(data) > 1) {
    ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
  } else {
    ord <- 1
  }
  df$geneName <- factor(df$geneName, levels = unique(df$geneName)[ord])
  # x axis title 
  x_lab <- levels(df$comparison)
  x_lab <- paste0(gsub("_.*", "", x_lab), "\n", gsub(".*_", "", x_lab), "/CTR")
  x_lab <- gsub("Astro", "Astrocyte", x_lab)
  df$comparison <- plyr::mapvalues(df$comparison, levels(df$comparison), x_lab)
  # Main plot
  p <- ggplot(df, aes(y=geneName, x=comparison, fill=log2FC)) + geom_tile(color = "black")
  p <- p + labs(x="", y="") + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0), position = "right")
  p <- p + scale_fill_gradientn(values = scales::rescale(c(min(df$log2FC), -0.1, 0, 0.1, max(df$log2FC))), colours = c("blue", "blue4", "black", "yellow4", "yellow"))
  p <- p + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), axis.text.y=element_text(color = "black"), plot.margin = unit(c(2,2,2,0), "pt"), axis.ticks.y.left = element_line(color = "black")) 
  p <- p + geom_text(aes(y=geneName, x=comparison, label=padj), color="grey60")
  # annotation
  annot <- annot_row[which(annot_row$geneName %in% levels(df$geneName)),]
  annot <- tidyr::gather(annot, GO, annotation, 2:ncol(annot))
  annot$Target <- factor(annot$geneName, levels = levels(df$geneName))
  x <- unique(annot$annotation)
  x <- x[order(x)]
  annot$annotation <- factor(annot$annotation, levels = x)
  my_color.sub <- my_color.sub[x]
  h <- ggplot(annot, aes(x = GO, y = Target, fill = annotation)) + geom_tile(color = "black") + scale_fill_manual(values = my_color.sub, na.value="white") + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + theme_minimal()
  h <- h + labs(x="", y="", title=suffix) + theme(axis.text.y = element_blank(), axis.text.x = element_text(color = "black", angle = 30, hjust = 1), axis.ticks.x = element_blank(), axis.ticks.y.left = element_blank(), axis.ticks.y.right = element_line(color = "black"), legend.position = "none", plot.margin = unit(c(0,5,0,10), "pt"))
  # 
  P <- cowplot::plot_grid(h, p, ncol = 2, rel_widths = c(0.6, 1), align = "h")
  return(P)
}

#####################################################################################
#####################################################################################
plotexplorer <- function(transformed_data, BY, color_vector, sampleTable, my_title, jj1, jj2, dimension) {
  # function to compute PCA with percentage and plot them
  # transformed_data: a data frame of transformed TPM values with colnames as samples and rownames as genes 
  # sampleTable: based on this to get what the group color should be based on. 
  # color_vector: color BY group
  # my_title: What should we name the plot? 
  # return PCA plot 
  pcaRes <- prcomp(t(transformed_data))
  pcaCoord <- pcaRes$x
  eigs <- pcaRes$sdev^2
  percentage <- c(eigs[jj1]/sum(eigs), eigs[jj2]/sum(eigs))
  percentage <- round(percentage*100, 1)
  percentage <- paste0("(",percentage, "%)")
  ## IF color_vector is missing 
  if (missing(color_vector)) {
    color_vector <- RColorBrewer::brewer.pal(8, "Accent")
  }
  pcaCoord <- as.data.frame(pcaCoord)
  sampleTemp <- sampleTable[row.names(pcaCoord),]
  pcaCoord$Group <- sampleTemp[, BY]
  pcaCoord$Batch <- sampleTemp[, "study_id"]
  ##### 3D or 2D 
  if (dimension=="2D") {
  pcaCoord2 <- pcaCoord[, c("Group", paste0("PC", jj1), paste0("PC", jj2))]
  colnames(pcaCoord2) <- c("Group", "firstPC", "secondPC")
  p <- ggplot(pcaCoord2, aes(firstPC, secondPC, color=Group)) + geom_point(alpha=0.99) + theme_bw() 
 # p <- p + geom_text(aes(firstPC, secondPC, label=row.names(pcaCoord)), nudge_x = 1, nudge_y = 5, size = 2)
  p <- p + ggtitle(my_title) + scale_color_manual(values=color_vector)
  p <- p + xlab(paste0("PC", jj1, " ", percentage[1])) + ylab(paste0("PC", jj2," ", percentage[2])) 
  return(p)
  } else {
    return(pcaCoord)
  }
}


#####################################################################################
#####################################################################################
disorder_genes.f <- function(merged_data_frame, fgseaRes) {
  # Function returns a vector of genes of all the terms belonging to the same cluster
  # merged_data_frame: curated data frame with `representative_term` and `other_terms` as hard-coded columns in the data frame. Data were curated after clusterProfiler::simplifed and emapplot
  # fgseaRes: results from gseGO
  X <- fgseaRes@result
  disorder_gene_list <- list()
  for (i in unique(merged_data_frame$representative_term)) {
    goi <- vector()
    for (term in merged_data_frame$other_terms[which(merged_data_frame$representative_term == i)]) {
      if (term %in% X$Description) {
      x <- X$core_enrichment[which(X$Description == term)]
      x <- strsplit(x, "\\/")[[1]]
      goi <- c(goi, x)
      }
    }
    goi <- unique(goi)
    disorder_gene_list[[i]] <- goi
  }
  return(disorder_gene_list)
}

#####################################################################################
#####################################################################################
addgrids3d <- function(x, y=NULL, z=NULL, grid = TRUE,
                       col.grid = "grey", lty.grid = par("lty"),
                       lab = par("lab"), lab.z = mean(lab[1:2]),
                       scale.y = 1, angle = 40,
                       xlim=NULL, ylim=NULL, zlim=NULL){
#' Add grids to a scatterplot3d
#' 
#' @description The goal of this function is to add grids on an existing
#'  plot created using the package scatterplot3d
#' @param x,y,z numeric vectors specifying the x, y, z coordinates of points.
#'  x can be a matrix or a data frame containing 3 columns corresponding to
#'  the x, y and z coordinates. In this case the arguments y and z are optional
#' @param grid specifies the facet(s) of the plot on which grids should be drawn.
#'  Possible values are the combination of "xy", "xz" or "yz".
#'  Example: grid = c("xy", "yz"). The default value is TRUE to add grids only on xy facet.
#' @param col.grid,lty.grid color and line type to be used for grids
#' @param lab a numerical vector of the form c(x, y, len).
#'  The values of x and y give the (approximate) number of tickmarks on the x and y axes.
#' @param lab.z the same as lab, but for z axis
#' @param scale.y of y axis related to x- and z axis
#' @param angle angle between x and y axis
#' @param "xlim, ylim, zlim" the x, y and z limits (min, max) of the plot.
#' 
#' @note
#' Users who want to extend an existing scatterplot3d graphic with the
#'  function addgrids3d, should consider to set the arguments scale.y, angle, ...,
#'  to the value used in scatterplot3d.
#' 
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' 
#' @example
#' library(scatterplot3d)
#' data(iris)
#' scatterplot3d(iris[, 1:3], pch = 16, grid=T, box=F)
#' addgrids3d(iris[, 1:3], grid = c("xy", "xz", "yz"))
  
  if(inherits(x, c("matrix", "data.frame"))){
    x <- as.data.frame(x)
    y <- unlist(x[,2])
    z <- unlist(x[,3])
    x <- unlist(x[,1])
  }
  
  p.lab <- par("lab")
  
  angle <- (angle%%360)/90
  yz.f <- scale.y * abs(if (angle < 1) angle else if (angle >3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if (angle < 2) 1 - angle else angle - 3)
  
  
  # x axis range
  x.range <- range(x[is.finite(x)], xlim)
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 *lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  x <- x/x.scal
  x.range <- range(x.prty)/x.scal
  x.max <- ceiling(x.range[2])
  x.min <- floor(x.range[1])
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2]/x.scal))
    x.min <- min(x.min, floor(xlim[1]/x.scal))
  }
  x.range <- range(x.min, x.max)
  
  # y axis range
  y.range <- range(y[is.finite(y)], ylim)
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 *lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  y <- (y - y.add)/y.scal
  y.max <- (max(y.prty) - y.add)/y.scal
  if (!is.null(ylim))
    y.max <- max(y.max, ceiling((ylim[2] - y.add)/y.scal))
  
  # Z axis range
  z.range <- range(z[is.finite(z)], zlim)
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 *lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  z <- z/z.scal
  z.range <- range(z.prty)/z.scal
  z.max <- ceiling(z.range[2])
  z.min <- floor(z.range[1])
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2]/z.scal))
    z.min <- min(z.min, floor(zlim[1]/z.scal))
  }
  z.range <- range(z.min, z.max)
  
  # Add grid
  if ("xy" %in% grid || grid == TRUE) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + 
               z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + 
               (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }
  
  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min, 
             i + (yx.f * y.max), yz.f * y.max + z.max, 
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.max + temp,temp1 + i , col = col.grid, lty = lty.grid)
    
  }
  
  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min,  
             x.min + (i * yx.f) ,i * yz.f + z.max,  
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.min, i , col = col.grid, lty = lty.grid)
  }
  
}
#####################################################################################
#####################################################################################
find_lindep_alias <- function(covariate, sampleChild, response, force=FALSE) {
  # function to systematically find alias in the linear model 
  # covariate: a vector all all covariates to be considered
  # sampleChild: a data frame containing information on the samples. requires: covariate %in% colnames(sampleChild)
  # respose : a vector with the samples' response. requires: names(response) == row.names(sampleChild)
  # force: TRUE/FALSE, whether the function should be run even when *condition* is not being considered. Default: FALSE
  # return a vector of all linearly dependent covariates
  ##---- Systematically go through all covaratiates, and remove lindep covs 
  ##---- Once removed, build the linear model again, and keep removing if necessary 
  ##########-----##########
  # If there is no "condition" in covariate - break 
  if ((!"condition" %in% covariate) & (force == TRUE)) {
    cat("Please check the vector of covariate to make sure it contains 'condition' ")
    break
  }
  covariate_to_keep <- c()
  repeat{
    # We value condition above all other covariates, so reorder so that condition always comes first 
    if ("condition" %in% covariate) {covariate <- c("condition", covariate[!covariate %in% "condition"])}
    if (length(covariate) == 0) {
      break
    } else {
    # Covariate to consider all other covs to be dependent on 
    dep_on <- covariate[1]
    # Create a data frame for our linear fit 
    bigFit <- cbind(sampleChild[, colnames(sampleChild) %in% covariate, drop = FALSE], response)
    formula <- as.formula(paste("response ~ ", paste(covariate, collapse = "+" ))) 
    test_tmp <- lm(formula, data = bigFit, na.action = na.omit)
    # Get alias 
    alias <- alias(test_tmp)

    if (is.null(alias$Complete)) {
      # This is all the dep_on & the last covariates when all the variables are not linear-codependent anymore
      covariate_to_keep <- c(covariate_to_keep, covariate)
      break
    } else {
      alias <- alias$Complete
      alias <- alias[,grep(dep_on, colnames(alias)), drop = FALSE]
      alias <- apply(alias, 1, function(x) sum(x!=0))
      alias <- names(alias)[alias != 0]
      # These are coupled with the factor levels, we need to set it back to our covariate
      covariate_to_remove <- c()
      for (xx in covariate) {
        if (sum(grep(xx, alias)) != 0) {
          covariate_to_remove <- c(covariate_to_remove, xx)
        }
      } # end for xx 
      covariate_to_remove <- unique(covariate_to_remove)
      # Remove dep_on & reset covariate
      covariate <- covariate[!covariate %in% c(dep_on, covariate_to_remove)]
      covariate_to_keep <- c(covariate_to_keep, dep_on)
    } # end if
    } # end big if
  } # end repeat
  
  ##########-----##########
  # The second part is to remove highly correlated covariates
  # Because we already remove all the not full-ranked covs, we should not run into the problem with car::vif 
  formula <- as.formula(paste("response ~ ", paste(covariate_to_keep, collapse = "+" ))) 
  bigFit <- cbind(sampleChild[, colnames(sampleChild) %in% covariate_to_keep, drop = FALSE], response)
  test_tmp <- lm(formula, data = bigFit, na.action = na.omit)
  if (length(covariate_to_keep) > 1) {
    vif_res <- as.data.frame(car::vif(test_tmp))
    covariate_to_keep <- row.names(vif_res)[which(vif_res[,1] < 5 | row.names(vif_res) == "condition")]
  }
  
  # At the end, we should have a vector of all the good dep_on, and the last covariate vector 
  return(covariate_to_keep)
}
