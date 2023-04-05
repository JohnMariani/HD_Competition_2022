Identification of important transcription factors in competition
================
Nguyen Huynh
01/03/2023

Load libraries and data

``` r
object <- readRDS("RDS/mergedUpdated.rds")
# Subset to gpc 
gpc <- subset(object, subset = (CellType == "GPC"))
# new name 
gpc$newGroup <- plyr::mapvalues(gpc$otherGroup,
                                unique(gpc$otherGroup), 
                                c("Alone_Young_WT", "Isograft_Young_WT", "Allograft_Young_WT", "Alone_Aged_WT", "Isograft_Aged_WT", "Alone_Aged_HD", "Allograft_Aged_HD"))
gpc$newGroup <- factor(gpc$newGroup, levels = c("Alone_Aged_HD", "Alone_Aged_WT" , "Alone_Young_WT",  "Allograft_Aged_HD", "Allograft_Young_WT", "Isograft_Aged_WT", "Isograft_Young_WT"))
# Environment factor: Whether the group is going to be in competition or not competing
fg <- c("Not_Competing", "In_Competition", "In_Competition", rep(c("Not_Competing", "In_Competition"), times = 2))
gpc$Environment <- plyr::mapvalues(gpc$otherGroup, unique(gpc$otherGroup), fg)
my_color = c("darkslategray1", "deepskyblue", "blue4", "seagreen1", "seagreen", "lightpink", "deeppink4")
# saveRDS(gpc, "RDS/gpc.RDS")
```

``` r
# DEG generated in script 03
allFC <- read.table("DE/allLogFC.txt")
res <- list(allFC[allFC$Comparison == "WT.Adult.CtrlRescue.vs.WT.Neonatal.CtrlRescue", ], 
            allFC[allFC$Comparison == "WT.Adult.Rescue.vs.HD.Neonatal.Rescue", ])
names(res) <- c("Isograft_Young_WT_vs_Isograft_Aged_WT","Allograft_Young_WT_vs_Allograft_Aged_HD")
```

``` r
# WGCNA results from 13_WGCNA_runTOM.Rmd
load("output/WGCNA_output.RData")
gene_module <- as.data.frame(dynamicColor)
MEG <- readRDS("output/MEG.RDS")
rm(filt, geneTree, MEs); gc()
```

``` r
# SCENIC results 
auc <- read.csv("pySCENIC/outs/auc.txt", row.names=1)
colnames(auc) <- gsub("\\..*", "", colnames(auc))
##-----
regs <- read_json("pySCENIC/outs/regulonsTarget.json", simplifyVector = T) # TF to target pairs
weight <- read_json("pySCENIC/outs/regulonsWeight.json", simplifyVector = T) # TF to gene2weight based on target order
names(regs) <- gsub("\\(.*", "", names(regs))
names(weight) <- gsub("\\(.*", "", names(weight))
### Make a data frame of all transcription factors and its respective targets and weights
reg_res <- data.frame(TF=character(), Target=character(), Weight=numeric())
for (i in 1:length(regs)) {
  tf <- names(regs)[i]
  numTarget <- length(regs[[i]])
  tmp <- data.frame(TF=rep(tf, numTarget), Target=as.character(regs[[i]]), Weight=as.numeric(weight[[i]]))
  reg_res <- rbind(reg_res, tmp)
}
```

## PART 1: TFs highly correlated to winning/losing

1.  Linear model to pick TFs. Focus only on clear winning/losing signal.
    So the isograft and allograft with new engraftment as “winning” and
    old engraftment as “losing”.

2.  Logistic regression to determine which TF activity gives the largest
    weight to winning/losing. Biological significance: All TFs are
    acting in concerts, but which ones have the most weight in
    establishing winning/losing identities? Winning/losing is assigned
    based on histology data.

``` r
gpc <- SetIdent(gpc, value = "newGroup")
gpc_subset_wl <- subset(gpc, subset = (newGroup %in% c("Allograft_Aged_HD", "Allograft_Young_WT", "Isograft_Aged_WT", "Isograft_Young_WT"))) 
subset_wl <- auc[row.names(auc) %in% colnames(gpc_subset_wl),] %>% as.matrix()
```

``` r
all.equal(row.names(subset_wl), row.names(gpc_subset_wl@meta.data))
wl <- ifelse(gpc_subset_wl$newGroup %in% c("Allograft_Young_WT", "Isograft_Young_WT"), 1, 0)

cv.lasso <- cv.glmnet(x=subset_wl, y = wl, alpha = 1, family = "binomial")
model <- glmnet(x = subset_wl, y = wl, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)
coef <- as.data.frame(as.matrix(coef(model)))
```

``` r
# Cohort of TFs that are downregulated in wining cells 
losing_tfs <- row.names(coef)[coef$s0 < 0 & row.names(coef) != "(Intercept)"]
# Cohort of TFs that are upregulated in winning cells 
winning_tfs <- row.names(coef)[coef$s0 > 0 & row.names(coef) != "(Intercept)"]
```

Mean expression between groups to keep TFs with low activity in HD and
low activity in old

``` r
gpc_auc <- auc[row.names(auc) %in% row.names(gpc@meta.data),]
all.equal(row.names(gpc@meta.data), row.names(gpc_auc))
winning_tfs_ttest <- c()
for (tf in winning_tfs) {
  sub_auc <- gpc_auc[, tf, drop = FALSE] 
  sub_auc$newGroup <- gpc$newGroup
  # test 1: new WT vs old HD
  mean_test_1_old <- sub_auc[sub_auc$newGroup == "Allograft_Aged_HD",][,1] %>% mean()
  mean_test_1_new <- sub_auc[sub_auc$newGroup == "Allograft_Young_WT",][,1] %>% mean()
  test1 <- mean_test_1_new - mean_test_1_old
  # test 2: new WT vs old WT
  mean_test_2_old <- sub_auc[sub_auc$newGroup == "Isograft_Aged_WT",][,1] %>% mean()
  mean_test_2_new <- sub_auc[sub_auc$newGroup == "Isograft_Young_WT",][,1] %>% mean()
  test2 <- mean_test_2_new - mean_test_2_old
  if (test1 > 0 & test2 > 0) {
    winning_tfs_ttest <- c(winning_tfs_ttest, tf)
  }
}
```

An extra filtering step to get winning TFs with high weights for the
above lasso but also whose targets are enriched for DEGs

``` r
RES <- res[c("Isograft_Young_WT_vs_Isograft_Aged_WT", "Allograft_Young_WT_vs_Allograft_Aged_HD")]
ranks <- lapply(RES, function(x) {z <- x$logFC; names(z) <- x$Gene; return(z)})
ranks <- lapply(ranks, function(x) x[order(x, decreasing = TRUE)])
ranks <- lapply(ranks, function(x) x[!is.na(x)]) # Remove those with NA log2FC
ranks <- lapply(ranks, function(x) x[!x == 0]) # Remove those with 0 log2FC
# 
fgseaRes <- lapply(ranks, function(x) fgsea(pathways = regs[winning_tfs_ttest], stats = x, minSize = 4, maxSize = 5000))
# Filter for signif 
fgseaRes_filt <- lapply(fgseaRes, function(x) x[which(x$padj < 1e-3 & x$NES > 0), ]) # Because these are winning signals, they are supposed to be enriched for upregulated genes since their AUC is positive. Note how SCENIC only do positive correlation right now. 
```

``` r
# Make a table object for plotting with ggplot later 
for (i in 1:length(fgseaRes_filt)) {
  fgseaRes_filt[[i]]$comparison <- rep(names(fgseaRes_filt)[i], nrow(fgseaRes_filt[[i]]))
}
fgsea.df <- Reduce(rbind, fgseaRes_filt)
# Add in total regulon size 
total_size <- c()
for (pw in fgsea.df$pathway) {
  total_size <- c(total_size, length(regs[[pw]]))
}
fgsea.df$total_size <- total_size
fgsea.df$propor <- fgsea.df$size/fgsea.df$total_size
# saveRDS(fgsea.df, "RDS/fgsea_filt_df.RDS")
```

List of final winning TFs

``` r
winning_tfs_filt <- unique(fgsea.df$pathway)
winning_tfs_filt
```

## PART 2: Interesting modules

``` r
# Background signifies what cells are already waiting in the environment when the young newly-engrafted cells are coming in. 
features_map <- data.frame(sample = c("Alone_Young_WT", "Isograft_Young_WT", "Allograft_Young_WT", "Alone_Aged_WT", "Isograft_Aged_WT", "Alone_Aged_HD", "Allograft_Aged_HD"), 
                           genotype = c("WT", "WT", "WT", "WT", "WT", "HD", "HD"), 
                           age = c("Young", "Young", "Young", "Aged", "Aged", "Aged", "Aged"), 
                           background = c("mouse", "WT", "HD", NA, NA, NA, NA))
gpc$genotype <- plyr::mapvalues(as.character(gpc$newGroup), features_map$sample, features_map$genotype)
gpc$age <- plyr::mapvalues(as.character(gpc$newGroup), features_map$sample, features_map$age)
gpc$background <- plyr::mapvalues(as.character(gpc$newGroup), features_map$sample, features_map$background)
```

Relative importance analysis with genotype and age as contributors
(background is not considered here because they don’t have all the
information)

``` r
# Matrix to test for multicolinearity 
all.equal(row.names(gpc@meta.data), row.names(MEG))
tmp_mat <- data.frame(ME = MEG[,1], genotype = gpc$genotype, age = gpc$age)
tmp_mod <- lm(ME ~ genotype + age, data = tmp_mat)
# VIF < 5 so all good to go
car::vif(tmp_mod)
```

``` r
# Build a model where ME is explained by the additive effect of both genotype and age
rsq  <- data.frame(row.names = gsub("ME", "", colnames(MEG)), R2=rep(0, ncol(MEG)), pval = rep(NA, ncol(MEG)))
proR <- data.frame(matrix(ncol=5, nrow=0))
colnames(proR) <- c("confidence", "lower", "upper", "covariate", "ME")
if(!exists("plots/Relative Covariate Importance")) {dir.create("plots/Relative Covariate Importance")}
if(!exists("tables/Relative Covariate Importance")) {dir.create("tables/Relative Covariate Importance")}
for (i in 1:ncol(MEG)) {
  me <- MEG[, i]
  me_name <- gsub("ME", "", colnames(MEG)[i])
  data_for_lm <- data.frame(ME = me, genotype = as.factor(gpc$genotype), age = as.factor(gpc$age))
  fit <- lm(ME ~ genotype + age, data = data_for_lm)
  # R-squared and p-val for the fit
  rsq[me_name, "R2"] <- summary(fit)$r.squared
  rsq[me_name, "pval"] <- lmp(fit)
  # Relaimpo
  cov_keep <- c("genotype", "age")
  tmp <- calculate_relative_importance(fit = fit, MEcolor = me_name, suffix = me_name)
  colnames(tmp) <- c("confidence", "lower", "upper", "covariate", "ME")
  proR <- rbind(proR, tmp)
}
# Adjust p-values for multiple comparisons 
rsq$padj <- p.adjust(as.vector(rsq$pval), method = "BH")
```

Visualize

``` r
# Rearrange proR by contribution to age, then set factor for ggplot
x <- proR %>% filter(covariate == "age")
x <- x[order(x$confidence, decreasing = TRUE),]
proR$ME <- factor(proR$ME, levels = x$ME)
col_pal <- RColorBrewer::brewer.pal(8, "Set3")[4:5]
# 
p3 <- ggplot(proR, aes(x=ME, y=confidence, fill = covariate)) + geom_col(position = position_dodge2(preserve = "single")) + scale_fill_manual(values = col_pal) + theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) + ylab(bquote('Proportion'~R^2)) 
# saveRDS(proR, "RDS/proR.RDS")
```

## PART 3: Modules harboring targets that are enriched in winning TF regulons

``` r
regs_filt <- list()
for (tf in winning_tfs_filt) {
  regs_filt[[tf]] <- reg_res$Target[reg_res$TF == tf]
}

module_to_regulons <- data.frame(ME = character(), TF = character(), pval = numeric(), intersection_size = numeric(), regulon_size = numeric(), module_size = numeric())
for (m in colnames(MEG)) {
  for (tf in names(regs_filt)) {
    modular_gene_member <- row.names(gene_module)[gene_module$dynamicColor == gsub("ME", "", m)]
    overlap_test <- newGeneOverlap(modular_gene_member, 
                                   regs_filt[[tf]], 
                                   genome.size = nrow(gene_module))
    overlap_test <- testGeneOverlap(overlap_test)
    module_to_regulons <- rbind(module_to_regulons, 
                                data.frame(ME = gsub("ME", "", m), 
                                           TF = tf, 
                                           pval = overlap_test@pval, 
                                           intersection_size = length(overlap_test@intersection),
                                           regulon_size = length(regs_filt[[tf]]),
                                           module_size = length(modular_gene_member)))
  }
}
module_to_regulons$padj <- p.adjust(as.vector(module_to_regulons$pval), method = "BH")
module_to_regulons$intRegulon_prop <- module_to_regulons$intersection_size/module_to_regulons$regulon_size
module_to_regulons$intModular_prop <- module_to_regulons$intersection_size/module_to_regulons$module_size
# filter
module_to_regulons_filt <- module_to_regulons[module_to_regulons$padj < 0.01 & module_to_regulons$intRegulon_prop > 0.2, ]
# saveRDS(module_to_regulons_filt, "RDS/module_to_regulons_filt.RDS")
```

## PART 4: Making of network tables

### Edges

``` r
tf.df <- reg_res[reg_res$TF %in% winning_tfs_filt, ]
# write.table(tf.df, "output/regulons_full.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# Prune for targets in prioritized modules 
picked_genes <- row.names(gene_module)[which(gene_module$dynamicColor %in% "blue")]
tf.df <- tf.df[tf.df$Target %in% picked_genes,]
# write.table(tf.df, "output/regulons_filt_module.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Prune for shared targets 
shared_list <- list()
for (tf in winning_tfs_filt) {
  shared_list[[tf]] <- tf.df$Target[tf.df$TF == tf]
}

counter <- data.frame(target = character(), counter = numeric()) 

for (i in tf.df$Target[tf.df$TF == "MYC"]) {
  counter_tmp <- length(unlist(lapply(shared_list, function(x) x[x==i])))
  counter <- rbind(counter, 
                   data.frame(target = i, counter = counter_tmp))
}
# MYC targets that are shared with at least one other TFs 
shared <- counter$target[counter$counter >= 2]
tf.df <- tf.df[tf.df$Target %in% shared, ]
# write.table(tf.df, "output/regulons_filt_shared_targets.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE) # 182 3
```

### Nodes

``` r
node <- data.frame(node = unique(c(tf.df$TF, tf.df$Target)))
# Module each node belongs to 
node$ME <- plyr::mapvalues(node$node, row.names(gene_module),gene_module$dynamicColor, warn_missing = FALSE)
# How much it is explained by age - 1-genotype
proR_age <- proR[proR$covariate == "age", ]
node$age_strength <- plyr::mapvalues(node$ME,proR_age$ME, proR_age$confidence, warn_missing = FALSE)
# Type of TF or targets 
node$Type <- ifelse(node$node %in% winning_tfs_filt, "TF", "Target")
# write.table(node, "output/regulons_filt_shared_targets_node.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
```

## END

## Softwares

``` r
sessionInfo()
```
