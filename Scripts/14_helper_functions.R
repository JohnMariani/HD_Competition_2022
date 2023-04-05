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