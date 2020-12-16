library(MASS)
library(dplyr)

# Adapted from monocle3, runs the model for each gene / pattern iteration
fit_helper <- function(thisGene, thisPattern, model_formula_string, annotDF){
  
  # sample x annotations, gene expression, and p weights 
  datFrame <- as.data.frame(merge(annotDF,thisGene,by=0,sort=FALSE))
  datFrame$patternWeights <- thisPattern
  
  # prepare model formula string
  model_formula_str <- paste("y", model_formula_string, sep = "")
  model_formula = stats::as.formula(model_formula_str)
  
  # run model
  tryCatch({
    FM_fit <- glm.nb(model_formula,
                     data = datFrame, link = log,
                     control = glm.control(maxit = 50))
    FM_summary = summary(FM_fit)
    df = list(model = FM_fit, model_summary = FM_summary)
  }, error = function(e){
    return(conditionMessage(e)) # return error messages; don't stop iterating
  })
}


# countsMatrix: genes x samples matrix with gene counts
# annotDF: samples x annotations matrix. includes condition(s) of interest
# projectedPs: patterns x sample matrix with projected p weights (learned with projectR)
# returns a list of model fits for each pattern. Each item in that list contains a fit
#      for each gene
getProjectionDrivers <- function(countsMatrix, annotDF, projectedPs, model_formula_str){
  
  # ensure that all genes in countsMatrix are expressed in at least one sample
  totalCounts <- apply(countsMatrix,1,sum)
  numZeroExpGenes <- sum(totalCounts==0)
  if(numZeroExpGenes>0)
    stop(paste("Number of genes with 0 counts across all samples:",numZeroExpGenes))
  
  # fit models for each pattern
  fits <- list()
  for (pattI in 1:(dim(projectedPs)[1])){
    print(paste0("running models for projection driver analysis: pattern ", pattI, " of ",dim(projectedPs)[1]))
    thisPatt <- projectedPs[pattI,]
    # fit a model for each gene
    fits[[pattI]] <- apply(countsMatrix,1,fit_helper,thisPattern=thisPatt,model_formula_string=model_formula_str,annotDF=annotDF)
  }
  names(fits) <- rownames(projectedPs)
  return(fits)
}