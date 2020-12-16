###########################################################
# Functions for parallelization and fitting of GLM models
###########################################################


#' Fit Models (cds)
#' This function wraps around monocle3::fit_model to parallelize glm models for a given pattern `pattern` across genes in a cds object
#' @param cds A monocle3 cell_data_set object
#' @param pattern A vector of projected pattern weights to use as an explanatory variable in the model fit for each gene.
#' @param cores Integer defining the number of cores to use for parallelization of model fitting across genes.
#' @param exp_family The regression family to use (default: 'negbinomial')
#' @param clean_model Boolean.  Passed through to monocle3::fit_model
#' @param verbose Boolean. Verbose output?
#'
#' @return Returns a list of fitted models (output similar to monocle3::fit_model)
#' @import monocle3
#' @import MASS
#' @export
#'
# #' @examples
fit_model_cds<-function(cds,pattern,cores,exp_family="negbinomial",clean_model=F,verbose=T){
  #TODO: Add 'pattern' to cds unless already there

  #TODO: As is this function assumes that 'cell_type' is a value in pData(cds).  We may want to just pass an argument for a model string if we are using values in pData(cds).

  # Create model formula string to include pattern of interest from 'pattern' argument
  model_str <- paste0("~0 + cell_type*",pattern)
  #exp_family <- "negbinomial" # set in args
  #ncores = min(length(genes), 16) # set in args
  glm_models <- monocle3::fit_models(cds = cds, model_formula_str = model_str, expression_family = exp_family, cores = cores,
                                       clean_model = clean_model, verbose = verbose)
  return(glm_models)
  }

# From Alina's code dump

#' fit_helper
#'
#' @param thisGene String of gene name
#' @param thisPattern String of pattern name
#' @param model_formula_string Model Formula String
#' @param annotDF Not sure what this argument is (alina?)
#'
#' @return
#' @importFrom stats as.formula
#' @import MASS
#' @import dplyr
#' @export
#'
# #' @examples
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
#' Title
#'
#' @param countsMatrix genes x samples matrix with gene counts
#' @param annotDF samples x annotations matrix. includes condition(s) of interest
#' @param projectedPs patterns x sample matrix with projected p weights (learned with projectR)
#' @param model_formula_str
#'
#' @return A list of model fits for each pattern. Each item in that list contains a fit for each gene
#' @export
#'
# #' @examples
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

#' extractCoeffHelper_glm_nb
#'
#' @param model A fitted GLM model
#' @param model_summary Not sure what this argument is offhand (Alina?)
#'
#' @return
#' @import tibble
#' @export
#'
# #' @examples
extractCoeffHelper_glm_nb <- function (model, model_summary)
{
  coef_mat <- model_summary[["coefficients"]]
  coef_mat <- apply(coef_mat, 2, function(x) {
    as.numeric(as.character(x))
  })
  row.names(coef_mat) = row.names(model_summary[["coefficients"]])
  colnames(coef_mat) = c("estimate", "std_err", "test_val",
                         "p_value")
  coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
  coef_mat$model_component = "count"
  return(coef_mat)
}


#' coeff_table_glmpd
#'
#' @param fits A List of model fits. output of getProjectionDrivers.R
#'
#' @return A list of coefficient tables
#' @import dplyr
#' @import tibble
#' @importFrom stats p.adjust
#' @import purrr
#' @export
#'
# #' @examples
coeff_table_glmpd <- function(fits){

  coeffTables <- list()
  # one coefficient table per pattern
  for (pattFitI in (1:length(fits)))
  {
    models <- fits[[pattFitI]]

    successfulModels <- models[unlist(lapply(models,length)) > 1]
    fm <- tibble::as_tibble(purrr::transpose(successfulModels))
    fm_named = fm %>% dplyr::mutate(id = names(successfulModels))
    M_f = fm_named %>% dplyr::mutate(terms = purrr::map2(.f = purrr::possibly(extractCoeffHelper_glm_nb,NA_real_), .x=model, .y= model_summary)) %>% tidyr::unnest(terms)
    fit_coefs = M_f %>% dplyr::group_by(model_component, term) %>%
      dplyr::mutate(q_value = stats::p.adjust(p_value)) %>%
      dplyr::ungroup()

    coeffTables[[pattFitI]] <- fit_coefs
  }
  names(coeffTables) <- names(fits)
  return(coeffTables)
}
