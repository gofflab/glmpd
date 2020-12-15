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
#' @examples
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
