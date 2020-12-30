#########################################
# Model Diagnostics
#########################################

#' Residual Plot
#'
#' @param model A model fit using fit_models
#' @param gene A string of gene name (must be in <where is this matching?>)
#'
#' @return
#' @import monocle3
#' @import ggplot2
#' @import MASS
#' @export
#'
#' @examples
#' <add example code here for how to use this function>
resid<-function(model, gene){
  resid_raw <- residuals.glm(model)
  resid_std <- rstandard(model)
  resid_part <- residuals.glm(model, type= "partial")
  var(resid_std)

  plot(resid_raw, resid_std)
}
