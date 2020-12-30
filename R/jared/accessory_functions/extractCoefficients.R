library(dplyr)

extractCoefficients <- function(model){

  #percent confidence interval for betas
    #This takes a super long time, what interpretations is this necessary for?

  #TODO: warning/error if parameters not in model
  #coef_conf_confint <- as.data.frame(confint(model, parm = params)) %>% tibble::rownames_to_column(var = "coefficient")
  #print(coef(model))

  #pre-empt failed fits to return NA and continue
  coefs <- NA

  try(coefs <- as.data.frame(stats::coef(model)) %>% tibble::rownames_to_column(var = "coefficient"))
  #coef_df <- left_join(coefs, coef_conf_confint, by = "coefficient")

  return(coefs)
  }


#function call
if(0){
  extractCoefficients(extracted_models[[1]], params = c("cellPattern37"))
}


