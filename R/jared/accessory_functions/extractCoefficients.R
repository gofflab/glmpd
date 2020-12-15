library(dplyr)

extractCoefficients <- function(model, params = names(coef(model))){

  #percent confidence interval for betas
    #This takes a super long time, what interpretations is this necessary for?

  #TODO: warning/error if parameters not in model
  coef_conf_confint <- as.data.frame(confint(model, parm = params)) %>% tibble::rownames_to_column(var = "coefficient")

  coefs <- as.data.frame(coef(summary(model))) %>% tibble::rownames_to_column(var = "coefficient")
  coef_df <- left_join(coefs, coef_conf_confint, by = "coefficient")

  return(coef_df)
  }


#function call
if(0){
  extractCoefficients(extracted_models[[1]], params = c("cellPattern37"))
}


