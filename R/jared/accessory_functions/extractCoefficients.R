library(dplyr)

extractCoefficients <- function(model){
  #percent confidence interval for betas
  #This takes a super long time, what interpretations is this necessary for?
  coef_conf_confint <- as.data.frame(confint(mod.glm, parm = "cell_typeRBC")) %>% tibble::rownames_to_column(var = "coefficient")
  
  coef_df <- data.frame(beta_estimates =coef(model)) %>% tibble::rownames_to_column(var = "coefficient")
  
  coef_df <- left_join(coef_df, coef_conf_confint, by = "coefficient")
}

extractCoefficients(mod.glm)
