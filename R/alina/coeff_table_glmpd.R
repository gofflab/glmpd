# Takes the output of getProjectionDrivers.R and returns tables of coefficients
# for each term in the model (one table per pattern)

# adapted from monocle3::coefficient_table

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

# fits: list of model fits. output of getProjectionDrivers.R
# returns a list of coefficient tables
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