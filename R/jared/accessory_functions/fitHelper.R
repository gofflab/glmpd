#Modified from Alina Spiegel 12/4/20
library(tictoc)
library(purrr)

#This function runs for each gene:pattern combo
fitHelper <- function(gene, cds, model_formula_string){
  tic(gene)
  #TODO: to be generalizable these either need to be grabbed by var name or all need to be there
  #TODO: Match ensembl ids to gene short name (how does monocle do it)
  
  #create data frame with all of pData and expression of our one gene. Gene is gene_short_name. SLOW
  cellDat <- cbind(as.data.frame(pData(cds)), data.frame(GeneExpression = exprs(cds)[which(fData(cds)$gene_short_name %in% gene),]))
  
  model_formula_str <- paste("GeneExpression", model_formula_string, sep = "")
  model_formula = stats::as.formula(model_formula_str)
  #FM_fit <- NA #pre-empt try() call by making failed results be NA
  safe_glm_nb <- safely(glm.nb)
  FM_fit <- safe_glm_nb(model_formula,
                   data = cellDat, link = log)
  #FM_summary = summary(FM_fit)
  #df = list(model = FM_fit, model_summary = FM_summary)
  time <- toc()
  return(FM_fit)
}


# fitHelper <- function(gene, cds, pattern){
#   tic(gene)
#   #TODO: to be generalizable these either need to be grabbed by var name or all need to be there
#   #TODO: Match ensembl ids to gene short name (how does monocle do it)
#   datFrame <- data.frame("PatternWeights"= pattern_df[,paste0("cellPattern",pattern)], 
#                          "GeneExpression" = exprs(cds)[which(fData(cds)$gene_short_name %in% gene),]+1, #This is not elegant
#                          "CellType" = pData(cds)[,"cell_type"])
#   
#   FM_fit <- NA #pre-empt try() call by making failed results be NA
#   try(FM_fit <- glm.nb(GeneExpression ~ 0+  PatternWeights * CellType,
#                        data = datFrame, link = log))
#   #FM_summary = summary(FM_fit)
#   #df = list(model = FM_fit, model_summary = FM_summary)
#   toc()
#   return(FM_fit)
# }

#fitHelper("Hbb-bs",lmmp, 37, model_formula_string = "GeneExpression ~ PatternWeights + CellType")

###Original Code follows
# fit_helper <- function(thisPattern, cds, model_formula_string, distr,start=NULL){
#   datFrame <- as.data.frame(merge(as.data.frame(pData(cds)),thisPattern,by=0,sort=FALSE))
#   datFrame$Genotype <- factor(datFrame$Genotype,levels=levels(datFrame$Genotype)[2:1])
#   model_formula_str <- paste("y", model_formula_string, sep = "")
#   model_formula = stats::as.formula(model_formula_str)
#   FM_fit <- speedglm(model_formula, data=datFrame, family=distr,start = start,trace=TRUE,maxit=100)
#   FM_summary = summary(FM_fit)
#   df = list(model = FM_fit, model_summary = FM_summary)
# } 
# # run a model for each of the patterns (gamma distribution)
# model_formula_str <- "~0 + assigned_cell_type + Genotype:assigned_cell_type + num_genes_expressed + batch + Sex"
# #model_formula_str <- "~num_genes_expressed"
# fits <- apply(pWeights,1,fit_helper,cds=dat_ALS,model_formula_string=model_formula_str,distr=Gamma("log"),start = rep(.001,27))
# fit_summary <- summary(fits)
# fm <- tibble::as_tibble(purrr::transpose(fits))
# fm <- fm %>%dplyr::mutate(status = purrr::map(.f = purrr::possibly(monocle3:::extract_model_status_helper, NA_real_), .x = model)) %>% tidyr::unnest(status)
# fit_coefs <- coefficient_table(fm)



