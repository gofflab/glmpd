###########################################################
# Functions for parallelization and fitting of GLM models
###########################################################

#' .fitGLM
#' This function parallelizes glm models:
#' Parallelizes across genes in a counts matrix
#' Parallelizes across patterns in a matrix of projected patterns
#' @param countsMatrix matrix of gene expression counts. genes x cells / samples.
#' @param annotDF A dataframe with annotations for each cell or sample in the gene expression dataset
#' @param geneDF (currently) optional dataframe with annotations for each gene in the dataset
#' @param model_formula_str A string specifying the model to fit, in the format "~ var1 + var2"
#' @param projected_patterns A matrix of projected pattern weights to use as an explanatory variable in the model fit for each gene. colnames are patterns, rownames are cells. Must have viable column names.
#' @param cores Integer defining the number of cores to use for parallelization of model fitting across genes.
#' @param exp_family The regression family to use (default: 'negbinomial')
#' @param clean_model Boolean.  Passed through to monocle3::fit_model
#' @param verbose Boolean. Verbose output for model fitting.
#' @param result Return full or summarized results. Summarized results are stripped of memory-intensive models and return coefficients from extractCoefficients().
#'
#' @return Returns a list of fitted models (output similar to monocle3::fit_model)
#' @import monocle3
#' @import MASS
#' @import Matrix
#' @import RhpcBLASctl
#' @import future
#' @import furrr
.fitGLM <- function(countsMatrix, annotDF, model_formula_str, projected_patterns, cores, geneDF= NULL){

  #from provided matrix
  pattern_names <- colnames(projected_patterns)

  #dependent on SCE structure
  genes <- rownames(countsMatrix)

  message(paste0("Fittings models for ", length(pattern_names), " patterns and ", length(genes), " genes"))
  if(sum(Matrix::rowSums(countsMatrix)== 0) > 0){
    warnings(paste0(sum(Matrix::rowSums(countsMatrix)== 0), " genes have zero expression and will not be successfully fit. It is recommended to remove them before running."))
  }

  #Not actually sure what these do... other than limit conflicts with furrr multicore operations
  RhpcBLASctl::omp_set_num_threads(1)
  RhpcBLASctl::blas_set_num_threads(1)

  #uses future parallelization structure
  #open multicore operations
  plan(multicore, workers = cores)
  full_glm_models <- furrr::future_map(pattern_names, function(pattern_name){
    message(paste0(date(),": working on pattern ",pattern_name))
    thisPatt <- projected_patterns[,pattern_name]

    #fit one pattern, all genes
    glm_models <- apply(countsMatrix,1,fit_helper,thisPattern=thisPatt,model_formula_string=model_formula_str,annotDF=annotDF)

    #tranpose lists to create a dataframe, each row is a gene
    #TODO: status extractor "success" or "fail"
    glm_models_df <- tibble::as_tibble(purrr::transpose(glm_models))

    #Add in gene fData to returned dataframe. Currently optional.
    if(!is.null(geneDF)){
      #TODO might need to assert length or gene names are the same. don't love the assumption they are the same.
      glm_models_df <- dplyr::bind_cols(geneDF, glm_models_df)
    }else{
      #TODO: JS: would names(glm_models) ever not be "genes" var from above? would be clearer
      glm_models_df <- glm_models_df %>%
        dplyr::mutate(gene_id = names(glm_models))
    }

    #if(result == "full_model"){
    #returns tibble with a column that includes full model for each gene. memory hog
    return(glm_models_df)
    #}

    #return coefficients, estimates, [anything else] to be easier on memory. Create a multitiered list.
    #top tiered list is pattern, second tiered list is gene, third tier is data.frame for those coefficients
    #TODO: this is essentially extractCoefficients(). replace with that function
    #else if(result == "summarized_model"){
    #  summarized_glm_model <- glm_model %>%
    #    monocle3::coefficient_table() %>%
    #    dplyr::select(-c(model, model_summary)) %>% as.data.frame()

    #  summarized_glm_model <- purrr:::map(genes, function(gene){
    #    summarized_glm_model %>% dplyr::filter(gene_id == gene)})

    #  names(summarized_glm_model) <- genes_of_interest

    #  return(summarized_glm_model)
    #}
  })

  #close multicore operations
  plan(sequential)

  names(full_glm_models) <- pattern_names

  return(full_glm_models)
}

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
                     data = datFrame, link = log, epsilon = 1e-3,
                     control = glm.control(maxit = 50))
    FM_summary = summary(FM_fit)
    df = list(model = FM_fit, model_summary = FM_summary)
  }, error = function(e){
    error_df <- list(model = NA, model_summary = NA)
    return(error_df) # return error messages; don't stop iterating (wondering if this is the best way to handle errors)
  })
}

#' fitGLMpd
#' Uses a GLM to identify the set of genes used in a pattern for a given cell type or condition
#' @param countsMatrix matrix of gene expression counts. genes x cells / samples.
#' @param annotDF A dataframe with annotations for each cell or sample in the gene expression dataset
#' @param model_formula_str A string specifying the model to fit, in the format "~ var1 + var2"
#' @param projected_patterns A matrix of projected pattern weights to use as an explanatory variable in the model fit for each gene. colnames are patterns, rownames are cells. Must have viable column names.
#' @param cores Integer defining the number of cores to use for parallelization of model fitting across genes.
#' @param exp_family The regression family to use (default: 'negbinomial')
#' @param clean_model Boolean.  Passed through to monocle3::fit_model
#' @param verbose Boolean. Verbose output for model fitting.
#' @param result Return full or summarized results. Summarized results are stripped of memory-intensive models and return coefficients from extractCoefficients().
#'
#' @return Returns a list of fitted models (output similar to monocle3::fit_model)
#' @import monocle3
#' @import SummarizedExperiment
#' @export
#'
# #' @examples
setGeneric("fitGLMpd", function(object, ..., verbose=TRUE) standardGeneric("fitGLMpd"),
           signature = "object")


setMethod("fitGLMpd",signature(object="cell_data_set"), function(object, model_formula_str, projected_patterns,cores=1, count_assay = "counts", geneDF = NULL){ #,exp_family="negbinomial",cores,clean_model=T,verbose=T, result = "full_model"){
  #TODO: also accept SingleCellExperiment (and SummarizedExperiment?)
  cds <- object

  countsMatrix <- assay(cds, count_assay)
  annotDF <- colData(cds) %>% as.data.frame()

  if(is.null(geneDF)){
    geneDF <- rowData(cds) %>% as.data.frame()
  }

  .fitGLM(countsMatrix, annotDF, model_formula_str, projected_patterns,cores, geneDF)
})

setMethod("fitGLMpd",signature(object="matrix"), function(object, annotDF, model_formula_str, projected_patterns,cores=1){ #,exp_family="negbinomial",cores,clean_model=T,verbose=T, result = "full_model"){

  countsMatrix <- object

  .fitGLM(countsMatrix, annotDF, model_formula_str, projected_patterns,cores)
})

###########################################################
# Functions for organizing model coefficients by genes
###########################################################


#' extractCoefficients
#' Extract coefficients from model fitting object. Each object of list requires a column exactly named "gene_id".
#' @param glm_models Named List (over patterns) of objects returned from "full_result" of fit_models_cds. Organizes each pattern by gene.
#' @param genes gene ids that exactly match the count matrix which was provided to fitGLMpd()
#' @return Returns a mutileveled list of fitted model coefficients. Hierarchy of organization: pattern (list), gene (list), coefficients (data.frame)
#' @import monocle3
#' @import dplyr
#' @import purrr
#' @export
#'
#TODO: assert genes length is the same as glm_models. Allow for "gene_id" to be set manually to allow other columns to match genes on.
#TODO: is there a way to get gene names from models instead of providing... limit misnaming error ... also assumes orders are the same as is
#Generates a multi-leveled list.
#First level: pattern (named list)
#Second level: gene (named list)
#Third level: beta coefficients, columns are: coefficient name, coefficient
extractCoefficients <- function(glm_models, genes){
  extracted_coefficients <- lapply(glm_models, function(glm_model){

    #Drop full models for memory effeciency. keeps parameters
    summarized_glm_model <- glm_model %>%
      monocle3::coefficient_table() %>%
      dplyr::select(-c(model, model_summary)) %>% as.data.frame()

    #organize with a list for each gene
    summarized_glm_model <- purrr:::map(genes, function(gene){
      summarized_glm_model %>% dplyr::filter(gene_id == gene)})

    names(summarized_glm_model) <- genes_of_interest

    return(summarized_glm_model)
  })
}

#' orderCoefficientsByGene
#' Order model coefficient estimates by genes, such as for visualization by heatmap
#' @param pattern_coefficient_list List (over patterns) of objects returned from extractCoefficients(), or "summarized_result" option of fit_models_cds().
#' @param model_terms_to_keep optional. limit returned dataframe to coefficients of interest.
#' @param filter_significance optional. Numeric value at which to filter q-value significance.
#' Genes with signficant coefficients for one or more  terms are returned. If not provided, all genes returned.
#' @param string optional. all model terms that contain this string will be used for signficance filtering.
#' @return Returns a list of fitted models (output similar to monocle3::fit_model)
#' @import dplyr
#' @import tibble
#' @export
#'
#Pass extracted coefficients to organize into a data frame with each gene as as a row ----
#TODO: Can make this more elegant with dplyr instead of lapply
#TODO: filter out failed genes
orderCoefficientsByGene <- function(pattern_coefficient_list, model_terms_to_keep = NULL, filter_significance, string){
  pattern_full_names <- names(pattern_coefficient_list)
  #loop thru each pattern
  params_list <- lapply(pattern_full_names, function(pattern_full_name){

    genes <- names(pattern_coefficient_list[[pattern_full_name]])
    #loop thru each gene
    param_df <- lapply(genes, function(gene){

      #pre-empt failed fits by creating a mini df with gene name
      #failed genes will have parameters filled in with NA during bind_rows() below
      param <- data.frame("gene_id" = gene)

      #extract all parameters of interest. "term" is parameter name, "estimate" is the estimated beta value
      try(param <- pattern_coefficient_list[[pattern_full_name]][[gene]] %>%
            dplyr::select(any_of(c("term","estimate","std_err","p_value","q_value","test_val"))) %>%
            tibble::column_to_rownames(var = "term") %>%
            t() %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "measure") %>%
            mutate("gene_id" = gene)
      )

      #if provided, restrict parameters to only these
      if(!is.null(model_terms_to_keep)){
        param <- param %>%
          dplyr::select(any_of(model_terms_to_keep))
      }

      #check if the minimum q-value for a gene is significant. If none are significant, return null to get rid of gene in df
      #
      if(!is.null(filter_significance)){
        #pre-empt failure (Why would this fail? gene doesn't exist?)
        min_q <- 1
        try(min_q <- param %>%
              filter(measure == "q_value") %>%
              dplyr::select(matches(string)) %>%
              min())


        if(min_q > filter_significance){
          #not significant
          return(NULL)
        }
      }

      return(param)

    })


    #generate for each pattern a df, with tibbles for each parameter. allows for filtering on qvalue etc
    tbl <- bind_rows(param_df) %>%
      group_by(measure) %>%
      nest() %>%
      ungroup()

    return(tbl)
  })

  names(params_list) <- pattern_full_names

  #list of dataframes, one df per pattern
  return(params_list)
}

#' organizeEstimates
#' Generate matrices of coefficients of interest for later plotting with heatmaps.
#' @param coefficients_list List (over patterns) of objects returned from orderCoefficientByGene()
#' @param terms_exact String. Names of terms to plot which exactly match the terms of the model.
#' @param terms_match String. specify beginning string to match terms to plot. Necessary for interaction terms.
#' @param feature which feature to pull from the model. Default ["estimate"] is beta. Options include "q_value"...
#' @param gene_name name of column which contains gene names. Default ["gene_id"]
#' @param tranpose boolean. Should matrices be tranposed. Default [F] returns genes as rows.
#' @return Returns a list of matrices.
#' @import dplyr
#' @import tibble
#' @import purrr
#' @export
#'
##coefficients_list should be output from orderCoefficientsByGene(). List for each pattern, containing a tibble with estimate, std error, pvalue...
##with tranpose = F (default), rownames are genes, colnames are model parameters
organizeEstimates <- function(coefficients_list, terms_exact, terms_match, feature = "estimate", gene_name = "gene_id", transpose = F){

  param_list <- purrr::map(names(coefficients_list), function(pattern_name){
    #contains tibbles for each measure, eg. estimate, std error, p-value.
    gene_coef_df<- coefficients_list[[pattern_name]] %>%
      dplyr::filter(measure == feature) %>%
      dplyr::select(data) %>%
      unnest(cols =data)

    #add gene_short_name
    #gene_coef_df <- left_join(gene_coef_df, gene_mapping)

    if(!is.null(rownames(gene_coef_df))){
      message("Warning: existing rownames may be removed")
      rownames(gene_coef_df) <- NULL
    }

    #remove genes that did not fit any model (ie rows that are all NA)
    rm.ind <- which(apply(gene_coef_df, 1, function(x) all(is.na(x))))
    if(length(rm.ind) > 0) {
      message(paste0("Removing ", length(rm.ind), " genes that did not have any successful fits"))
      gene_coef_df <- gene_coef_df[-c(rm.ind),]
    }

    #organize coefficients for these parameters, match parameter name exactly
    exact_parameters <- purrr::map(terms_exact, function(term){
      exact <- gene_coef_df %>%
        dplyr::select(c(all_of(term),all_of(gene_name))) %>%
        column_to_rownames(var = gene_name) %>%
        as.matrix()

      exact <- exact[,order(colnames(exact))]
      if(transpose){
        exact <- t(exact)
      }

    })
    names(exact_parameters) <- terms_exact

    #organize coefficients for these parameters, match parameter name with starting string
    matched_parameters <- purrr:::map(terms_match, function(term){
      match <- gene_coef_df %>% dplyr::select(c(starts_with(term), gene_name)) %>% column_to_rownames(var = gene_name) %>% as.matrix()

      match <- match[,order(colnames(match))]
      if(transpose){
        match <- t(match)
      }
    })
    names(matched_parameters) <- terms_match

    #for each pattern, return lists with matrix of coefficients for each provided parameter
    parameters_list <- c(exact_parameters, matched_parameters)

  })
  names(param_list) <- names(coefficients_list)
  return(param_list)
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
    models_df <- fits[[pattFitI]]

    # successfulModels <- models[unlist(lapply(models,length)) > 1]
    # fm <- tibble::as_tibble(purrr::transpose(successfulModels))
    # fm_named = fm %>% dplyr::mutate(id = names(successfulModels))
    #
    M_f = models_df %>% dplyr::mutate(terms = purrr::map2(.f = purrr::possibly(extractCoeffHelper_glm_nb,NA_real_), .x=model, .y= model_summary)) %>% tidyr::unnest(terms)
    fit_coefs = M_f %>% dplyr::group_by(model_component, term) %>%
      dplyr::mutate(q_value = stats::p.adjust(p_value)) %>%
      dplyr::ungroup()

    coeffTables[[pattFitI]] <- fit_coefs
  }
  names(coeffTables) <- names(fits)
  return(coeffTables)
}
