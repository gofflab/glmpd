library(monocle3)
library(tidyr)
library(tibble)
library(dplyr)
library(MASS)
source("./R/jared/accessory_functions/extractCoefficients.R")



# Load in required datsets and pweights --------------------------
#dataset to project into.
lmmp <- readRDS("./data/jared/TC_LMMP.rds")
#lmmp <- lmmp[fData(lmmp)$num_cells_expressed > 5,]
lmmp
genes_to_use <- rownames(lmmp)
#Remove any patterns learned on this dataset so we don't get confused
previous.pattern.ind <- stringr::str_detect(colnames(pData(lmmp)), pattern = "Pattern")
pData(lmmp) <- pData(lmmp)[,!previous.pattern.ind]

#Previously learned gene weights from NMF on logscaled expression counts
learned_weights <- read.csv(file = "./data/jared/lmmp_6mo_11-1_pattern_gene_weights.csv",
                            row.names = 1)
npattern <- dim(learned_weights)[2]



load.existing.projection.values <- 1

if(load.existing.projection.values){
  #first column is cell barcodes, rest of columns are patterns
  transferred_cell_wts <- read.csv("./data/jared/TC_LMMP_in_6mo_LMMP_11-25_patterns.csv", row.names = 1)
}else{
  set.seed(539)
  target_cell_wts_list <- lapply(1:npattern, function(i){
    wts <- as.matrix(learned_weights[,paste0("cellPattern",i), drop = F])
    target_cell_wts <- t(projectR::projectR(data = as.matrix(log10(exprs(lmmp)+1)),
                                            loadings = as.matrix(wts)))
  })

  transferred_cell_wts <- do.call(cbind, target_cell_wts_list)
  colnames(transferred_cell_wts) = paste0("cellPattern",1:npattern)
  transferred_cell_wts <- as.data.frame(transferred_cell_wts) %>% rownames_to_column(var = "cell_barcode")

  #write.csv(transferred_cell_wts,"/home/jared/projects/projection_drivers/ENS/results/TC_LMMP_in_6mo_LMMP_11-25_patterns.csv", quote = F)}
}

#Bind pattern weights to object, make sure cells are exactly the same. Could left_merge but pData not being a true data.frame is an issue
#assertthat::assert_that(sum(rownames(pData(lmmp)) == rownames(transferred_cell_wts)) == length(rownames(pData(lmmp))),
#                        msg = "Cells in projection weights matrix are not the same as cells in cds object")
#pData(lmmp) <- cbind(pData(lmmp), transferred_cell_wts)




learned_weights <- as.data.frame(learned_weights) %>% rownames_to_column(var = "gene_id")
learned_weights <- learned_weights[learned_weights[,"gene_id"] %in% genes_to_use,]


# Model fittings ----
genes_of_interest <- c("Hbb-bs", "Hba-a1","Hba-a2", "Malat1")
genes <- c(fData(lmmp)$gene_short_name[1:1], genes_of_interest)

#in the future iterate over patterns
pattern_of_interest <- c(37,38) #RBC

lmmp_sub <- lmmp[fData(lmmp)$gene_short_name %in% genes,]

projected_wts <- transferred_cell_wts[,paste0("cellPattern",pattern_of_interest)]
# "cell_type" is a columns in pData
model_str <- paste0("~0 + cell_type*patternWeight")
exp_family <- "negbinomial"
ncores = min(length(genes), 16)

system.time(glm_models <- fit_model_cds(cds = lmmp_sub, model_formula_str = model_str, projected_patterns = projected_wts,
                                        exp_family = exp_family, cores = ncores, clean_model = F, verbose = T))

#TODO: if a model fails, make sure the genes are labelled appropriately
#TODO: make wrapper function here?
extracted_models <- unlist(as.list(glm_models[,"model"]), recursive = F)
names(extracted_models) <- genes

#for each model (one gene one pattern) grab beta estimates of interest, p-values for it being non-zero, confidence intervals
extracted_coefficients <- lapply(extracted_models, FUN = extractCoefficients, params = c("cellPattern37", "cell_typeRBC", "cell_typeRBC:cellPattern37"))


#This is the way to implement outside of the monocle3 framework------
do.manual.fit <- 0
if(do.manual.fit){

  system.time(for(pattern in pattern_of_interest){
    model_str <- paste0("~0 + cell_type*cellPattern",pattern)
    glm.nb.out.list <- lapply(genes, fitHelper, lmmp, model_str)
    names(glm.nb.out.list) <- genes
  })


#Reverses the nesting, gives two main heading lists (result and error) which each have a list (or NULL) for every gene
glm.nb.out.list <- glm.nb.out.list %>% purrr::transpose()

#Either NULL or a successful model, full model returned from glm.nb()
successful_models <- pluck(glm.nb.out.list, "result")

#either NULL or an error message from glm.nb()
failed_models_msg <- pluck(glm.nb.out.list, "error")
}
###-----------





