library(monocle3)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(MASS)
source("/home/jared/projects/projection_drivers/ENS/scripts/accessory_functions/fitHelper.R")

# Load in required datsets and pweights --------------------------
#dataset to project into.
lmmp <- readRDS("/home/jared/ENS/Timecourse_ENS/TC_LMMP.rds")
#lmmp <- lmmp[fData(lmmp)$num_cells_expressed > 5,]
lmmp
genes_to_use <- rownames(lmmp)
#Remove any patterns learned on this dataset so we don't get confused
previous.pattern.ind <- stringr::str_detect(colnames(pData(lmmp)), pattern = "Pattern")
pData(lmmp) <- pData(lmmp)[,!previous.pattern.ind]

#Previously learned gene weights from NMF on logscaled expression counts
learned_weights <- read.csv(file = "/home/jared/ENS/6mo_LMMP/results/NMF/lmmp/lmmp_6mo_11-1_pattern_gene_weights.csv",
                            row.names = 1)
npattern <- dim(learned_weights)[2]



load.existing.projection.values <- 1

if(load.existing.projection.values){
  #first column is cell barcodes, rest of columns are patterns
  transferred_cell_wts <- read.csv("/home/jared/projects/projection_drivers/ENS/results/TC_LMMP_in_6mo_LMMP_11-25_patterns.csv", row.names = 1)
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

#Bind pattern weights to object
pData(lmmp) <- cbind(pData(lmmp), transferred_cell_wts)




learned_weights <- as.data.frame(learned_weights) %>% rownames_to_column(var = "gene_id")
learned_weights <- learned_weights[learned_weights[,"gene_id"] %in% genes_to_use,]


# Model fittings ----
genes_of_interest <- c("Hbb-bs", "Hba-a1","Hba-a2", "Malat1")
genes <- c(fData(lmmp)$gene_short_name[1:1], genes_of_interest)

pattern_of_interest <- c(37) #RBC 

lmmp_sub <- lmmp[fData(lmmp)$gene_short_name %in% genes,]

model_str <- paste0("~0 + cell_type*cellPattern",pattern_of_interest)
exp_family <- "negbinomial"
ncores = min(length(genes), 16)

system.time(glm_models <- fit_models(cds = lmmp_sub, model_formula_str = model_str, expression_family = exp_family, cores = ncores,
                                     clean_model = F, verbose = T))

#TODO: if a model fails, make sure the genes are labelled appropriately
extracted_models <- unlist(as.list(glm_models[,"model"]), recursive = F)
names(extracted_models) <- genes

do.manual.fit <- 0
if(do.manual.fit){
##TODO: Loop or lapply this over patterns (mapply?)
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

