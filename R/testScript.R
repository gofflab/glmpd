# Test the getProjectionDrivers.R functions using data from the ALS analysis

#source('getProjectionDrivers.R')
#source('coeff_table_glmpd.R')
source('model_fitting.R')

# Directories
dataDir <- '../inst/extdata/'

# load the counts mat, sample annotations, and projected p weights
countsMat <- readRDS(paste(dataDir,'humanCountsMat.rds',sep=''))
annotDF <- readRDS(paste(dataDir,'humanAnnotations.rds',sep=''))
POIs <- readRDS(paste(dataDir,'sexPatts_0_9thresh.rds',sep='')) # patterns of interest

# limit to genes expressed in at least one sample
totCounts <- apply(countsMat,1,sum)
countsMat_exp <- countsMat[totCounts>0,]

# order factors to obtain desired baseline categories in annotations map
annotDF$group <- as.factor(as.character(annotDF$group))
annotDF$group <- factor(annotDF$group,levels=levels(annotDF$group)[c(2,1,3)])

# start by trying a couple patterns and a couple genes
POIs_test <- POIs[7:8,] # 7 is female, 8 is male
countsMat_exp_test <- as.matrix(countsMat_exp[c("XIST","IL22","TAGLN"),])

# designate model parameters
model_formula_str <- "~  patternWeights*sex + tissueType + patternWeights:tissueType + group + patternWeights:group"

#fits <- getProjectionDrivers(countsMat_exp_test,annotDF,POIs_test,model_formula_str)
fits <- fitGLMpd(countsMat_exp_test,annotDF=annotDF,model_formula_str=model_formula_str,projected_patterns = t(POIs_test))
coeffTables <- coeff_table_glmpd(fits)
