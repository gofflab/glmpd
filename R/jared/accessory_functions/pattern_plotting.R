library(monocle3)
library(viridis)
library(ggplot2)
source("/home/jared/ENS/Timecourse_ENS/scripts/accessory_functions/monocle_mods.R")

#Define function to create a plot for a given pattern
plotCellPatterns = function (pattern_number, cds_obj, red_method = "UMAP", do.clip = NULL) {
  if(is.null(do.clip)){
    plot_cells_mod(cds_obj, reduction_method = red_method, color_cells_by = paste0("cellPattern",pattern_number), cell_size =1) 
    #  + ggtitle(paste0("Pattern ",pattern_number)) + 
    #    font("title", size= 40)
  }else{   
    plot_cells_mod(cds_obj, reduction_method = red_method, color_cells_by = paste0("cellPattern",pattern_number), cell_size =1) +
      scale_color_viridis(limits = quantile(pData(cds_obj)[,paste0("cellPattern",pattern_number)], do.clip) , oob = scales::squish,
                          guide_legend(title = paste0("cellPattern",pattern_number)))
  }
}

###Example usage - taken out of context might not run as is
#print to png to minimize file size. clip so that outlier cells do not dominate color scaling
if(0){
weighted_emb <- lapply(1:npattern,
                       FUN = plot_cell_patterns,
                       cds_obj = MENS,
                       red_method <- "UMAP",
                       do.clip <- c(0.02,.98))
png("/home/jared/ENS/Timecourse_ENS/plots/NMF/MENS/MENS_A_NMF_Patterns.png", height= 2000, width = 2000)
print(do.call(ggarrange, weighted_emb[1:length(condition_patterns)]))
dev.off()
}

#Function to plot a continuous variable (eg pattern usage) over a binned variable such as time.
#provide feature (continuous) and bin_by/color_by (categorical, discrete) as a string that matches a column in pData.
plotPatternUsageByCondition <- function(feature, cds, conditions,  bin_by, color_by){
  ggplot(as.data.frame(pData(cds)), aes_string(x = bin_by, y = feature)) +
    geom_boxplot() +
    ggtitle(feature) +
    xlab("cell pattern weight")
}

###Example usage - taken out of context might not run as is
if(0){
  condition_patterns<- as.list(1:npattern)
  myplots <- lapply(paste0("cellPattern",condition_patterns), 
                    FUN= plotPatternUsageByCondition,
                    cds = cds,
                    bin_by = "age")
  pdf("/home/jared/ENS/Timecourse_ENS/plots/NMF/MENS/MENS_patterns_over_age.pdf", height= 20, width = 20)
  print(do.call(ggarrange, myplots[1:length(condition_patterns)]))
  dev.off()
}