#####################################################
# Visualizations around model fitting and summaries
#####################################################


#####################################################
# Visualizations for projected pattern weights
#####################################################


#TODO:Jared usually calls with this lapply but that is not user friendly, should make the default
#TODO:Pass through arguments with ...

#' plotCellPatterns
#' This function provides a wrapper around monocle3's plotting function to allow for clipping of color scales, to limit the effect of outliers on visualization.
#' Implemented to plot continuous variables in pData that share common prefixes and unique appended ids... eg "pattern1" "pattern2" "patternZ"
#' @param pattern_suffix A vector of names (or unique parts of names - used with pattern_prefix) of continuous values in pData to color cells by
#' @param cds A monocle3 cell_data_set object
#' @param pattern_prefix Used with pattern_suffix, if all variables to be plotted have a common prefix. only implemented like this to allow for lapply
#' @param red_method Passed through monocle::plot_cells
#' @param do.clip A vector of two numbers on the scale [0,1] indicating minimum and maximum percentiles at which to clip the color scale.
#'
#'
#' @return Returns a plot (or list of plots) with cells colored by pattern weight
#' @import monocle3
#' @import viridis
#' @export
#'
# #' @examples
plotCellPatterns = function (pattern_suffix, cds, pattern_prefix=NULL, red_method = "UMAP", do.clip = c(0,1)) {
  plot_cells(cds, reduction_method = red_method, color_cells_by = paste0(cell_prefix,pattern_suffix), cell_size =1) +
    scale_color_viridis(limits = quantile(pData(cds)[,paste0(pattern_prefix,pattern_suffix)], do.clip) , oob = scales::squish,
                        guide_legend(title = paste0(pattern_prefix,pattern_suffix)))

}


#TODO:Jared usually calls with this lapply but that is not user friendly, should make the default
#TODO:Pass through arguments with ...

#' #Function to plot a continuous variable (eg pattern usage) over a categorical variable such as time.
#' Provide feature (continuous) and bin_by (categorical, discrete) as strings that matches a column in pData.
#' @param feature A names or vector of names of continuous values in pData to plot the values
#' @param cds A monocle3 cell_data_set object
#' @param bin_by A categorical value in pData to group cells by
#'
#' @return Returns a boxplot (or list of plots) of cells projected pattern weights, for groups of interest
#' @import monocle3
#' @import ggplot2
#' @export
#'
# #' @examples

plotPatternUsageByCondition <- function(feature, cds, conditions,  bin_by = NULL){

  ggplot(as.data.frame(pData(cds)), aes_string(x = bin_by, y = feature)) +
    geom_boxplot() +
    ggtitle(feature) +
    xlab(feature)
}


#TODO: Add Roxygen2 parameters
#Extract one parameter for all patterns for heatmap plotting
#coefficient_list is the output from orderCoefficientByGene
coefficientHeatmap <- function(coefficients_list, coefficient_to_plot){

  pattern_full_names <- names(coefficients_list)
  param_df <- lapply(pattern_full_names, function(pattern_full_name){
    #pull each gene name and the coefficient to be plotted
    param <- coefficients_list[[pattern_full_name]][,c("gene_short_name",coefficient_to_plot)]
    #indexing isn't very elegant
    colnames(param)[2] <- paste(colnames(param)[2],pattern_full_name, sep = "-")
    return(param)
  })

  #combine each list (;attern) to one data frame. Columns are coefficients (one per pattern), rows are gene
  param_merged_df <- param_df %>% purrr::reduce(left_join, by = "gene_short_name")

  #Create heatmap for all genes, all patterns, one parameter
  param_mat <- as.matrix(param_merged_df %>% tibble::column_to_rownames(var = "gene_short_name"))
  pl <- ComplexHeatmap::Heatmap(param_mat, name = paste0("beta coefficient for ", coefficient_to_plot))

  return(pl)
}


