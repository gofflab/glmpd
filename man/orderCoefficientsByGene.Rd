% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fitting.R
\name{orderCoefficientsByGene}
\alias{orderCoefficientsByGene}
\title{orderCoefficientsByGene
Order model coefficient estimates by genes, such as for visualization by heatmap}
\usage{
orderCoefficientsByGene(
  pattern_coefficient_list,
  model_terms_to_keep = NULL,
  filter_significance,
  string
)
}
\arguments{
\item{pattern_coefficient_list}{List (over patterns) of objects returned from extractCoefficients(), or "summarized_result" option of fit_models_cds().}

\item{model_terms_to_keep}{optional. limit returned dataframe to coefficients of interest.}

\item{filter_significance}{optional. Numeric value at which to filter q-value significance.
Genes with signficant coefficients for one or more  terms are returned. If not provided, all genes returned.}

\item{string}{optional. all model terms that contain this string will be used for signficance filtering.}
}
\value{
Returns a list of fitted models (output similar to monocle3::fit_model)
}
\description{
orderCoefficientsByGene
Order model coefficient estimates by genes, such as for visualization by heatmap
}
