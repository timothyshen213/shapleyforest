#' Female Mice Liver Expression Data
#'
#' This dataset includes gene expression levels from liver tissue
#' of female mice. The data originates from research conducted by 
#' Ghazalpour et al (\url{https://doi.org/10.1186/1471-2105-9-559}).
#' It is also widely used as an WGCNA tutorial. However, at this moment,
#' the documentations for WGCNA no longer exist. See the original
#' WGCNA paper for more information 
#' (\url{https://doi.org/10.1186/1471-2105-9-559}).
#'
#' @docType data
#' @keywords datasets
#' @name femData
#' @usage data(femData)
#' @format A data frame with 66 rows and 3600 columns
NULL

#' Female Mice Clinical Trait Data
#'
#' This dataset goes a long with the Female Mice Liver Expression
#' Data. It contains clinical trait measurements of female mice.
#' It shares the same observations as Female Mice Liver Expression
#' Data and provides possible response variables for liver expression.
#' 
#' @docType data
#' @keywords datasets
#' @name traitData
#' @usage data(traitData)
#' @format A data frame with 38 rows and 361 columns
NULL

#' Cleaned and Preprocess Female Mice Liver Expression Data
#' 
#' Removes outlier genes through WGCNA preprocessing of the Female
#' Liver Expression Data. Weight (g) variable from Clinical Trait
#' Data is also included. See Female Mice Liver Expression Vignette
#' for more information.
#'
#' @docType data
#' @keywords datasets
#' @name Liver_Exp
#' @usage data(Liver_Exp)
#' @format A data frame with 66 rows and 3601 columns
NULL
