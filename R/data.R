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

#' ALL Leukemia Gene Expression Data
#'
#' A list containing the top 500 most-variable probe-sets from the Acute
#' Lymphoblastic Leukemia (ALL) microarray dataset (Chiaretti et al., 2004),
#' along with the corresponding cell-lineage outcome.
#'
#' The full dataset (128 patients, 12 625 probe-sets) was obtained from the
#' \code{ALL} Bioconductor package.  Features were ranked by across-sample
#' variance and the top 500 retained.  See \code{data-raw/prep_leukemia.R}
#' for the exact preprocessing steps.
#'
#' @docType data
#' @keywords datasets
#' @name leukemia
#' @usage data(leukemia)
#' @format A list with two elements:
#' \describe{
#'   \item{X}{A numeric matrix with 128 rows (patients) and 500 columns
#'     (Affymetrix HGU95Av2 probe-sets, already log\eqn{_2}-scaled by the
#'     original authors).}
#'   \item{y}{An integer vector of length 128: \code{1} (B-cell lineage) or
#'     \code{0} (T-cell lineage), derived from the \code{BT} column of the
#'     original ALL phenotype data.}
#' }
#' @references Chiaretti, S. et al. (2004). Gene expression profile of adult
#'   T-cell acute lymphocytic leukemia identifies distinct subsets of patients
#'   with different response to therapy and survival. \emph{Blood}, 103(7),
#'   2771–2778. \doi{10.1182/blood-2003-09-3243}
NULL
