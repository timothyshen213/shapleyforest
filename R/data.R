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
#' along with the corresponding molecular-subtype outcome.
#'
#' Starting from the \code{ALL} Bioconductor package (128 patients, 12 625
#' probe-sets), samples are restricted to the B-lineage cases carrying one of
#' two balanced molecular subtypes — the BCR/ABL (Philadelphia-chromosome, Ph+)
#' fusion versus cytogenetically normal (NEG) — yielding 79 patients (37
#' BCR/ABL, 42 NEG).  This is the standard *balanced* ALL benchmark used in the
#' Bioconductor "Case Studies", replacing the heavily imbalanced B-vs-T lineage
#' split.  Features were then ranked by across-sample variance and the top 500
#' retained.  See \code{data-raw/prep_leukemia.R} for the exact preprocessing
#' steps.
#'
#' @docType data
#' @keywords datasets
#' @name leukemia
#' @usage data(leukemia)
#' @format A list with two elements:
#' \describe{
#'   \item{X}{A numeric matrix with 79 rows (patients) and 500 columns
#'     (Affymetrix HGU95Av2 probe-sets, already log\eqn{_2}-scaled by the
#'     original authors).}
#'   \item{y}{An integer vector of length 79: \code{1} (BCR/ABL subtype) or
#'     \code{0} (NEG, cytogenetically normal), derived from the
#'     \code{mol.biol} column of the original ALL phenotype data.}
#' }
#' @references Chiaretti, S. et al. (2004). Gene expression profile of adult
#'   T-cell acute lymphocytic leukemia identifies distinct subsets of patients
#'   with different response to therapy and survival. \emph{Blood}, 103(7),
#'   2771–2778. \doi{10.1182/blood-2003-09-3243}
NULL

#' World Development Indicators (2020)
#'
#' A list containing World Bank development indicators and life expectancy for
#' 181 countries in 2020.  Unlike the other bundled datasets — which are all
#' gene-expression panels — this is a non-biological, genuinely
#' high-dimensional (p > n), and highly-correlated tabular dataset.  It is
#' included to demonstrate Bonsai Forest on data whose features (GDP per
#' capita, literacy rate, CO2 emissions, urban population share, ...) are
#' intuitive to a general audience, while still posing the same challenge as
#' genomic data: many more correlated, partially redundant measurements than
#' observations competing to explain the outcome.
#'
#' @docType data
#' @keywords datasets
#' @name wdi
#' @usage data(wdi)
#' @format A list with three elements:
#' \describe{
#'   \item{X}{A numeric matrix with 181 rows (countries) and 501 columns —
#'     World Bank development indicators for 2020, filtered to those present
#'     for at least 80\% of countries and with residual gaps median-imputed.
#'     Row names are country short names.}
#'   \item{y}{A numeric vector of length 181: life expectancy at birth
#'     (total, in years).  All life-table-derived series (life expectancy,
#'     mortality, survival, and death counts) were removed from \code{X} to
#'     avoid tautological prediction.}
#'   \item{labels}{A named character vector mapping each indicator code in
#'     \code{colnames(X)} to its human-readable World Bank indicator name.}
#' }
#' @source The World Development Indicators bulk database
#'   (\url{https://databank.worldbank.org/}).  See \code{data-raw/prep_wdi.R}
#'   for the exact download and preprocessing steps.
#' @references World Bank. \emph{World Development Indicators}. Washington, DC:
#'   The World Bank Group.
NULL
