#' Screening Step Parameters
#'
#' Creates a parameter object controlling the screening (recursive feature
#' elimination) step of \code{\link{sf}}.
#'
#' @param drop_fraction  Fraction of features dropped at each RFE iteration.
#'                       Default \code{0.25}.
#' @param keep_fraction  Fraction of features retained per module after
#'                       screening. Default \code{0.50}.
#' @param mtry_factor    Multiplier applied to \eqn{\sqrt{p}} to derive
#'                       \code{mtry} for each screening random forest.
#'                       Default \code{1}.
#' @param min_ntree      Minimum number of trees per screening forest.
#'                       Default \code{500}.
#' @param ntree_factor   Trees per forest = \code{max(p * ntree_factor, min_ntree)}.
#'                       Default \code{1}.
#'
#' @return An object of class \code{screen_control}.
#' @export
#' @seealso \code{\link{select_control}}, \code{\link{sf}}
screen_control <- function(drop_fraction = 0.25,
                            keep_fraction = 0.50,
                            mtry_factor   = 1,
                            min_ntree     = 500L,
                            ntree_factor  = 1L) {
  obj <- list(
    drop_fraction = as.numeric(drop_fraction),
    keep_fraction = as.numeric(keep_fraction),
    mtry_factor   = as.numeric(mtry_factor),
    min_ntree     = as.integer(min_ntree),
    ntree_factor  = as.integer(ntree_factor)
  )
  class(obj) <- "screen_control"
  obj
}


#' Selection Step Parameters
#'
#' Creates a parameter object controlling the shadow-stability selection step
#' of \code{\link{sf}}.
#'
#' @param drop_fraction      Fraction of features dropped per RFE iteration
#'                           during survivor-pool pre-reduction. Default
#'                           \code{0.10}.
#' @param mtry_factor        \code{mtry} multiplier for shadow-stability random
#'                           forests. Default \code{1}.
#' @param min_ntree          Minimum trees per forest. Default \code{500}.
#' @param ntree_factor       Trees per forest scaling factor. Default \code{1}.
#' @param shadow_mode        Shadow stability pooling strategy. One of
#'                           \code{"split"} (correlated and independent pools
#'                           evaluated separately) or \code{"within_module"}
#'                           (one shadow run per module). Default \code{"split"}.
#' @param n_boots            Bootstrap replicates for shadow stability.
#'                           Default \code{50}.
#' @param pi_thr             Stability frequency threshold for correlated-pool
#'                           features. A feature is selected if its shadow-win
#'                           rate across bootstraps exceeds \code{pi_thr}.
#'                           Default \code{0.60}.
#' @param pi_thr_indep       Stability threshold for the independent pool.
#'                           \code{NULL} (default) uses an automatic elbow
#'                           detection rule.
#' @param shadow_percentile  Percentile of the shadow-feature distribution used
#'                           as the importance threshold within each bootstrap.
#'                           Default \code{95}.
#' @param indep_modules      Character vector of module names treated as
#'                           the independent pool (e.g. \code{"M10_indep"}).
#'                           \code{NULL} disables split-pool logic.
#' @param use_dml_residual   Logical. If \code{TRUE}, a cross-fitted random
#'                           forest is used to residualize \code{y} before the
#'                           independent-pool stability step, reducing confounding
#'                           from correlated modules. Default \code{FALSE}.
#' @param dml_n_folds        Number of cross-fitting folds for DML
#'                           residualization. Default \code{5}.
#' @param dml_ntree          Trees per fold-forest in DML residualization.
#'                           Default \code{300}.
#' @param dml_scope          Scope of DML residualization. \code{"indep"}
#'                           (default) applies DML only to the independent pool.
#'                           \code{"all_modules"} applies per-module DML to
#'                           every module (removes the contribution of all other
#'                           modules before that module's stability run).
#'
#' @return An object of class \code{select_control}.
#' @export
#' @seealso \code{\link{screen_control}}, \code{\link{sf}}
select_control <- function(drop_fraction    = 0.10,
                            mtry_factor      = 1,
                            min_ntree        = 500L,
                            ntree_factor     = 1L,
                            shadow_mode      = "split",
                            n_boots          = 50L,
                            pi_thr           = 0.60,
                            pi_thr_indep     = NULL,
                            shadow_percentile = 95,
                            indep_modules    = NULL,
                            use_dml_residual = FALSE,
                            dml_n_folds      = 5L,
                            dml_ntree        = 300L,
                            dml_scope        = "indep") {

  shadow_mode <- match.arg(shadow_mode, c("split", "within_module"))
  dml_scope   <- match.arg(dml_scope,   c("indep", "all_modules"))

  obj <- list(
    drop_fraction     = as.numeric(drop_fraction),
    number_selected   = 1L,          # placeholder; no K-trim in current version
    mtry_factor       = as.numeric(mtry_factor),
    min_ntree         = as.integer(min_ntree),
    ntree_factor      = as.integer(ntree_factor),
    n_boots           = as.integer(n_boots),
    pi_thr            = as.numeric(pi_thr),
    pi_thr_indep      = if (is.null(pi_thr_indep)) NULL else as.numeric(pi_thr_indep),
    evidence_min_pct  = 50,
    shadow_mode       = shadow_mode,
    shadow_percentile = as.numeric(shadow_percentile),
    indep_modules     = if (is.null(indep_modules)) NULL else as.character(indep_modules),
    use_dml_residual  = as.logical(use_dml_residual),
    dml_n_folds       = as.integer(dml_n_folds),
    dml_ntree         = as.integer(dml_ntree),
    dml_scope         = dml_scope,
    # internal/legacy — kept for interface compatibility
    split_mode           = "k",
    split_pct_select     = 0.25,
    c_corr               = 3L,
    use_rfe_trim         = TRUE,
    null_floor_adjust    = FALSE,
    k_corr_pool          = NULL,
    k_ind_pool           = NULL,
    ind_selection_mode   = "shadow",
    dual_pass2           = FALSE,
    k_ind_pass2          = NULL,
    ind_bypass_pass2     = FALSE,
    use_interventional_shap = FALSE,
    interventional_bg_size  = 100L,
    percent_flag         = FALSE
  )
  class(obj) <- "select_control"
  obj
}


#' WGCNA Parameter Organization
#'
#' Creates a parameter object for the WGCNA co-expression network step used
#' by \code{\link{wsf}}.
#'
#' @param power        Soft-thresholding power passed to
#'                     \code{\link[WGCNA]{blockwiseModules}}. Default \code{6}.
#' @param min_features Minimum number of features per module. Modules smaller
#'                     than this trigger a warning. Default \code{20}.
#' @param ...          Additional arguments forwarded to
#'                     \code{\link[WGCNA]{blockwiseModules}}.
#'
#' @return An object of class \code{WGCNA_control}.
#' @export
WGCNA_control <- function(power = 6, min_features = 20, ...) {
  obj <- list(
    power        = power,
    min_features = min_features,
    extra_args   = list(...)
  )
  class(obj) <- "WGCNA_control"
  obj
}
