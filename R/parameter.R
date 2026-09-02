#' Screening Step Parameters
#'
#' Creates a parameter object controlling the screening (recursive feature
#' elimination) step of \code{\link{mf}}.
#'
#' @param drop_fraction  Fraction of features dropped at each RFE iteration.
#'                       Default \code{0.25}.
#' @param keep_fraction  Fraction of features retained per module after
#'                       screening. Default \code{0.50}.
#' @param mtry_factor    Multiplier applied to \eqn{\sqrt{p}} to derive
#'                       \code{mtry} for each screening random forest.
#'                       Default \code{1}.
#' @param mtry_rule      How \code{mtry} is derived from the feature count
#'                       \code{p}. \code{"auto"} (default) uses \code{p/3} for
#'                       classification and \code{sqrt(p)} for regression (the
#'                       historical rule); \code{"sqrt"} uses \code{sqrt(p)} for
#'                       both; \code{"p_over_3"} uses \code{p/3} for both.
#' @param min_ntree      Minimum number of trees per screening forest.
#'                       Default \code{500}.
#' @param ntree_factor   Trees per forest = \code{max(p * ntree_factor, min_ntree)}.
#'                       Default \code{1}.
#'
#' @return An object of class \code{screen_control}.
#' @export
#' @seealso \code{\link{select_control}}, \code{\link{mf}}
screen_control <- function(drop_fraction = 0.25,
                            keep_fraction = 0.50,
                            mtry_factor   = 1,
                            mtry_rule     = "auto",
                            min_ntree     = 500L,
                            ntree_factor  = 1) {
  mtry_rule <- match.arg(mtry_rule, c("auto", "sqrt", "p_over_3"))
  obj <- list(
    drop_fraction = as.numeric(drop_fraction),
    keep_fraction = as.numeric(keep_fraction),
    mtry_factor   = as.numeric(mtry_factor),
    mtry_rule     = mtry_rule,
    min_ntree     = as.integer(min_ntree),
    ntree_factor  = as.numeric(ntree_factor)
  )
  class(obj) <- "screen_control"
  obj
}


#' Selection Step Parameters
#'
#' Creates a parameter object controlling the shadow-stability selection step
#' of \code{\link{mf}}.
#'
#' @param drop_fraction      Fraction of features dropped per RFE iteration
#'                           during survivor-pool pre-reduction. Default
#'                           \code{0.10}.
#' @param mtry_factor        \code{mtry} multiplier for shadow-stability random
#'                           forests. Default \code{1}.
#' @param min_ntree          Minimum trees per forest. Default \code{500}.
#' @param ntree_factor       Trees per forest scaling factor. Default \code{1}.
#' @param shadow_mode        Shadow stability pooling strategy. One of
#'                           \code{"split"} (correlated and unassigned pools
#'                           evaluated separately) or \code{"within_module"}
#'                           (one shadow run per module). Default \code{"split"}.
#' @param n_boots            Bootstrap replicates for shadow stability.
#'                           Default \code{50}.
#' @param pi_thr             Stability frequency threshold for correlated-pool
#'                           features. A feature is selected if its shadow-win
#'                           rate across bootstraps exceeds \code{pi_thr}.
#'                           Default \code{0.60}.
#' @param pi_thr_unassigned       Stability threshold for the unassigned pool.
#'                           \code{NULL} (default) uses an automatic elbow
#'                           detection rule.
#' @param shadow_percentile  Percentile of the shadow-feature distribution used
#'                           as the importance threshold within each bootstrap.
#'                           Default \code{95}.
#' @param unassigned_modules      Character vector of module names treated as
#'                           the unassigned pool (e.g. \code{"M10_unassigned"}).
#'                           \code{NULL} disables split-pool logic.
#' @param use_cfres   Logical. If \code{TRUE}, a cross-fitted random
#'                           forest is used to residualize \code{y} before the
#'                           unassigned-pool stability step, reducing confounding
#'                           from correlated modules. Default \code{FALSE}.
#' @param cfres_n_folds        Number of cross-fitting folds for cfRes
#'                           residualization. Default \code{5}.
#' @param cfres_ntree          Trees per fold-forest in cfRes residualization.
#'                           Default \code{300}.
#' @param cfres_scope          Scope of cfRes residualization. \code{"unassigned"}
#'                           (default) applies cfRes only to the unassigned pool.
#'                           \code{"all_modules"} applies per-module cfRes to
#'                           every module (removes the contribution of all other
#'                           modules before that module's stability run).
#'
#' @param mtry_rule        How \code{mtry} is derived from the feature count.
#'                         One of \code{"auto"} (default), \code{"sqrt"},
#'                         \code{"p_over_3"}. See \code{\link{screen_control}}.
#' @param mtry_on_real     Logical. If \code{TRUE}, \code{mtry} on the augmented
#'                         \eqn{[real | shadow]} matrix is \code{2 * rule(p_real)}
#'                         so real-feature sampling is not diluted by the shadows.
#'                         A statistical correction, not a speed-up.
#'                         Default \code{FALSE}.
#' @param threshold_mode   Selection cutoff rule. \code{"pi_thr"} (default) uses
#'                         the fixed \code{pi_thr} for the correlated pool and
#'                         elbow detection for the unassigned pool;
#'                         \code{"elbow"} applies the adaptive elbow rule to
#'                         \emph{both} pools.
#' @param early_stop_boots Logical. If \code{TRUE}, stop the bootstrap loop once
#'                         selection frequencies converge. Default \code{FALSE}.
#' @param early_stop_tol   Convergence tolerance for \code{early_stop_boots}:
#'                         stop when the largest change in any feature's
#'                         frequency between batches falls below this.
#'                         Default \code{0.01}.
#' @param early_stop_check_every Bootstraps per batch; convergence is checked
#'                         between batches. Default \code{10}.
#'
#' @return An object of class \code{select_control}.
#' @export
#' @seealso \code{\link{screen_control}}, \code{\link{mf}}
select_control <- function(drop_fraction    = 0.10,
                            mtry_factor      = 1,
                            mtry_rule        = "auto",
                            mtry_on_real     = FALSE,
                            min_ntree        = 500L,
                            ntree_factor     = 1,
                            shadow_mode      = "split",
                            n_boots          = 50L,
                            pi_thr           = 0.60,
                            pi_thr_unassigned     = NULL,
                            threshold_mode   = "pi_thr",
                            shadow_percentile = 95,
                            unassigned_modules    = NULL,
                            use_cfres = FALSE,
                            cfres_n_folds      = 5L,
                            cfres_ntree        = 300L,
                            cfres_scope        = "unassigned",
                            early_stop_boots = FALSE,
                            early_stop_tol   = 0.01,
                            early_stop_check_every = 10L) {

  shadow_mode    <- match.arg(shadow_mode,    c("split", "within_module"))
  cfres_scope      <- match.arg(cfres_scope,      c("unassigned", "all_modules"))
  mtry_rule      <- match.arg(mtry_rule,      c("auto", "sqrt", "p_over_3"))
  threshold_mode <- match.arg(threshold_mode, c("pi_thr", "elbow"))
  if (ntree_factor <= 0) stop("ntree_factor must be > 0", call. = FALSE)

  obj <- list(
    drop_fraction     = as.numeric(drop_fraction),
    mtry_factor       = as.numeric(mtry_factor),
    mtry_rule         = mtry_rule,
    mtry_on_real      = as.logical(mtry_on_real),
    min_ntree         = as.integer(min_ntree),
    ntree_factor      = as.numeric(ntree_factor),
    shadow_mode       = shadow_mode,
    n_boots           = as.integer(n_boots),
    pi_thr            = as.numeric(pi_thr),
    pi_thr_unassigned      = if (is.null(pi_thr_unassigned)) NULL else as.numeric(pi_thr_unassigned),
    threshold_mode    = threshold_mode,
    shadow_percentile = as.numeric(shadow_percentile),
    unassigned_modules     = if (is.null(unassigned_modules)) NULL else as.character(unassigned_modules),
    use_cfres  = as.logical(use_cfres),
    cfres_n_folds       = as.integer(cfres_n_folds),
    cfres_ntree         = as.integer(cfres_ntree),
    cfres_scope         = cfres_scope,
    early_stop_boots  = as.logical(early_stop_boots),
    early_stop_tol    = as.numeric(early_stop_tol),
    early_stop_check_every = as.integer(early_stop_check_every)
  )
  class(obj) <- "select_control"
  obj
}


#' WGCNA Parameter Organization
#'
#' Creates a parameter object for the WGCNA co-expression network step used
#' by \code{\link{wmf}}.
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
