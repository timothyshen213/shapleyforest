#' Shapley Forest Feature Selection
#'
#' Runs the Shapley Forest algorithm: module-stratified recursive feature
#' elimination (screening), followed by shadow-stability selection, and a
#' single-shot SHAP importance pass on the stable feature set.
#'
#' @export
#' @param X                 A data frame where each column is a feature.
#' @param y                 Response vector (numeric for regression; factor for
#'                          classification).
#' @param module_membership A named character vector assigning each feature in
#'                          \code{X} to a module. Names must match
#'                          \code{colnames(X)}.
#' @param screen_params     Screening parameters. See
#'                          \code{\link{screen_control}}.
#' @param select_params     Selection parameters. See
#'                          \code{\link{select_control}}.
#' @param final_ntree       Number of trees in the final random forest.
#'                          Default \code{500}.
#' @param nodesize          Minimum terminal node size. Defaults to \code{1}
#'                          for classification and \code{5} for regression.
#' @param num_processors    Number of parallel threads for random forests.
#'                          Default \code{1}.
#' @param max_modsize       Maximum module size before blockwise pre-reduction
#'                          is applied (Python backend only). Default \code{5000}.
#' @param auto_initial      Controls the initial screening output.
#'                          \code{"1"} — save CSV and stop;
#'                          \code{"2"} — save CSV and continue;
#'                          \code{"3"} — stop without saving;
#'                          \code{"4"} — continue without saving (default).
#' @param verbose           Verbosity level. \code{0} = silent; \code{1} =
#'                          progress messages; \code{2} = detailed output.
#'                          Default \code{1}.
#' @param seed              Integer random seed for reproducibility. Defaults to
#'                          current system time.
#' @param backend           Computation backend. \code{"python"} (default) uses
#'                          a Python scikit-learn / SHAP engine via
#'                          \pkg{reticulate} and provides exact TreeSHAP values
#'                          plus an interaction matrix.  \code{"R"} uses a
#'                          pure-R \pkg{ranger} implementation — no Python
#'                          installation required. Choose the final importance
#'                          method with \code{r_shap}.
#' @param r_shap            Final SHAP / importance method used when
#'                          \code{backend = "R"}. Ignored for the Python
#'                          backend. One of:
#'   \describe{
#'     \item{\code{"permutation"}}{Ranger permutation variable importance
#'       (default). No extra packages needed. Returns a scalar VIM per
#'       feature; per-observation SHAP values and SHAP-based plots are
#'       unavailable.}
#'     \item{\code{"fastshap"}}{Approximate Shapley values via
#'       \pkg{fastshap} (\code{nsim = 50} permutations). Requires the
#'       \pkg{fastshap} package. No interaction matrix.}
#'     \item{\code{"treeshap"}}{Exact TreeSHAP values via the R
#'       \pkg{treeshap} package. Requires \pkg{treeshap}. Also computes a
#'       \eqn{p \times p} interaction matrix, enabling
#'       \code{detect_interaction()} and
#'       \code{plot_potential_interactions()}.}
#'   }
#'
#' @return An object of class \code{\link{shapley_forest}}.
#'
#' @references
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#'   Newey, W., & Robins, J. (2018). Double/debiased machine learning for
#'   treatment and structural parameters. \emph{The Econometrics Journal},
#'   21(1), C1–C68.
#'
#' Lundberg, S. M., & Lee, S. I. (2017). A unified model for interpreting
#'   predictions. \emph{NeurIPS}, 30.
#'
#' @importFrom ranger ranger
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   sf_setup(condaenv = "sfenv")
#'
#'   data(iris)
#'   X   <- iris[, 1:4]
#'   y   <- iris$Species
#'   mem <- setNames(c("A","A","B","B"), colnames(X))
#'
#'   res <- sf(X, y, module_membership = mem,
#'             screen_params = screen_control(min_ntree = 100),
#'             select_params = select_control(n_boots = 25))
#'   print(res)
#'   plot_importance(res)
#' }
sf <- function(X, y,
               module_membership,
               screen_params   = screen_control(),
               select_params   = select_control(),
               final_ntree     = 500L,
               nodesize        = NULL,
               num_processors  = 1L,
               max_modsize     = 5000L,
               auto_initial    = "4",
               verbose         = 1L,
               seed            = as.integer(Sys.time()),
               backend         = "python",
               r_shap          = "permutation") {

  backend <- match.arg(backend, c("python", "R"))
  r_shap  <- match.arg(r_shap,  c("permutation", "fastshap", "treeshap"))

  # ── lazy-load Python backend ─────────────────────────────────────────────────
  if (backend == "python") .ensure_python()

  set.seed(seed)

  # ── validation ───────────────────────────────────────────────────────────────
  if (!is.data.frame(X))
    stop("X must be a data.frame.", call. = FALSE)
  if (!(is.vector(y) || is.factor(y)))
    stop("y must be a numeric vector (regression) or factor (classification).",
         call. = FALSE)
  if (is.null(names(module_membership)))
    names(module_membership) <- colnames(X)[seq_along(module_membership)]

  CLASSIFICATION <- is.factor(y)
  if (is.null(nodesize))
    nodesize <- if (CLASSIFICATION) 1L else 5L

  if (verbose == 0L) options(warn = -1)

  sc <- screen_params
  sl <- select_params

  module_membership_screen <- module_membership[names(module_membership) %in% colnames(X)]
  module_list_screen       <- unique(module_membership_screen)

  runtime <- list(Screen = NA_real_, Selection = NA_real_, Final_RF = NA_real_)

  # ══════════════════════════════════════════════════════════════════════════════
  # R backend — pure ranger, no Python required
  # ══════════════════════════════════════════════════════════════════════════════
  if (backend == "R") {
    return(.sf_R_backend(
      X = X, y = y,
      module_membership        = module_membership,
      module_membership_screen = module_membership_screen,
      module_list_screen       = module_list_screen,
      screen_params            = screen_params,
      select_params            = select_params,
      final_ntree              = final_ntree,
      nodesize                 = nodesize,
      num_processors           = num_processors,
      auto_initial             = auto_initial,
      verbose                  = verbose,
      seed                     = seed,
      CLASSIFICATION           = CLASSIFICATION,
      runtime                  = runtime,
      r_shap                   = r_shap
    ))
  }

  # ── prepare inputs for Python ────────────────────────────────────────────────
  X_mat         <- as.matrix(X); mode(X_mat) <- "double"
  y_py          <- if (CLASSIFICATION) as.integer(y) else as.numeric(y)
  feature_names <- colnames(X)
  module_assign <- as.character(module_membership_screen[feature_names])
  module_list_r <- as.character(module_list_screen)

  # ── Phase 1: Screening ───────────────────────────────────────────────────────
  if (verbose >= 1L) cat("Screening ...\n")
  t0 <- proc.time()

  screen_py <- .sf_py$run_screen_rfe(
    X_mat              = X_mat,
    y                  = y_py,
    feature_names      = feature_names,
    module_assignments = module_assign,
    module_list        = module_list_r,
    drop_fraction      = sc$drop_fraction,
    keep_fraction      = sc$keep_fraction,
    mtry_factor        = sc$mtry_factor,
    ntree_factor       = as.numeric(sc$ntree_factor),
    min_ntree          = as.integer(sc$min_ntree),
    mtry_rule          = sc$mtry_rule,
    nodesize           = as.integer(nodesize),
    classification     = CLASSIFICATION,
    seed               = as.integer(seed),
    n_jobs             = as.integer(num_processors),
    max_modsize        = as.integer(max_modsize),
    verbose            = as.logical(verbose >= 2L)
  )

  runtime$Screen <- (proc.time() - t0)["elapsed"]

  # ── Reconstruct survivor_results from Python output ──────────────────────────
  survivor_results <- vector("list", length(module_list_screen))
  for (i in seq_along(module_list_screen)) {
    mod_name <- as.character(module_list_screen[i])
    feats    <- as.character(screen_py$survivor_features[[i]])
    vims     <- as.numeric(screen_py$survivor_vims[[i]])
    surv_df  <- data.frame(features    = feats,
                            mod_varlist = vims,
                            stringsAsFactors = FALSE)
    mask <- (screen_py$initial_module == mod_name)
    init_df <- if (any(mask)) {
      data.frame(
        "Module"   = screen_py$initial_module[mask],
        "Feature"  = screen_py$initial_feat[mask],
        "VIM"      = as.numeric(screen_py$initial_vim[mask]),
        "Survivor" = as.logical(screen_py$initial_survivor[mask]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    } else NULL
    survivor_results[[i]] <- list(survivor = surv_df, screen_df = init_df)
  }
  initial_screen <- do.call(rbind, lapply(survivor_results, `[[`, "screen_df"))

  # ── Initial screening output ──────────────────────────────────────────────────
  if (!is.null(auto_initial) && auto_initial %in% c("1", "2")) {
    write.csv(initial_screen, "initial_screen.csv", row.names = FALSE)
    if (verbose >= 1L) cat("Initial screen saved to 'initial_screen.csv'.\n")
    if (auto_initial == "1")
      stop("Execution stopped after initial screen (auto_initial='1').", call. = FALSE)
  }
  if (!is.null(auto_initial) && auto_initial == "3")
    stop("Execution stopped after initial screen (auto_initial='3').", call. = FALSE)

  # ── Assemble survivor pool ────────────────────────────────────────────────────
  survivor_list  <- lapply(survivor_results, `[[`, "survivor")
  names(survivor_list) <- module_list_screen

  survivors         <- do.call(rbind, survivor_list)
  survivors         <- as.data.frame(survivors, stringsAsFactors = FALSE)
  survivors[, 2]    <- as.numeric(survivors[, 2])
  names(survivors)  <- c("featureID", "Permutation VIM")
  X_surv <- X[, names(X) %in% survivors[, 1], drop = FALSE]

  # ── Phase 2: Shadow-stability selection ──────────────────────────────────────
  if (verbose >= 1L) cat("Selecting stable features ...\n")
  t0 <- proc.time()

  X_surv_mat       <- as.matrix(X_surv); mode(X_surv_mat) <- "double"

  surv_feat_to_mod <- unlist(lapply(names(survivor_list), function(mod) {
    df <- survivor_list[[mod]]
    if (is.null(df) || nrow(df) == 0L) return(character(0))
    setNames(rep(mod, nrow(df)), as.character(df[, 1L]))
  }))
  module_assigns_ordered <- surv_feat_to_mod[colnames(X_surv)]
  module_assigns_ordered[is.na(module_assigns_ordered)] <- "unknown"

  module_ev <- vapply(names(survivor_list), function(mod) {
    df <- survivor_list[[mod]]
    if (is.null(df) || nrow(df) == 0L) return(0.0)
    mean(as.numeric(df[, 2L]), na.rm = TRUE)
  }, numeric(1L))

  select_py <- .sf_py$run_select_rfe(
    X_surv_mat              = X_surv_mat,
    y                       = y_py,
    feature_names_surv      = as.list(colnames(X_surv)),
    module_assignments_surv = as.list(as.character(module_assigns_ordered)),
    drop_fraction           = sl$drop_fraction,
    number_selected         = 1L,
    mtry_factor             = sl$mtry_factor,
    ntree_factor            = as.numeric(sl$ntree_factor),
    min_ntree               = as.integer(sl$min_ntree),
    nodesize                = as.integer(nodesize),
    classification          = CLASSIFICATION,
    seed                    = as.integer(seed),
    n_jobs                  = as.integer(num_processors),
    max_modsize             = as.integer(max_modsize),
    n_boots                 = as.integer(sl$n_boots),
    pi_thr                  = as.numeric(sl$pi_thr),
    shadow_mode             = sl$shadow_mode,
    shadow_percentile       = as.numeric(sl$shadow_percentile),
    pi_thr_indep            = if (!is.null(sl$pi_thr_indep))
                                as.numeric(sl$pi_thr_indep) else NULL,
    threshold_mode          = sl$threshold_mode,
    mtry_rule               = sl$mtry_rule,
    mtry_on_real            = as.logical(sl$mtry_on_real),
    early_stop_boots        = as.logical(sl$early_stop_boots),
    early_stop_tol          = as.numeric(sl$early_stop_tol),
    early_stop_check_every  = as.integer(sl$early_stop_check_every),
    indep_modules           = if (!is.null(sl$indep_modules))
                                as.list(sl$indep_modules) else NULL,
    use_dml_residual        = as.logical(sl$use_dml_residual),
    dml_n_folds             = as.integer(sl$dml_n_folds),
    dml_ntree               = as.integer(sl$dml_ntree),
    dml_scope               = sl$dml_scope,
    verbose                 = as.logical(verbose >= 2L)
  )

  runtime$Selection <- (proc.time() - t0)["elapsed"]

  # ── Unpack Python selection output ───────────────────────────────────────────
  stable_features <- as.character(select_py$final_features)
  n_stable        <- length(stable_features)

  if (verbose >= 1L)
    cat(sprintf("Stable features: %d\n", n_stable))

  stability_freq_vec <- if (n_stable > 0L)
    setNames(as.numeric(select_py$stability_freq), stable_features)
  else
    setNames(numeric(0), character(0))

  shap_matrix <- if (n_stable > 0L && length(select_py$shap_matrix) > 0L) {
    matrix(
      as.numeric(unlist(select_py$shap_matrix)),
      nrow     = nrow(X_surv),
      ncol     = n_stable,
      byrow    = TRUE,   # Python .tolist() is row-major; fill by row to preserve
                         # per-observation alignment (else the matrix is scrambled)
      dimnames = list(NULL, stable_features)
    )
  } else {
    matrix(NA_real_, nrow = nrow(X_surv), ncol = 0L)
  }

  # ── TreeSHAP interaction matrix (p × p mean |SHAP interaction|) ──────────────
  interaction_matrix <- if (n_stable > 0L && length(select_py$interaction_matrix) > 0L) {
    m <- matrix(
      as.numeric(unlist(select_py$interaction_matrix)),
      nrow     = n_stable,
      ncol     = n_stable,
      byrow    = TRUE,   # row-major from Python (symmetric, but keep consistent)
      dimnames = list(stable_features, stable_features)
    )
    m
  } else {
    NULL
  }

  selection_list <- lapply(seq_along(select_py$selection_list_feats), function(i) {
    data.frame(
      feature_name        = as.character(select_py$selection_list_feats[[i]]),
      variable_importance = round(as.numeric(select_py$selection_list_vims[[i]]), 4L),
      stringsAsFactors    = FALSE
    )
  })

  # ── Per-pool shadow-stability data (for plot_stability_elbow) ───────────────
  stability_data <- if (!is.null(select_py$pool_stability) &&
                         length(select_py$pool_stability) > 0L) {
    pool_rows <- lapply(select_py$pool_stability, function(p) {
      n_p <- length(p$feature_names)
      if (n_p == 0L) return(NULL)
      er  <- if (!is.null(p$elbow_rank)) as.integer(p$elbow_rank) else NA_integer_
      data.frame(
        pool        = rep(as.character(p$pool), n_p),
        rank        = seq_len(n_p),
        feature     = as.character(p$feature_names),
        freq        = as.numeric(p$freqs),
        selected    = as.logical(p$selected),
        threshold   = rep(as.numeric(p$threshold), n_p),
        elbow_rank  = rep(er, n_p),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), pool_rows))
  } else {
    NULL
  }

  # ── Phase 3: Final ranger RF ─────────────────────────────────────────────────
  if (verbose >= 1L) cat("Fitting final model ...\n")
  t0 <- proc.time()

  feature_list <- data.frame(
    feature_name        = stable_features,
    variable_importance = round(as.numeric(select_py$final_vims), 4L),
    module_membership   = as.character(module_membership[stable_features]),
    stringsAsFactors    = FALSE
  )

  final_X    <- X[, stable_features, drop = FALSE]
  current_p  <- ncol(final_X)

  if (current_p == 0L) {
    # empty stable set: return a valid (predictor-less) fit instead of letting
    # ranger error on a zero-covariate design
    if (verbose >= 1L)
      cat("No stable features selected; skipping final model.\n")
    final_rf <- NULL
  } else {
  final_mtry <- {
    base <- switch(sl$mtry_rule,
      sqrt     = sqrt(current_p),
      p_over_3 = current_p / 3,
      if (CLASSIFICATION) current_p / 3 else sqrt(current_p))   # "auto"
    max(1L, min(ceiling(sl$mtry_factor * base), current_p))
  }

  final_rf <- ranger::ranger(
    x             = final_X,
    y             = y,
    mtry          = final_mtry,
    num.trees     = as.integer(final_ntree),
    importance    = "none",
    min.node.size = as.integer(nodesize),
    probability   = CLASSIFICATION,
    verbose       = FALSE,
    num.threads   = as.integer(num_processors),
    seed          = as.integer(seed)
  )
  }

  final_module_membership <- data.frame(
    feature_name = names(final_X),
    module       = as.character(module_membership[names(final_X)]),
    stringsAsFactors = FALSE
  )

  # ── Build final_SHAP data frame (from Python TreeSHAP) ───────────────────────
  mean_shap <- setNames(as.numeric(select_py$final_vims), stable_features)

  shap_final_list <- data.frame(
    feature_name        = stable_features,
    variable_importance = round(mean_shap, 4L),
    module_membership   = as.character(module_membership[stable_features]),
    stringsAsFactors    = FALSE
  )

  # Add per-observation SHAP CIs if matrix is available
  if (ncol(shap_matrix) > 0L) {
    feat_ord  <- stable_features
    valid_cols <- feat_ord[feat_ord %in% colnames(shap_matrix)]
    shap_ord  <- shap_matrix[, valid_cols, drop = FALSE]
    n_obs     <- nrow(shap_ord)

    abs_mat      <- abs(shap_ord)
    mean_abs     <- colMeans(abs_mat)
    se_abs       <- apply(abs_mat,  2L, sd) / sqrt(n_obs)
    mean_signed  <- colMeans(shap_ord)
    se_signed    <- apply(shap_ord, 2L, sd) / sqrt(n_obs)

    shap_final_list$mean_abs_shap    <- round(mean_abs[feat_ord],                              4L)
    shap_final_list$se_abs_shap      <- round(se_abs[feat_ord],                                4L)
    shap_final_list$ci_lower_abs     <- round(pmax(0, mean_abs[feat_ord] - 1.96 * se_abs[feat_ord]), 4L)
    shap_final_list$ci_upper_abs     <- round(mean_abs[feat_ord] + 1.96 * se_abs[feat_ord],   4L)
    shap_final_list$mean_signed_shap <- round(mean_signed[feat_ord],                           4L)
    shap_final_list$se_signed_shap   <- round(se_signed[feat_ord],                             4L)
    shap_final_list$ci_lower_signed  <- round(mean_signed[feat_ord] - 1.96 * se_signed[feat_ord], 4L)
    shap_final_list$ci_upper_signed  <- round(mean_signed[feat_ord] + 1.96 * se_signed[feat_ord], 4L)
    shap_final_list$stability_freq   <- round(stability_freq_vec[feat_ord],                    4L)

    # ── global effect direction: Spearman(feature value, SHAP value) ──────────
    # Mean signed SHAP ~ 0 (contributions cancel); direction lives in how a
    # feature's SHAP tracks its value. Spearman is rank-based (captures monotone
    # direction without assuming linearity). |rho| < 0.1 flags non-monotone /
    # ambiguous direction (a single sign is misleading — use a dependence plot).
    val_ord  <- final_X[, valid_cols, drop = FALSE]
    dir_corr <- vapply(seq_along(valid_cols), function(k) {
      xj <- as.numeric(val_ord[[k]]); sj <- shap_ord[, k]
      if (stats::sd(xj) == 0 || stats::sd(sj) == 0) return(0)
      cc <- suppressWarnings(stats::cor(xj, sj, method = "spearman"))
      if (is.na(cc)) 0 else cc
    }, numeric(1L))
    names(dir_corr) <- valid_cols
    direction <- sign(dir_corr)
    shap_final_list$direction           <- as.integer(direction[feat_ord])
    shap_final_list$dir_corr            <- round(dir_corr[feat_ord], 4L)
    shap_final_list$signed_importance   <- round(direction[feat_ord] * mean_abs[feat_ord], 4L)
    shap_final_list$direction_ambiguous <- abs(dir_corr[feat_ord]) < 0.1
  }

  runtime$Final_RF <- (proc.time() - t0)["elapsed"]

  # shap_obj: $shapley_values = signed n×p SHAP matrix
  #           $interaction_matrix = p×p mean |SHAP interaction| (TreeSHAP exact)
  shap_obj <- list(
    shapley_values    = shap_matrix,
    interaction_matrix = interaction_matrix
  )
  class(shap_obj) <- "explanation"

  if (verbose >= 1L) cat("Done.\n")
  options(warn = 1)

  shapley_forest(
    final_rf          = final_rf,
    final_X           = final_X,
    module_membership = final_module_membership,
    WGCNA_object      = NULL,
    survivor_list     = survivor_list,
    selection_list    = selection_list,
    final_shap        = shap_final_list,
    shap_obj          = shap_obj,
    runtime           = runtime,
    feature_list      = feature_list,
    stability_freq    = stability_freq_vec,
    stability_data    = stability_data
  )
}


#' Shapley Forest with WGCNA Module Detection
#'
#' Runs a Weighted Gene Co-expression Network Analysis (WGCNA) to derive
#' module membership, then passes the result to \code{\link{sf}}.
#'
#' @inheritParams sf
#' @param WGCNA_params WGCNA parameters. See \code{\link{WGCNA_control}}.
#'
#' @return An object of class \code{\link{shapley_forest}} with an additional
#'   \code{$WGCNA_object} slot containing the raw WGCNA output.
#'
#' @export
#' @importFrom ranger ranger
#'
#' @references
#' Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
#'   co-expression network analysis. \emph{Statistical Applications in Genetics
#'   and Molecular Biology}, 4(1), Article 17.
wsf <- function(X, y,
                WGCNA_params  = WGCNA_control(p = 6),
                screen_params = screen_control(),
                select_params = select_control(),
                final_ntree   = 500L,
                nodesize      = NULL,
                num_processors = 1L,
                max_modsize   = 5000L,
                auto_initial  = "4",
                verbose       = 1L,
                seed          = as.integer(Sys.time()),
                backend       = "python",
                r_shap        = "permutation") {

  if (!requireNamespace("WGCNA", quietly = TRUE))
    stop("Package 'WGCNA' is required for wsf(). Install it first.", call. = FALSE)
  # WGCNA's blockwiseModules calls cor() with extra args (weights.x, cosine, ...)
  # that only WGCNA::cor understands.  When WGCNA is only namespace-loaded
  # (not library()-attached) the lookup resolves to stats::cor and errors.
  # attachNamespace() puts WGCNA on the search path without startup noise.
  if (!"WGCNA" %in% .packages())
    attachNamespace("WGCNA")

  if (!(is.vector(y) || is.factor(y)))
    stop("y must be a numeric vector or factor.", call. = FALSE)

  # Convert X to numeric for WGCNA
  integer_cols <- sapply(X, is.integer)
  if (any(integer_cols))
    X[, integer_cols] <- lapply(X[, integer_cols, drop = FALSE], as.numeric)
  if (!all(sapply(X, is.numeric)))
    stop("All columns of X must be numeric for WGCNA.", call. = FALSE)

  # Run WGCNA
  wgcna_ctrl <- WGCNA_params
  WGCNA_args <- c(list(datExpr = X, power = wgcna_ctrl$power),
                  wgcna_ctrl$extra_args)
  if (verbose == 0L) {
    invisible(capture.output(suppressWarnings(suppressMessages({
      bwise <- do.call(WGCNA::blockwiseModules, WGCNA_args)
    }))))
  } else {
    bwise <- do.call(WGCNA::blockwiseModules, WGCNA_args)
  }
  module_membership <- bwise$colors

  out <- sf(
    X = X, y = y,
    module_membership = module_membership,
    screen_params     = screen_params,
    select_params     = select_params,
    final_ntree       = final_ntree,
    nodesize          = nodesize,
    num_processors    = num_processors,
    max_modsize       = max_modsize,
    auto_initial      = auto_initial,
    verbose           = verbose,
    seed              = seed,
    backend           = backend,
    r_shap            = r_shap
  )

  out$WGCNA_object <- bwise
  out
}
