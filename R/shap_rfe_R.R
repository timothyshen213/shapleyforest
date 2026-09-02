# Pure-R backend for mossyforest — no Python / reticulate required.
#
# Screening RFE and shadow-stability selection both use ranger + permutation
# importance.  For the final importance step the user chooses one of three
# methods via the `r_shap` argument passed down from mf():
#
#   "permutation" — ranger permutation VIM (scalar per feature, no extra pkg)
#   "fastshap"    — approximate Shapley values via the fastshap package
#                   (interaction_matrix = NULL)
#   "treeshap"    — exact TreeSHAP via the treeshap package; also computes an
#                   interaction matrix identical in structure to the Python
#                   backend's (diagonal = main effects, off-diagonal = pairwise)
#
# Internal helpers — not exported.


# ── Top-level orchestrator called by mf() ────────────────────────────────────

.mf_R_backend <- function(X, y,
                           module_membership,
                           module_membership_screen,
                           module_list_screen,
                           screen_params,
                           select_params,
                           final_ntree,
                           nodesize,
                           num_processors,
                           auto_initial,
                           verbose,
                           seed,
                           CLASSIFICATION,
                           runtime,
                           r_shap) {

  sc <- screen_params
  sl <- select_params
  if (is.null(nodesize)) nodesize <- if (CLASSIFICATION) 1L else 5L
  if (verbose == 0L) options(warn = -1)

  # ── Phase 1: Screening ───────────────────────────────────────────────────────
  if (verbose >= 1L) cat(sprintf("Screening [R / %s] ...\n", r_shap))
  t0 <- proc.time()

  screen_out <- .mf_R_screen(
    X                 = X,
    y                 = y,
    module_membership = module_membership_screen,
    screen_params     = sc,
    nodesize          = nodesize,
    num_processors    = num_processors,
    seed              = seed,
    verbose           = verbose >= 2L,
    r_shap            = r_shap
  )

  runtime$Screen <- (proc.time() - t0)["elapsed"]

  survivor_list  <- screen_out$survivor_list
  initial_screen <- screen_out$initial_screen

  # ── Initial screening output ─────────────────────────────────────────────────
  if (!is.null(auto_initial) && auto_initial %in% c("1", "2")) {
    utils::write.csv(initial_screen, "initial_screen.csv", row.names = FALSE)
    if (verbose >= 1L) cat("Initial screen saved to 'initial_screen.csv'.\n")
    if (auto_initial == "1")
      stop("Execution stopped after initial screen (auto_initial='1').", call. = FALSE)
  }
  if (!is.null(auto_initial) && auto_initial == "3")
    stop("Execution stopped after initial screen (auto_initial='3').", call. = FALSE)

  # ── Assemble survivor pool ───────────────────────────────────────────────────
  survivors_df        <- do.call(rbind, survivor_list)
  survivors_df        <- as.data.frame(survivors_df, stringsAsFactors = FALSE)
  survivors_df[, 2L]  <- as.numeric(survivors_df[, 2L])
  names(survivors_df) <- c("featureID", "Permutation VIM")
  X_surv <- X[, names(X) %in% survivors_df[, 1L], drop = FALSE]

  surv_feat_to_mod <- unlist(lapply(names(survivor_list), function(mod) {
    df <- survivor_list[[mod]]
    if (is.null(df) || nrow(df) == 0L) return(character(0L))
    setNames(rep(mod, nrow(df)), as.character(df[, 1L]))
  }))
  module_assigns_surv <- surv_feat_to_mod[colnames(X_surv)]
  module_assigns_surv[is.na(module_assigns_surv)] <- "unknown"

  # ── Phase 2: Shadow-stability selection ─────────────────────────────────────
  if (verbose >= 1L) cat("Selecting stable features [R] ...\n")
  t0 <- proc.time()

  sel_out <- .mf_R_select(
    X_surv                 = X_surv,
    y                      = y,
    survivor_list          = survivor_list,
    select_params          = sl,
    module_membership_surv = module_assigns_surv,
    nodesize               = nodesize,
    num_processors         = num_processors,
    seed                   = seed,
    verbose                = verbose >= 1L,
    r_shap                 = r_shap
  )

  runtime$Selection <- (proc.time() - t0)["elapsed"]

  stable_features    <- sel_out$stable_features
  n_stable           <- length(stable_features)
  stability_freq_vec <- sel_out$stability_freq
  selection_list_r   <- sel_out$selection_list

  if (verbose >= 1L) cat(sprintf("Stable features: %d\n", n_stable))

  # ── Build stability_data (same format as Python backend) ────────────────────
  stability_data <- if (length(sel_out$pool_stability) > 0L) {
    pool_rows <- lapply(sel_out$pool_stability, function(ps) {
      np <- length(ps$feature_names)
      if (np == 0L) return(NULL)
      er <- if (is.null(ps$elbow_rank) || is.na(ps$elbow_rank))
              NA_integer_ else as.integer(ps$elbow_rank)
      data.frame(
        pool       = rep(as.character(ps$pool), np),
        rank       = seq_len(np),
        feature    = as.character(ps$feature_names),
        freq       = as.numeric(ps$freqs),
        selected   = as.logical(ps$selected),
        threshold  = rep(as.numeric(ps$threshold), np),
        elbow_rank = rep(er, np),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), pool_rows))
  } else {
    NULL
  }

  # ── Phase 3: Final ranger RF ─────────────────────────────────────────────────
  if (verbose >= 1L)
    cat(sprintf("Fitting final model [R / %s] ...\n", r_shap))
  t0 <- proc.time()

  final_X   <- if (n_stable > 0L)
    X[, stable_features, drop = FALSE]
  else
    X[, character(0L), drop = FALSE]

  current_p  <- ncol(final_X)
  final_mtry <- if (current_p > 0L) {
    .compute_mtry_R(current_p, sl$mtry_factor, CLASSIFICATION, sl$mtry_rule)
  } else 1L

  # Always fit with permutation importance so variable.importance is available
  final_rf <- if (current_p > 0L) {
    ranger::ranger(
      x             = final_X,
      y             = y,
      mtry          = final_mtry,
      num.trees     = as.integer(final_ntree),
      importance    = "permutation",
      min.node.size = as.integer(nodesize),
      probability   = CLASSIFICATION,
      verbose       = FALSE,
      num.threads   = as.integer(num_processors),
      seed          = as.integer(seed)
    )
  } else {
    NULL
  }

  # ── Compute SHAP values (method depends on r_shap) ───────────────────────────
  shap_result <- if (!is.null(final_rf) && n_stable > 0L) {
    switch(r_shap,
      permutation = .mf_R_shap_permutation(final_rf, final_X, stable_features,
                                            verbose >= 1L),
      fastshap    = .mf_R_shap_fastshap(final_rf, final_X, y, stable_features,
                                         seed, verbose >= 1L),
      treeshap    = .mf_R_shap_treeshap(final_rf, final_X, stable_features,
                                         verbose >= 1L)
    )
  } else {
    list(
      shap_matrix        = matrix(NA_real_, nrow = nrow(X), ncol = 0L),
      interaction_matrix = NULL
    )
  }

  shap_matrix        <- shap_result$shap_matrix
  interaction_matrix <- shap_result$interaction_matrix

  # ── Importance scores (VIM or mean|SHAP|) ───────────────────────────────────
  # For permutation: use ranger's permutation VIM (scalar per feature).
  # For fastshap / treeshap: use mean |SHAP| per feature (matches Python output).
  perm_vims <- if (!is.null(final_rf) && n_stable > 0L)
    final_rf$variable.importance[stable_features]
  else
    setNames(rep(NA_real_, n_stable), stable_features)

  use_shap_vim <- r_shap %in% c("fastshap", "treeshap") &&
    !all(is.na(shap_matrix)) && ncol(shap_matrix) > 0L

  vim_for_report <- if (use_shap_vim)
    setNames(colMeans(abs(shap_matrix), na.rm = TRUE), stable_features)
  else
    as.numeric(perm_vims)

  # ── feature_list & final_SHAP data frame ────────────────────────────────────
  feature_list <- data.frame(
    feature_name        = stable_features,
    variable_importance = round(vim_for_report, 4L),
    module_membership   = as.character(module_membership[stable_features]),
    stringsAsFactors    = FALSE
  )

  final_module_membership <- data.frame(
    feature_name = stable_features,
    module       = as.character(module_membership[stable_features]),
    stringsAsFactors = FALSE
  )

  shap_final_list <- data.frame(
    feature_name        = stable_features,
    variable_importance = round(vim_for_report, 4L),
    module_membership   = as.character(module_membership[stable_features]),
    stringsAsFactors    = FALSE
  )

  # Add per-observation SHAP statistics when a SHAP matrix is available
  if (n_stable > 0L && !all(is.na(shap_matrix)) && ncol(shap_matrix) > 0L) {
    abs_mat     <- abs(shap_matrix)
    n_obs       <- nrow(shap_matrix)
    mean_abs    <- colMeans(abs_mat,     na.rm = TRUE)
    se_abs      <- apply(abs_mat,  2L, stats::sd, na.rm = TRUE) / sqrt(n_obs)
    mean_signed <- colMeans(shap_matrix, na.rm = TRUE)
    se_signed   <- apply(shap_matrix, 2L, stats::sd, na.rm = TRUE) / sqrt(n_obs)

    shap_final_list$mean_abs_shap    <- round(mean_abs[stable_features],  4L)
    shap_final_list$se_abs_shap      <- round(se_abs[stable_features],    4L)
    shap_final_list$ci_lower_abs     <- round(
      pmax(0, mean_abs[stable_features] - 1.96 * se_abs[stable_features]), 4L)
    shap_final_list$ci_upper_abs     <- round(
      mean_abs[stable_features] + 1.96 * se_abs[stable_features], 4L)
    shap_final_list$mean_signed_shap <- round(mean_signed[stable_features], 4L)
    shap_final_list$se_signed_shap   <- round(se_signed[stable_features],   4L)
    shap_final_list$ci_lower_signed  <- round(
      mean_signed[stable_features] - 1.96 * se_signed[stable_features], 4L)
    shap_final_list$ci_upper_signed  <- round(
      mean_signed[stable_features] + 1.96 * se_signed[stable_features], 4L)
    shap_final_list$stability_freq   <- round(
      stability_freq_vec[stable_features], 4L)

    # ── global effect direction: Spearman(feature value, SHAP value) ──────────
    # Only meaningful with a real per-observation SHAP matrix (fastshap /
    # treeshap); permutation VIM is scalar-only, so there is nothing to
    # correlate against and direction is skipped for r_shap = "permutation".
    # Mean signed SHAP ~ 0 (contributions cancel); direction lives in how a
    # feature's SHAP tracks its value. Spearman is rank-based (captures
    # monotone direction without assuming linearity). |rho| < 0.1 flags
    # non-monotone / ambiguous direction (a single sign is misleading).
    #
    # Classification caveat: for treeshap this classification path currently
    # errors upstream (treeshap::ranger.unify() does not support ranger
    # probability forests), so direction is only reachable for regression
    # today. For fastshap, the SHAP matrix explains the probability of the
    # LOWER-sorted factor level (ranger's predictions column 1) — the
    # opposite convention from the Python backend, which explains the
    # HIGHER-sorted class. Direction sign is therefore not comparable
    # between backends for fastshap + classification.
    if (r_shap %in% c("fastshap", "treeshap")) {
      dir_corr <- vapply(stable_features, function(feat) {
        xj <- as.numeric(final_X[[feat]]); sj <- shap_matrix[, feat]
        if (anyNA(sj) || stats::sd(xj) == 0 || stats::sd(sj) == 0) return(0)
        cc <- suppressWarnings(stats::cor(xj, sj, method = "spearman"))
        if (is.na(cc)) 0 else cc
      }, numeric(1L))
      direction <- sign(dir_corr)
      shap_final_list$direction           <- as.integer(direction[stable_features])
      shap_final_list$dir_corr            <- round(dir_corr[stable_features], 4L)
      shap_final_list$signed_importance   <- round(
        direction[stable_features] * mean_abs[stable_features], 4L)
      shap_final_list$direction_ambiguous <- abs(dir_corr[stable_features]) < 0.1
    }
  }

  runtime$Final_RF <- (proc.time() - t0)["elapsed"]

  shap_obj <- list(
    shapley_values     = shap_matrix,
    interaction_matrix = interaction_matrix
  )
  class(shap_obj) <- "explanation"

  if (verbose >= 1L) cat("Done.\n")
  options(warn = 1)

  # ── Rename survivor_list columns to standard names ───────────────────────────
  survivor_list_named <- setNames(
    lapply(names(survivor_list), function(mod) {
      df <- survivor_list[[mod]]
      if (is.null(df) || nrow(df) == 0L) return(df)
      names(df) <- c("featureID", "Permutation VIM")
      df
    }),
    names(survivor_list)
  )

  mossy_forest(
    final_rf          = final_rf,
    final_X           = final_X,
    module_membership = final_module_membership,
    WGCNA_object      = NULL,
    survivor_list     = survivor_list_named,
    selection_list    = selection_list_r,
    final_shap        = shap_final_list,
    shap_obj          = shap_obj,
    runtime           = runtime,
    feature_list      = feature_list,
    stability_freq    = stability_freq_vec,
    stability_data    = stability_data
  )
}


# ── SHAP method 1: ranger permutation VIM ────────────────────────────────────
# Permutation importance is a single scalar per feature, not per-observation.
# We return an NA shap_matrix so that plot_importance / SHAP-based plots are
# simply not available, while the variable_importance column in the output
# data frame carries the permutation VIM.

.mf_R_shap_permutation <- function(final_rf, final_X, stable_features, verbose) {
  if (verbose) cat("  Using permutation VIM (no per-observation SHAP values).\n")
  list(
    shap_matrix        = matrix(NA_real_, nrow = nrow(final_X),
                                ncol     = length(stable_features),
                                dimnames = list(NULL, stable_features)),
    interaction_matrix = NULL
  )
}


# ── SHAP method 2: fastshap approximate SHAP ─────────────────────────────────

.mf_R_shap_fastshap <- function(final_rf, final_X, y, stable_features,
                                 seed, verbose) {
  if (!requireNamespace("fastshap", quietly = TRUE))
    stop("Install the 'fastshap' package to use r_shap = 'fastshap'.",
         call. = FALSE)

  if (verbose) cat("  Computing SHAP values via fastshap (approximate) ...\n")

  CLASSIFICATION <- is.factor(y)
  pred_fn <- if (CLASSIFICATION)
    function(object, newdata)
      predict(object, data = newdata)$predictions[, 1L, drop = TRUE]
  else
    function(object, newdata)
      predict(object, data = newdata)$predictions

  set.seed(seed)
  shap_mat <- tryCatch(
    fastshap::explain(
      object        = final_rf,
      feature_names = stable_features,
      X             = final_X,
      pred_wrapper  = pred_fn,
      nsim          = 50L
    ),
    error = function(e) {
      message("fastshap failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(shap_mat)) {
    message("Falling back to NA SHAP matrix. Check fastshap installation.")
    return(list(
      shap_matrix        = matrix(NA_real_, nrow = nrow(final_X),
                                  ncol     = length(stable_features),
                                  dimnames = list(NULL, stable_features)),
      interaction_matrix = NULL
    ))
  }

  out <- as.matrix(shap_mat)
  colnames(out) <- stable_features
  list(shap_matrix = out, interaction_matrix = NULL)
}


# ── SHAP method 3: treeshap exact TreeSHAP + interaction matrix ──────────────

.mf_R_shap_treeshap <- function(final_rf, final_X, stable_features, verbose) {
  if (!requireNamespace("treeshap", quietly = TRUE))
    stop("Install the 'treeshap' package to use r_shap = 'treeshap'.",
         call. = FALSE)

  if (verbose)
    cat("  Computing exact TreeSHAP values via treeshap (+ interaction matrix) ...\n")

  unified <- tryCatch(
    treeshap::ranger.unify(final_rf, final_X),
    error = function(e)
      stop("treeshap::ranger.unify() failed: ", conditionMessage(e), call. = FALSE)
  )

  ts <- tryCatch(
    treeshap::treeshap(unified, final_X, interactions = TRUE, verbose = FALSE),
    error = function(e) {
      message("treeshap with interactions failed, retrying without: ",
              conditionMessage(e))
      tryCatch(
        treeshap::treeshap(unified, final_X, interactions = FALSE, verbose = FALSE),
        error = function(e2)
          stop("treeshap::treeshap() failed: ", conditionMessage(e2), call. = FALSE)
      )
    }
  )

  # SHAP matrix (n × p) — treeshap names columns by feature name already
  shap_raw    <- as.matrix(ts$shaps)
  shap_cols   <- colnames(shap_raw)

  # Subset / reorder to exactly stable_features (treeshap may add bias col etc.)
  keep <- intersect(stable_features, shap_cols)
  shap_matrix <- shap_raw[, keep, drop = FALSE]
  # If treeshap dropped some features, pad with NA columns
  missing <- setdiff(stable_features, keep)
  if (length(missing) > 0L) {
    pad <- matrix(NA_real_, nrow = nrow(shap_matrix), ncol = length(missing),
                  dimnames = list(NULL, missing))
    shap_matrix <- cbind(shap_matrix, pad)
  }
  shap_matrix <- shap_matrix[, stable_features, drop = FALSE]

  # Interaction matrix (p × p) — mean |SHAP interaction| across observations.
  # treeshap stores interactions as [p, p, n] (features × features × obs).
  interaction_matrix <- if (!is.null(ts$interactions)) {
    int_arr <- ts$interactions          # [p, p, n]
    int_mat <- apply(abs(int_arr), c(1L, 2L), mean, na.rm = TRUE)   # [p, p]
    p_dim   <- ncol(shap_matrix)
    if (nrow(int_mat) == p_dim && ncol(int_mat) == p_dim)
      rownames(int_mat) <- colnames(int_mat) <- stable_features
    int_mat
  } else {
    NULL
  }

  list(shap_matrix = shap_matrix, interaction_matrix = interaction_matrix)
}


# ── Elbow detection ───────────────────────────────────────────────────────────

.find_elbow_R <- function(freqs) {
  n <- length(freqs)
  if (n <= 2L) return(1L)
  # Perpendicular distance from each point to the line joining (1, f1)-(n, fn)
  p1x <- 1.0;        p1y <- freqs[1L]
  p2x <- as.double(n); p2y <- freqs[n]
  dx  <- p2x - p1x
  dy  <- p2y - p1y
  len <- sqrt(dx^2 + dy^2)
  if (len < .Machine$double.eps) return(1L)
  dist <- abs(dy * seq_len(n) - dx * freqs + p2x * p1y - p2y * p1x) / len
  which.max(dist)
}


# ── mtry rule (mirrors the Python package / backend) ─────────────────────────
.compute_mtry_R <- function(n_features, mtry_factor, classification,
                             rule = "auto") {
  base <- switch(rule,
    sqrt     = sqrt(n_features),
    p_over_3 = n_features / 3,
    if (classification) n_features / 3 else sqrt(n_features))   # "auto"
  max(1L, min(ceiling(mtry_factor * base), n_features))
}

# mtry on the augmented [real | shadow] matrix (2 * p_real cols)
.compute_mtry_aug_R <- function(p_real, mtry_factor, classification,
                                 rule = "auto", on_real = FALSE) {
  n_aug <- 2L * p_real
  if (isTRUE(on_real))
    return(max(1L, min(2L * .compute_mtry_R(p_real, mtry_factor, classification, rule), n_aug)))
  .compute_mtry_R(n_aug, mtry_factor, classification, rule)
}

.apply_stable_threshold_R <- function(names_sorted, freqs_sorted,
                                       pi_thr, use_elbow = FALSE) {
  if (length(names_sorted) == 0L)
    return(list(stable    = character(0L),
                threshold = 0.0,
                elbow_rank = NA_integer_))

  freqs <- as.numeric(freqs_sorted)

  if (use_elbow) {
    ei     <- .find_elbow_R(freqs)
    ef     <- freqs[ei]
    stable <- names_sorted[freqs >= ef & freqs > 0.0]
    return(list(stable     = stable,
                threshold  = ef,
                elbow_rank = as.integer(ei)))
  } else {
    stable <- names_sorted[freqs >= as.numeric(pi_thr)]
    return(list(stable     = stable,
                threshold  = as.numeric(pi_thr),
                elbow_rank = NA_integer_))
  }
}


# ── Phase 1: Screening (R backend) ───────────────────────────────────────────

# ── TreeSHAP importance for a ranger forest (mean |SHAP| per feature) ─────────
# Used when r_shap = "treeshap" so SCREENING and STABILITY SELECTION rank by
# TreeSHAP (mirroring the Python engine) instead of ranger permutation VIM.
# Requires the 'treeshap' package; supports regression + binary classification.
.vim_treeshap_ranger <- function(rf, Xdata, feats) {
  if (!requireNamespace("treeshap", quietly = TRUE))
    stop("Install the 'treeshap' package to use r_shap = 'treeshap' for ",
         "screening / stability selection.", call. = FALSE)
  unified <- treeshap::ranger.unify(rf, Xdata)
  ts      <- treeshap::treeshap(unified, Xdata, interactions = FALSE, verbose = FALSE)
  ma      <- colMeans(abs(as.matrix(ts$shaps)), na.rm = TRUE)
  out     <- setNames(rep(0.0, length(feats)), feats)
  common  <- intersect(feats, names(ma))
  out[common] <- ma[common]
  out
}


.mf_R_screen <- function(X, y, module_membership, screen_params,
                          nodesize, num_processors, seed, verbose,
                          r_shap = "permutation") {
  sc             <- screen_params
  CLASSIFICATION <- is.factor(y)
  use_ts         <- identical(r_shap, "treeshap")
  module_list    <- unique(module_membership)
  # NOTE: unlike the Python backend, the pure-R backend applies NO blockwise
  # pre-reduction (max_modsize is ignored); modules are screened whole.

  survivor_list <- vector("list", length(module_list))
  names(survivor_list) <- module_list

  for (mod in module_list) {
    feats  <- names(module_membership)[module_membership == mod]
    X_mod  <- X[, feats, drop = FALSE]
    n_mod  <- length(feats)

    if (verbose) cat(sprintf("  Screening module '%s' (%d features) ...\n", mod, n_mod))

    target_keep   <- max(1L, ceiling(sc$keep_fraction * n_mod))
    current_feats <- feats

    # RFE loop
    while (length(current_feats) > target_keep) {
      cp      <- length(current_feats)
      n_trees <- max(as.integer(sc$min_ntree),
                     as.integer(sc$ntree_factor * cp))
      mtry_r  <- .compute_mtry_R(cp, sc$mtry_factor, CLASSIFICATION, sc$mtry_rule)

      rf_step <- ranger::ranger(
        x             = X_mod[, current_feats, drop = FALSE],
        y             = y,
        mtry          = mtry_r,
        num.trees     = n_trees,
        importance    = "permutation",
        min.node.size = as.integer(nodesize),
        probability   = CLASSIFICATION,
        verbose       = FALSE,
        num.threads   = as.integer(num_processors),
        seed          = as.integer(seed)
      )

      vims   <- if (use_ts)
        .vim_treeshap_ranger(rf_step, X_mod[, current_feats, drop = FALSE], current_feats)
      else rf_step$variable.importance
      n_drop <- max(1L, ceiling(sc$drop_fraction * cp))
      drop_i <- order(vims)[seq_len(n_drop)]
      current_feats <- setdiff(current_feats, current_feats[drop_i])
    }

    # Final importance fit on survivors
    cp      <- length(current_feats)
    n_trees <- max(as.integer(sc$min_ntree),
                   as.integer(sc$ntree_factor * cp))
    mtry_r  <- .compute_mtry_R(cp, sc$mtry_factor, CLASSIFICATION, sc$mtry_rule)

    rf_fin <- ranger::ranger(
      x             = X_mod[, current_feats, drop = FALSE],
      y             = y,
      mtry          = mtry_r,
      num.trees     = n_trees,
      importance    = "permutation",
      min.node.size = as.integer(nodesize),
      probability   = CLASSIFICATION,
      verbose       = FALSE,
      num.threads   = as.integer(num_processors),
      seed          = as.integer(seed)
    )

    vims_fin <- if (use_ts)
      .vim_treeshap_ranger(rf_fin, X_mod[, current_feats, drop = FALSE], current_feats)
    else rf_fin$variable.importance

    surv_df <- data.frame(
      features    = current_feats,
      mod_varlist = as.numeric(vims_fin[current_feats]),
      stringsAsFactors = FALSE
    )

    screen_df <- data.frame(
      Module   = rep(mod, n_mod),
      Feature  = feats,
      VIM      = ifelse(feats %in% current_feats,
                        as.numeric(vims_fin[feats]), NA_real_),
      Survivor = feats %in% current_feats,
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )

    survivor_list[[mod]] <- list(survivor = surv_df, screen_df = screen_df)
  }

  initial_screen <- do.call(rbind, lapply(survivor_list, `[[`, "screen_df"))
  survivor_dfs   <- lapply(survivor_list, `[[`, "survivor")

  list(survivor_list = survivor_dfs, initial_screen = initial_screen)
}


# ── cfRes cross-fitted residualization (R backend, mirrors Python) ─────────────
.cross_fit_residuals_R <- function(X_corr, y_num, k_folds = 5L, ntree = 300L,
                                    mtry_factor = 1, nodesize = 5L,
                                    seed = 42L, num_processors = 1L) {
  n <- length(y_num); p_corr <- ncol(X_corr)
  set.seed(seed)
  perm     <- sample.int(n)
  fold_ids <- integer(n)
  for (k in seq_len(k_folds)) fold_ids[perm[seq(k, n, by = k_folds)]] <- k
  y_hat <- numeric(n)
  mtry  <- max(1L, min(p_corr, ceiling(mtry_factor * sqrt(p_corr))))
  for (k in seq_len(k_folds)) {
    tr <- which(fold_ids != k); te <- which(fold_ids == k)
    rf <- ranger::ranger(x = X_corr[tr, , drop = FALSE], y = y_num[tr],
                         num.trees = as.integer(ntree), mtry = mtry,
                         min.node.size = as.integer(nodesize),
                         probability = FALSE, verbose = FALSE,
                         num.threads = as.integer(num_processors),
                         seed = as.integer(seed + k * 1999L))
    y_hat[te] <- stats::predict(rf, data = X_corr[te, , drop = FALSE])$predictions
  }
  y_num - y_hat
}


# ── Phase 2: Shadow-stability selection (R backend) ──────────────────────────

.mf_R_select <- function(X_surv, y, survivor_list, select_params,
                          module_membership_surv,
                          nodesize, num_processors, seed, verbose,
                          r_shap = "permutation") {
  sl             <- select_params
  CLASSIFICATION <- is.factor(y)
  use_ts         <- identical(r_shap, "treeshap")
  n              <- nrow(X_surv)
  all_feats      <- colnames(X_surv)
  shadow_pct     <- sl$shadow_percentile / 100.0

  y_num     <- if (CLASSIFICATION) as.numeric(as.integer(y)) else as.numeric(y)
  unassigned_set <- if (!is.null(sl$unassigned_modules)) as.character(sl$unassigned_modules) else character(0L)

  # cfRes residualized target for the unassigned pool (scope = "unassigned")
  y_for_ind <- y_num
  if (isTRUE(sl$use_cfres) && identical(sl$cfres_scope, "unassigned") &&
      length(unassigned_set) > 0L) {
    corr_feats <- all_feats[!(module_membership_surv[all_feats] %in% unassigned_set)]
    if (length(corr_feats) > 0L) {
      if (verbose) cat(sprintf("  cfRes: residualizing y by %d corr features (%d folds)\n",
                               length(corr_feats), sl$cfres_n_folds))
      y_for_ind <- .cross_fit_residuals_R(
        X_surv[, corr_feats, drop = FALSE], y_num,
        k_folds = sl$cfres_n_folds, ntree = sl$cfres_ntree,
        mtry_factor = sl$mtry_factor, nodesize = nodesize,
        seed = seed + 88317L, num_processors = num_processors)
    }
  }

  # ── inner: one pool's shadow-bootstrap loop ────────────────────────────────
  run_pool <- function(pool_feats, y_pool, classif_pool, pool_seed) {
    np <- length(pool_feats)
    freq_map <- setNames(integer(np), pool_feats)
    prev_freqs <- NULL; n_done <- 0L
    for (b in seq_len(sl$n_boots)) {
      bseed <- pool_seed + b * 97L
      set.seed(bseed)
      idx <- sample.int(n, max(1L, floor(n / 2L)), replace = FALSE)
      Xb  <- X_surv[idx, pool_feats, drop = FALSE]
      yb  <- y_pool[idx]
      Xb_aug   <- Xb
      shad_nms <- paste0("__shad__", pool_feats)
      for (fi in seq_along(pool_feats))
        Xb_aug[[shad_nms[fi]]] <- sample(Xb[[pool_feats[fi]]])
      cp_aug  <- ncol(Xb_aug)
      n_trees <- max(as.integer(sl$min_ntree), as.integer(sl$ntree_factor * cp_aug))
      mtry_r  <- .compute_mtry_aug_R(np, sl$mtry_factor, classif_pool,
                                     sl$mtry_rule, sl$mtry_on_real)
      rf_b <- ranger::ranger(
        x = Xb_aug, y = yb, mtry = mtry_r, num.trees = n_trees,
        importance = "permutation", min.node.size = as.integer(nodesize),
        probability = classif_pool, verbose = FALSE,
        num.threads = as.integer(num_processors), seed = as.integer(bseed))
      vim_b    <- if (use_ts)
        .vim_treeshap_ranger(rf_b, Xb_aug, colnames(Xb_aug))
      else rf_b$variable.importance
      real_vim <- vim_b[pool_feats]; shad_vim <- vim_b[shad_nms]
      shad_thr <- stats::quantile(shad_vim, probs = shadow_pct, na.rm = TRUE)
      winners  <- pool_feats[!is.na(real_vim) & real_vim > shad_thr]
      freq_map[winners] <- freq_map[winners] + 1L
      n_done <- n_done + 1L
      if (isTRUE(sl$early_stop_boots) &&
          (b %% as.integer(sl$early_stop_check_every)) == 0L) {
        cur <- freq_map / n_done
        if (!is.null(prev_freqs) &&
            max(abs(cur - prev_freqs)) < as.numeric(sl$early_stop_tol)) break
        prev_freqs <- cur
      }
    }
    fr  <- freq_map / max(1L, n_done)
    ord <- order(fr, decreasing = TRUE)
    list(names = pool_feats[ord], freqs = as.numeric(fr[ord]))
  }

  # ── build pool specs (name, feats, is_unassigned) ───────────────────────────────
  shadow_mode <- sl$shadow_mode
  if (!shadow_mode %in% c("within_module", "split")) shadow_mode <- "split"
  pool_specs <- list()
  if (shadow_mode == "within_module") {
    grp <- split(all_feats, module_membership_surv[all_feats])
    for (nm in names(grp))
      pool_specs[[length(pool_specs) + 1L]] <-
        list(name = nm, feats = grp[[nm]], is_unassigned = nm %in% unassigned_set)
  } else if (length(unassigned_set) > 0L) {
    corr_feats  <- all_feats[!(module_membership_surv[all_feats] %in% unassigned_set)]
    unassigned_feats <- all_feats[  module_membership_surv[all_feats] %in% unassigned_set]
    if (length(corr_feats)  > 0L)
      pool_specs[[length(pool_specs)+1L]] <- list(name="corr",  feats=corr_feats,  is_unassigned=FALSE)
    if (length(unassigned_feats) > 0L)
      pool_specs[[length(pool_specs)+1L]] <- list(name="unassigned", feats=unassigned_feats, is_unassigned=TRUE)
  } else {
    pool_specs[[1L]] <- list(name = "all", feats = all_feats, is_unassigned = FALSE)
  }

  pool_stability <- list(); all_stable <- character(0L)
  for (si in seq_along(pool_specs)) {
    ps <- pool_specs[[si]]
    if (length(ps$feats) == 0L) next

    # cfRes-residualized target -> CONTINUOUS -> fit this pool as REGRESSION
    y_pool <- y; classif_pool <- CLASSIFICATION; use_resid <- FALSE
    if (isTRUE(sl$use_cfres)) {
      if (identical(sl$cfres_scope, "all_modules")) {
        other <- setdiff(all_feats, ps$feats)
        if (length(other) > 0L) {
          y_pool <- .cross_fit_residuals_R(
            X_surv[, other, drop = FALSE], y_num,
            k_folds = sl$cfres_n_folds, ntree = sl$cfres_ntree,
            mtry_factor = sl$mtry_factor, nodesize = nodesize,
            seed = seed + 88317L + si * 997L, num_processors = num_processors)
          classif_pool <- FALSE; use_resid <- TRUE
        }
      } else if (ps$is_unassigned) {           # scope = "unassigned"
        y_pool <- y_for_ind; classif_pool <- FALSE; use_resid <- TRUE
      }
    }

    if (verbose)
      cat(sprintf("  Shadow stability - pool '%s' (%d feats, %d boots)%s ...\n",
                  ps$name, length(ps$feats), sl$n_boots,
                  if (use_resid) " [y=y_resid -> regression]" else ""))

    rp <- run_pool(ps$feats, y_pool, classif_pool, seed + si)

    use_elbow_pool <- identical(sl$threshold_mode, "elbow") || ps$is_unassigned
    thr_val <- if (ps$is_unassigned && !is.null(sl$pi_thr_unassigned))
                 as.numeric(sl$pi_thr_unassigned) else as.numeric(sl$pi_thr)
    if (ps$is_unassigned && is.null(sl$pi_thr_unassigned)) use_elbow_pool <- TRUE
    res <- .apply_stable_threshold_R(rp$names, rp$freqs, thr_val,
                                      use_elbow = use_elbow_pool)
    all_stable <- c(all_stable, res$stable)
    pool_stability[[length(pool_stability) + 1L]] <- list(
      pool = ps$name, feature_names = rp$names, freqs = rp$freqs,
      threshold = res$threshold, selected = rp$names %in% res$stable,
      elbow_rank = res$elbow_rank)
  }

  stable_features    <- unique(all_stable)
  stability_freq_map <- setNames(numeric(length(stable_features)), stable_features)
  for (ps in pool_stability) {
    idx_s <- ps$feature_names %in% stable_features
    stability_freq_map[ps$feature_names[idx_s]] <- ps$freqs[idx_s]
  }
  selection_list <- list(data.frame(
    feature_name = all_feats, variable_importance = NA_real_,
    stringsAsFactors = FALSE))

  list(stable_features = stable_features, stability_freq = stability_freq_map,
       selection_list = selection_list,
       pool_stability = Filter(Negate(is.null), pool_stability))
}
