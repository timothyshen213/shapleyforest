# Pure-R backend for shapleyforest — no Python / reticulate required.
#
# Screening RFE and shadow-stability selection both use ranger + permutation
# importance.  For the final importance step the user chooses one of three
# methods via the `r_shap` argument passed down from sf():
#
#   "permutation" — ranger permutation VIM (scalar per feature, no extra pkg)
#   "fastshap"    — approximate Shapley values via the fastshap package
#                   (interaction_matrix = NULL)
#   "treeshap"    — exact TreeSHAP via the treeshap package; also computes an
#                   interaction matrix identical in structure to the Python
#                   backend's (diagonal = main effects, off-diagonal = pairwise)
#
# Internal helpers — not exported.


# ── Top-level orchestrator called by sf() ────────────────────────────────────

.sf_R_backend <- function(X, y,
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

  screen_out <- .sf_R_screen(
    X                 = X,
    y                 = y,
    module_membership = module_membership_screen,
    screen_params     = sc,
    nodesize          = nodesize,
    num_processors    = num_processors,
    seed              = seed,
    verbose           = verbose >= 2L
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

  sel_out <- .sf_R_select(
    X_surv                 = X_surv,
    y                      = y,
    survivor_list          = survivor_list,
    select_params          = sl,
    module_membership_surv = module_assigns_surv,
    nodesize               = nodesize,
    num_processors         = num_processors,
    seed                   = seed,
    verbose                = verbose >= 1L
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
      permutation = .sf_R_shap_permutation(final_rf, final_X, stable_features,
                                            verbose >= 1L),
      fastshap    = .sf_R_shap_fastshap(final_rf, final_X, y, stable_features,
                                         seed, verbose >= 1L),
      treeshap    = .sf_R_shap_treeshap(final_rf, final_X, stable_features,
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

  shapley_forest(
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

.sf_R_shap_permutation <- function(final_rf, final_X, stable_features, verbose) {
  if (verbose) cat("  Using permutation VIM (no per-observation SHAP values).\n")
  list(
    shap_matrix        = matrix(NA_real_, nrow = nrow(final_X),
                                ncol     = length(stable_features),
                                dimnames = list(NULL, stable_features)),
    interaction_matrix = NULL
  )
}


# ── SHAP method 2: fastshap approximate SHAP ─────────────────────────────────

.sf_R_shap_fastshap <- function(final_rf, final_X, y, stable_features,
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

.sf_R_shap_treeshap <- function(final_rf, final_X, stable_features, verbose) {
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

.sf_R_screen <- function(X, y, module_membership, screen_params,
                          nodesize, num_processors, seed, verbose) {
  sc             <- screen_params
  CLASSIFICATION <- is.factor(y)
  module_list    <- unique(module_membership)

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

      vims   <- rf_step$variable.importance
      n_drop <- max(1L, floor(sc$drop_fraction * cp))
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

    vims_fin <- rf_fin$variable.importance

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


# ── Phase 2: Shadow-stability selection (R backend) ──────────────────────────

.sf_R_select <- function(X_surv, y, survivor_list, select_params,
                          module_membership_surv,
                          nodesize, num_processors, seed, verbose) {
  sl             <- select_params
  CLASSIFICATION <- is.factor(y)
  n              <- nrow(X_surv)
  all_feats      <- colnames(X_surv)
  shadow_pct     <- sl$shadow_percentile / 100.0   # 95 -> 0.95

  # Pool grouping
  shadow_mode <- sl$shadow_mode
  if (!shadow_mode %in% c("within_module", "split")) shadow_mode <- "split"

  if (shadow_mode == "within_module") {
    pools          <- split(all_feats, module_membership_surv[all_feats])
    use_elbow_pool <- TRUE
  } else {
    pools          <- list(all = all_feats)
    use_elbow_pool <- FALSE
  }

  pool_stability <- vector("list", length(pools))
  pool_nms       <- names(pools)
  all_stable     <- character(0L)

  for (pi in seq_along(pool_nms)) {
    pool_nm    <- pool_nms[pi]
    pool_feats <- pools[[pool_nm]]
    np         <- length(pool_feats)
    if (np == 0L) next

    if (verbose)
      cat(sprintf("  Shadow stability — pool '%s' (%d features, %d boots) ...\n",
                  pool_nm, np, sl$n_boots))

    freq_map   <- setNames(integer(np), pool_feats)
    prev_freqs <- NULL
    n_done     <- 0L

    for (b in seq_len(sl$n_boots)) {
      bseed <- seed + b * 97L
      set.seed(bseed)

      idx  <- sample.int(n, max(1L, floor(n / 2L)), replace = FALSE)
      Xb   <- X_surv[idx, pool_feats, drop = FALSE]
      yb   <- y[idx]

      # Augment with permuted shadow columns
      Xb_aug   <- Xb
      shad_nms <- paste0("__shad__", pool_feats)
      for (fi in seq_along(pool_feats))
        Xb_aug[[shad_nms[fi]]] <- sample(Xb[[pool_feats[fi]]])

      cp_aug  <- ncol(Xb_aug)
      n_trees <- max(as.integer(sl$min_ntree),
                     as.integer(sl$ntree_factor * cp_aug))
      mtry_r  <- .compute_mtry_aug_R(np, sl$mtry_factor, CLASSIFICATION,
                                     sl$mtry_rule, sl$mtry_on_real)

      rf_b <- ranger::ranger(
        x             = Xb_aug,
        y             = yb,
        mtry          = mtry_r,
        num.trees     = n_trees,
        importance    = "permutation",
        min.node.size = as.integer(nodesize),
        probability   = CLASSIFICATION,
        verbose       = FALSE,
        num.threads   = as.integer(num_processors),
        seed          = as.integer(bseed)
      )

      vim_b    <- rf_b$variable.importance
      real_vim <- vim_b[pool_feats]
      shad_vim <- vim_b[shad_nms]

      shad_thr <- stats::quantile(shad_vim, probs = shadow_pct, na.rm = TRUE)
      winners  <- pool_feats[!is.na(real_vim) & real_vim > shad_thr]
      freq_map[winners] <- freq_map[winners] + 1L
      n_done <- n_done + 1L

      # adaptive early stopping: halt once frequencies have converged
      if (isTRUE(sl$early_stop_boots) &&
          (b %% as.integer(sl$early_stop_check_every)) == 0L) {
        cur_freqs <- freq_map / n_done
        if (!is.null(prev_freqs) &&
            max(abs(cur_freqs - prev_freqs)) < as.numeric(sl$early_stop_tol)) break
        prev_freqs <- cur_freqs
      }
    }

    freqs_raw <- freq_map / max(1L, n_done)
    ord       <- order(freqs_raw, decreasing = TRUE)
    names_s   <- pool_feats[ord]
    freqs_s   <- freqs_raw[ord]

    if (identical(sl$threshold_mode, "elbow")) use_elbow_pool <- TRUE
    res <- .apply_stable_threshold_R(names_s, freqs_s,
                                      sl$pi_thr,
                                      use_elbow = use_elbow_pool)
    all_stable <- c(all_stable, res$stable)

    pool_stability[[pi]] <- list(
      pool          = pool_nm,
      feature_names = names_s,
      freqs         = as.numeric(freqs_s),
      threshold     = res$threshold,
      selected      = names_s %in% res$stable,
      elbow_rank    = res$elbow_rank
    )
  }

  stable_features    <- unique(all_stable)
  stability_freq_map <- setNames(numeric(length(stable_features)), stable_features)
  for (ps in pool_stability) {
    if (is.null(ps)) next
    idx_s <- ps$feature_names %in% stable_features
    stability_freq_map[ps$feature_names[idx_s]] <- ps$freqs[idx_s]
  }

  selection_list <- list(data.frame(
    feature_name        = all_feats,
    variable_importance = NA_real_,
    stringsAsFactors    = FALSE
  ))

  list(
    stable_features = stable_features,
    stability_freq  = stability_freq_map,
    selection_list  = selection_list,
    pool_stability  = Filter(Negate(is.null), pool_stability)
  )
}
