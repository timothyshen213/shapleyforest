"""
sf_python_backend_v18.py
------------------------
ShapleyForest v18 — stable union + single-shot joint TreeSHAP.

Replaces pass-2 SHAP RFE (iterative elimination) with a single RF fit on
the entire stable union, then TreeSHAP on all features simultaneously.

Key change vs v14/v17
---------------------
  v14/v17: stable_union -> iterative RFE (many RF refits) -> exactly K features
  v18:     stable_union -> ONE RF + TreeSHAP -> return ALL stable features
                           with mean|SHAP| and full signed SHAP matrix

Why single-shot is better for independent features
---------------------------------------------------
  In iterative RFE, removing a corr feature donates its predictive credit to
  surviving corr features (their SHAPs inflate).  The last corr feature
  standing collects the entire corr module credit Sigma^2.  An ind feature
  needs sigma^2 > Sigma^2 to rank above it -- nearly impossible when the
  large corr module dominates.

  In a joint model, corr features SHARE credit: each gets ~Sigma^2/n_corr
  in SHAP.  The ind feature keeps its full sigma^2.  It wins high rank when
  sigma^2 > Sigma^2/n_corr -- a far weaker condition, satisfied when the
  stable union is of reasonable size.  No path dependency, no K trim, no
  bypass heuristic needed.

Output
------
  Returns ALL stable features (|U| varies) with:
    final_features  list[str]         stable features sorted by mean|SHAP|
    final_vims      list[float]       mean|SHAP| per feature
    stability_freq  list[float]       shadow boot frequency per feature
    shap_matrix     list[list[float]] signed n x |U| SHAP matrix for CIs

  Empty stable union -> empty lists (not padded to any K).

DML layer (optional)
---------------------
  use_dml_residual=True  cross-fitted corr residualization (as in v17).
  Applied BEFORE shadow stability to boost ind feature frequencies into U.
  Final single-shot RF always uses original y for interpretable SHAPs.

Shadow modes
------------
  "split"         one corr pool, one ind pool, separate shadow runs
  "within_module" one shadow run per module
  (dual_pool removed -- pass-2 architecture not needed here)

Shadow-stability algorithm (unchanged from v14/v17)
----------------------------------------------------
  n_boots half-sample bootstrap rounds:
    1. Subsample n//2 rows.
    2. Permute each real column -> X_shadow.
    3. X_aug = [X_real | X_shadow]  (n//2 x 2p).
    4. Fit RF; compute SHAP on X_aug.
    5. shadow_threshold = percentile(mean|SHAP| shadow cols, shadow_percentile).
    6. selected_b = { j : mean|SHAP|_real(j) > shadow_threshold }
  stability_freq(j) = |{ b : j in selected_b }| / n_boots

  Thresholding:
    corr pool: stable = { j : freq(j) >= pi_thr }
    ind  pool: elbow of sorted freq curve (or pi_thr_indep if given)

Screen RFE and run_final_rf_shap unchanged from v12/v14.
"""

import math
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers  (identical to v12)
# ─────────────────────────────────────────────────────────────────────────────

# ── v17: cross-fitted corr residualization ────────────────────────────────────

def _cross_fit_residuals(X_corr, y, k_folds=5, ntree=300,
                          mtry_factor=1.0, nodesize=5, seed=42, n_jobs=1):
    """
    Phase 2a — Robinson / DML residualization.

    Fits RF on corr-pool survivors via K-fold cross-fitting: each obs's
    predicted value ĝ(X_corr_i) comes from a model trained on the other
    K-1 folds — no obs influences its own residual.

    y_resid = y - ĝ_cv(X_corr)

    After this step, y_resid ≈ θᵀ·X_ind + ε (independent signal + noise),
    dramatically boosting the ind pool's effective SNR.

    Parameters
    ----------
    X_corr : ndarray (n, p_corr)   corr-pool feature matrix
    y      : ndarray (n,)           original continuous outcome
    k_folds: int                    cross-fitting folds (default 5)
    ntree  : int                    trees per fold-RF (default 300)

    Returns
    -------
    y_resid : ndarray (n,)   residuals; same dtype/shape as y
    """
    n       = len(y)
    p_corr  = X_corr.shape[1]
    rng     = np.random.default_rng(int(seed))

    # assign fold ids by cyclic permutation → balanced folds
    perm     = rng.permutation(n)
    fold_ids = np.empty(n, dtype=int)
    for k in range(k_folds):
        fold_ids[perm[k::k_folds]] = k

    y_hat   = np.zeros(n, dtype=np.float64)
    mtry    = max(1, min(p_corr,
                         int(math.ceil(mtry_factor * math.sqrt(p_corr)))))

    for k in range(k_folds):
        tr = np.where(fold_ids != k)[0]
        te = np.where(fold_ids == k)[0]

        rf = RandomForestRegressor(
            n_estimators     = int(ntree),
            max_features     = mtry,
            min_samples_leaf = int(nodesize),
            random_state     = (int(seed) + k * 1999) % (2**31 - 1),
            n_jobs           = int(n_jobs),
        )
        rf.fit(X_corr[tr, :], y[tr].astype(np.float64))
        y_hat[te] = rf.predict(X_corr[te, :])

    y_resid = y.astype(np.float64) - y_hat
    r2_corr = float(1.0 - np.var(y_resid) / (np.var(y) + 1e-12))
    # DML residualization diagnostics suppressed in release build
    return y_resid

def _build_rf(n_estimators, max_features, min_samples_leaf,
              random_state, n_jobs, classification):
    max_features = max(1, int(max_features))
    kw = dict(
        n_estimators     = int(n_estimators),
        max_features     = max_features,
        min_samples_leaf = int(min_samples_leaf),
        random_state     = int(random_state),
        n_jobs           = int(n_jobs),
    )
    return RandomForestClassifier(**kw) if classification else RandomForestRegressor(**kw)


def _shap_mean_abs(rf, X, classification):
    import shap
    explainer = shap.TreeExplainer(rf)
    shap_raw  = explainer.shap_values(X)

    if isinstance(shap_raw, list):
        if classification and len(rf.classes_) == 2:
            sv = np.abs(shap_raw[1])
        elif classification:
            sv = np.mean([np.abs(s) for s in shap_raw], axis=0)
        else:
            sv = np.abs(shap_raw[0])
    else:
        if shap_raw.ndim == 3:
            if classification and len(rf.classes_) == 2:
                sv = np.abs(shap_raw[:, :, 1])
            else:
                sv = np.mean(np.abs(shap_raw), axis=2)
        else:
            sv = np.abs(shap_raw)

    return np.mean(sv, axis=0)   # shape (2p,) when called with augmented X



def _shap_interaction_matrix(rf, X, classification):
    """
    Compute exact TreeSHAP interaction values and return the p×p mean-|interaction|
    summary matrix.

    Uses ``TreeExplainer.shap_interaction_values(X)`` (Lundberg et al. 2018).

    Diagonal entry [j, j] = mean_i |main_effect(j, obs_i)|  (same as mean |SHAP|).
    Off-diagonal   [j, k] = mean_i |interaction(j, k, obs_i)|  (symmetric).

    Parameters
    ----------
    rf             : fitted sklearn RF
    X              : ndarray (n, p)   feature matrix (stable features only)
    classification : bool

    Returns
    -------
    ndarray (p, p)   mean absolute interaction matrix (symmetric, non-negative)
    """
    import shap
    explainer       = shap.TreeExplainer(rf)
    interaction_raw = explainer.shap_interaction_values(X)

    # shap_interaction_values shapes vary by version / problem type:
    #   regression          -> ndarray (n, p, p)
    #   binary clf (list)   -> [class0 (n,p,p), class1 (n,p,p)]
    #   multiclass (list)   -> [c0 (n,p,p), c1 (n,p,p), ...]
    #   newer shap (ndarray)-> (n, p, p, n_classes)  [4-D]
    if isinstance(interaction_raw, list):
        if classification and len(rf.classes_) == 2:
            iv = np.array(interaction_raw[1])           # positive class (n, p, p)
        elif classification:
            iv = np.mean([np.abs(np.array(s))
                          for s in interaction_raw], axis=0)  # avg across classes
        else:
            iv = np.array(interaction_raw[0])           # regression (n, p, p)
    else:
        arr = np.array(interaction_raw)
        if arr.ndim == 4:                               # (n, p, p, n_classes)
            if classification and arr.shape[3] == 2:
                iv = arr[:, :, :, 1]
            elif classification:
                iv = np.mean(np.abs(arr), axis=3)
            else:
                iv = arr[:, :, :, 0]
        else:
            iv = arr                                    # (n, p, p)

    return np.mean(np.abs(iv), axis=0)                 # (p, p) mean |interaction|


def _shap_full_signed(rf, X, classification):
    import shap
    explainer = shap.TreeExplainer(rf)
    shap_raw  = explainer.shap_values(X)

    if isinstance(shap_raw, list):
        if classification and len(rf.classes_) == 2:
            return shap_raw[1]
        elif classification:
            return np.mean(shap_raw, axis=0)
        else:
            return shap_raw[0]
    else:
        if shap_raw.ndim == 3:
            if classification and len(rf.classes_) == 2:
                return shap_raw[:, :, 1]
            else:
                return np.mean(shap_raw, axis=2)
        else:
            return shap_raw


def _permutation_importance(rf):
    return rf.feature_importances_


def _compute_mtry(n_features, mtry_factor, classification, rule="auto"):
    """
    Features considered per split.
      rule="auto"     - p/3 (classification) / sqrt(p) (regression)  [historical]
      rule="sqrt"     - sqrt(p) for both
      rule="p_over_3" - p/3 for both
    """
    if rule == "sqrt":
        base = math.sqrt(n_features)
    elif rule == "p_over_3":
        base = n_features / 3.0
    else:
        base = (n_features / 3.0) if classification else math.sqrt(n_features)
    return min(max(1, math.ceil(mtry_factor * base)), n_features)


def _compute_mtry_augmented(p_real, mtry_factor, classification,
                            rule="auto", on_real=False):
    """
    mtry on the augmented [real | shadow] matrix (2*p_real cols).
      on_real=False - mtry = rule(2p)          [historical]
      on_real=True  - mtry = 2 * rule(p_real)  [keeps real-feature sampling
                      undiluted by the shadows; statistical, not a speedup]
    """
    n_aug = 2 * p_real
    if on_real:
        return min(max(1, 2 * _compute_mtry(p_real, mtry_factor, classification, rule)), n_aug)
    return _compute_mtry(n_aug, mtry_factor, classification, rule)


def _compute_ntree(n_features, ntree_factor, min_ntree):
    return max(int(n_features * ntree_factor), int(min_ntree))


def _map_labels(y):
    uniq = sorted(np.unique(y))
    lmap = {v: i for i, v in enumerate(uniq)}
    return np.array([lmap[v] for v in y], dtype=int)


def _find_elbow_idx(freqs_sorted):
    """
    Perpendicular-distance elbow detection on a descending frequency array.

    Normalises both axes to [0,1], draws the line from the first to the
    last point, and returns the 0-based index of the point that is
    furthest from that line (the elbow).

    The stable set should include all features from index 0 up to and
    including the elbow index (i.e. freq >= freqs_sorted[elbow_idx]).
    """
    n = len(freqs_sorted)
    if n <= 1:
        return 0
    freqs = np.array(freqs_sorted, dtype=np.float64)
    rng   = freqs[0] - freqs[-1]
    if rng < 1e-10:          # all flat -> take nothing meaningful
        return 0
    x  = np.linspace(0.0, 1.0, n)
    yr = (freqs - freqs[-1]) / rng      # normalise to [0,1]
    # Perpendicular distance from point i to the chord A=(x[0],yr[0])->B=(x[-1],yr[-1]).
    # Correct formula: |(yr[-1]-yr[0])*(x_i-x[0]) - (x[-1]-x[0])*(yr_i-yr[0])| / |AB|
    # This gives 0 at both endpoints (correct) and finds the farthest interior point.
    dx = x[-1]  - x[0]                 # = 1
    dy = yr[-1] - yr[0]                 # = -1
    dists = np.abs(dy * (x - x[0]) - dx * (yr - yr[0])) / np.sqrt(dx**2 + dy**2)
    return int(np.argmax(dists))


# ─────────────────────────────────────────────────────────────────────────────
# Core RFE loop  (identical to v12 — used in screen and pass-2 trim)
# ─────────────────────────────────────────────────────────────────────────────

def _rfe_loop(X_mod, feat_names, y, target,
              drop_fraction, mtry_factor, ntree_factor, min_ntree,
              nodesize, classification, seed, n_jobs, use_shap=True,
              shap_bg=None, mtry_rule="auto"):
    """
    shap_bg : ndarray (n_bg, p_full) or None
        If provided AND use_shap, switches the SHAP call to interventional
        TreeSHAP using X_bg (subsetted to the current surviving columns each
        iteration).  Columns of shap_bg must be aligned with feat_names on
        entry; we subset alongside current_X as features are dropped.
    """
    feat_names    = list(feat_names)
    current_X     = X_mod.copy()
    current_names = list(feat_names)
    current_bg    = shap_bg.copy() if shap_bg is not None else None
    is_initial    = True
    initial_info  = None
    survivors     = None
    survivor_vims = None
    all_steps     = []   # one entry per RFE iteration: (names_sorted, vims_sorted)

    while len(current_names) >= target:
        n_cur = len(current_names)
        mtry  = _compute_mtry(n_cur, mtry_factor, classification, mtry_rule)
        ntree = _compute_ntree(n_cur, ntree_factor, min_ntree)

        rf = _build_rf(ntree, mtry, int(nodesize), seed, n_jobs, classification)
        rf.fit(current_X, y.astype(int) if classification else y.astype(float))

        if use_shap:
            vim = _shap_mean_abs(rf, current_X, classification)
        else:
            vim = _permutation_importance(rf)

        order        = np.argsort(vim)[::-1]
        vim_sorted   = vim[order]
        names_sorted = [current_names[i] for i in order]
        reduction    = math.ceil(n_cur * drop_fraction)

        # record this iteration (all features ranked at this step)
        all_steps.append((list(names_sorted), vim_sorted.tolist()))

        if is_initial:
            n_keep_init   = n_cur - reduction if (n_cur - reduction) > target else target
            survived_init = ([True] * n_keep_init + [False] * (n_cur - n_keep_init))
            initial_info  = {
                'feature' : names_sorted,
                'vim'     : vim_sorted.tolist(),
                'survivor': survived_init,
            }
            is_initial = False

        if n_cur - reduction > target:
            keep_n        = n_cur - reduction
            keep_idx      = list(order[:keep_n])
            current_X     = current_X[:, keep_idx]
            current_names = names_sorted[:keep_n]
            if current_bg is not None:
                current_bg = current_bg[:, keep_idx]
        else:
            survivors     = names_sorted[:target]
            survivor_vims = vim_sorted[:target].tolist()
            break

    if survivors is None:
        survivors     = current_names[:target]
        survivor_vims = [0.0] * len(survivors)

    return survivors, survivor_vims, initial_info, all_steps


# ─────────────────────────────────────────────────────────────────────────────
# Blockwise pre-reduction  (identical to v12)
# ─────────────────────────────────────────────────────────────────────────────

def _blockwise_reduce(X_mod, feat_names, y, max_block_size, keep_fraction,
                      drop_fraction, mtry_factor, ntree_factor, min_ntree,
                      nodesize, classification, seed, n_jobs, use_shap=True,
                      mtry_rule="auto"):
    feat_names = list(feat_names)
    rng = np.random.default_rng(seed)

    while X_mod.shape[1] > max_block_size:
        n_feat   = X_mod.shape[1]
        perm     = rng.permutation(n_feat)
        n_blocks = math.ceil(n_feat / max_block_size)
        blocks   = np.array_split(perm, n_blocks)

        kept_cols = []
        for block in blocks:
            X_blk      = X_mod[:, block]
            names_blk  = [feat_names[i] for i in block]
            n_blk      = len(block)
            target_blk = max(1, math.ceil(n_blk * keep_fraction))
            surv_b, _, _, _ = _rfe_loop(
                X_blk, names_blk, y, target_blk,
                drop_fraction, mtry_factor, ntree_factor, min_ntree,
                nodesize, classification, seed, n_jobs, use_shap,
                mtry_rule=mtry_rule
            )
            name_to_orig = {feat_names[i]: i for i in block}
            kept_cols.extend([name_to_orig[f] for f in surv_b])

        kept_cols  = sorted(set(kept_cols))
        X_mod      = X_mod[:, kept_cols]
        feat_names = [feat_names[i] for i in kept_cols]

    return X_mod, feat_names


# ─────────────────────────────────────────────────────────────────────────────
# Pass-2 RFE trim  (identical to v12 — plain SHAP ranking)
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# v14 CORE: shadow-feature stability selection on one feature pool
# ─────────────────────────────────────────────────────────────────────────────

def _boruta_shadow_stability(X_pool, feat_names, y,
                              n_boots,
                              mtry_factor, ntree_factor, min_ntree,
                              nodesize, classification, seed, n_jobs,
                              use_shap=True, shadow_percentile=95,
                              mtry_rule="auto", mtry_on_real=False,
                              early_stop_boots=False, early_stop_tol=0.01,
                              early_stop_check_every=10):
    """
    Shadow-feature stability selection on a single feature pool.

    Per bootstrap b:
      1. Subsample n//2 rows of X_pool and y.
      2. Column-permute each real feature independently -> X_shadow
         (same shape as X_pool_sub; destroys any real association with y).
      3. X_aug = [X_real_sub | X_shadow_sub]  (n//2 × 2p).
      4. Fit RF on X_aug; compute mean|SHAP| on X_aug (2p values).
      5. shadow_threshold = np.percentile(mean|SHAP|[p:], shadow_percentile)
         (only the p shadow columns contribute to the null).
      6. selected_b = { j ∈ {0..p-1} : mean|SHAP|[j] > shadow_threshold }

    stability_freq(j) = |selected_b| / n_boots

    NOTE: thresholding (pi_thr or elbow) is NOT applied here — callers
    apply the appropriate cutoff for their pool type (corr vs ind).

    Returns
    -------
    names_sorted : list[str]   all pool features, sorted by freq desc
    freqs_sorted : list[float] parallel stability frequencies
    """
    feat_names = list(feat_names)
    p = X_pool.shape[1]
    n = X_pool.shape[0]

    if p == 0:
        return [], []

    counts  = np.zeros(p, dtype=np.float64)
    n_valid = 0

    # Per-bootstrap deterministic seeds (order-independent). Mirrors the
    # shapleyforest_py package so the two implementations stay matched.
    boot_seeds = np.random.SeedSequence(int(seed)).generate_state(int(n_boots))
    _prev_freqs = [None]   # early-stopping state

    for b in range(int(n_boots)):
        boot_seed = int(boot_seeds[b])
        rng = np.random.default_rng(boot_seed)

        idx = rng.choice(n, size=n // 2, replace=False)
        X_real_sub = X_pool[idx, :]
        y_sub      = y[idx]

        if classification and len(np.unique(y_sub)) < 2:
            continue

        # ── build shadow matrix: permute each column independently (C-level) ──
        n_sub = X_real_sub.shape[0]
        X_aug = np.empty((n_sub, 2 * p), dtype=X_pool.dtype)
        X_aug[:, :p] = X_real_sub
        X_aug[:, p:] = rng.permuted(X_real_sub, axis=0)
        n_aug = 2 * p

        mtry_aug = _compute_mtry_augmented(p, mtry_factor, classification,
                                           mtry_rule, mtry_on_real)
        ntree    = _compute_ntree(n_aug, ntree_factor, min_ntree)

        rf = _build_rf(ntree, mtry_aug, int(nodesize), boot_seed, n_jobs, classification)
        try:
            rf.fit(X_aug, y_sub.astype(int) if classification else y_sub.astype(float))
        except Exception:
            continue

        # ── importance on augmented X (length 2p) ─────────────────────────────
        if use_shap:
            try:
                vim_aug = _shap_mean_abs(rf, X_aug, classification)
            except Exception:
                vim_aug = _permutation_importance(rf)
        else:
            vim_aug = _permutation_importance(rf)

        vim_real   = vim_aug[:p]    # first p cols = real features
        vim_shadow = vim_aug[p:]    # last  p cols = shadow features

        # ── shadow threshold ──────────────────────────────────────────────────
        threshold = np.percentile(vim_shadow, float(shadow_percentile))

        # ── count features that beat their own shadow null ────────────────────
        selected_b = np.where(vim_real > threshold)[0]
        counts[selected_b] += 1.0
        n_valid += 1

        # adaptive early stopping: halt once frequencies have converged
        if early_stop_boots and n_valid > 0 and (b + 1) % int(early_stop_check_every) == 0:
            cur_freqs = counts / n_valid
            if _prev_freqs[0] is not None and \
               np.max(np.abs(cur_freqs - _prev_freqs[0])) < float(early_stop_tol):
                break
            _prev_freqs[0] = cur_freqs

    if n_valid == 0:
        return feat_names, [0.0] * p

    freqs = counts / n_valid
    order = np.argsort(freqs)[::-1]

    names_sorted = [feat_names[i] for i in order]
    freqs_sorted = freqs[order].tolist()

    return names_sorted, freqs_sorted



def _apply_stable_threshold(names_sorted, freqs_sorted, pi_thr,
                              use_elbow=False, label="pool"):
    """
    Apply stable-set threshold to shadow-stability outputs.

    use_elbow=True  (independent pools):
        Use perpendicular-distance elbow detection on the frequency curve.
        Keeps all features with freq >= freqs_sorted[elbow_idx].
        Falls back to top-1 if the curve is flat.

    use_elbow=False (correlated pools):
        Standard Meinshausen-Bühlmann threshold: freq >= pi_thr.

    Returns
    -------
    stable      : list[str]   selected feature names
    threshold   : float       frequency cutoff applied
    elbow_rank  : int or None 1-based rank of the elbow point (use_elbow only)
    """
    if not names_sorted:
        return [], 0.0, None

    freqs = np.array(freqs_sorted, dtype=np.float64)

    if use_elbow:
        elbow_idx  = _find_elbow_idx(freqs_sorted)
        elbow_freq = float(freqs[elbow_idx])
        stable = [f for f, fr in zip(names_sorted, freqs_sorted)
                  if fr >= elbow_freq and fr > 0.0]
        return stable, elbow_freq, int(elbow_idx) + 1   # 1-based rank
    else:
        thr    = float(pi_thr)
        stable = [f for f, fr in zip(names_sorted, freqs_sorted)
                  if fr >= thr]
        return stable, thr, None


# ─────────────────────────────────────────────────────────────────────────────
# Public API — Phase 1: Screen RFE  (identical to v12)
# ─────────────────────────────────────────────────────────────────────────────

def run_screen_rfe(X_mat, y, feature_names, module_assignments, module_list,
                   drop_fraction, keep_fraction,
                   mtry_factor, ntree_factor, min_ntree,
                   nodesize, classification, seed, n_jobs, max_modsize,
                   shap_model="full", mtry_rule="auto", verbose=True):
    X_mat              = np.array(X_mat, dtype=np.float64)
    y                  = np.array(y).ravel()
    feature_names      = [str(f) for f in feature_names]
    module_assignments = [str(m) for m in module_assignments]
    module_list        = [str(m) for m in module_list]
    classification     = bool(classification)
    drop_fraction      = float(drop_fraction)
    keep_fraction      = float(keep_fraction)
    mtry_factor        = float(mtry_factor)
    ntree_factor       = float(ntree_factor)
    min_ntree          = int(min_ntree)
    nodesize           = int(nodesize)
    seed               = int(seed)
    n_jobs             = int(n_jobs)
    max_modsize        = int(max_modsize)
    use_shap           = (shap_model == "full")

    if classification:
        y = _map_labels(y)
    else:
        y = y.astype(np.float64)

    feat_to_col = {f: i for i, f in enumerate(feature_names)}
    mod_to_cols = {}
    for feat, mod in zip(feature_names, module_assignments):
        mod_to_cols.setdefault(mod, []).append(feat_to_col[feat])

    all_survivor_features = []
    all_survivor_vims     = []
    init_feats    = []
    init_vims     = []
    init_survivors = []
    init_modules  = []

    for mod in module_list:
        cols = mod_to_cols.get(mod, [])
        if not cols:
            all_survivor_features.append([])
            all_survivor_vims.append([])
            continue

        X_mod    = X_mat[:, cols]
        feat_mod = [feature_names[c] for c in cols]
        n_mod    = len(cols)
        target   = max(1, math.ceil(n_mod * keep_fraction))

        if n_mod > max_modsize:
            X_mod, feat_mod = _blockwise_reduce(
                X_mod, feat_mod, y, max_modsize, 0.9,
                drop_fraction, mtry_factor, ntree_factor, min_ntree,
                nodesize, classification, seed, n_jobs, use_shap,
                mtry_rule=mtry_rule
            )
            n_mod  = X_mod.shape[1]
            target = max(1, math.ceil(n_mod * keep_fraction))

        survivors, vims, init_info, _ = _rfe_loop(
            X_mod, feat_mod, y, target,
            drop_fraction, mtry_factor, ntree_factor, min_ntree,
            nodesize, classification, seed, n_jobs, use_shap,
            mtry_rule=mtry_rule
        )

        all_survivor_features.append(survivors)
        all_survivor_vims.append(vims)

        if init_info is not None:
            n_rows = len(init_info['feature'])
            init_feats.extend(init_info['feature'])
            init_vims.extend(init_info['vim'])
            init_survivors.extend(init_info['survivor'])
            init_modules.extend([mod] * n_rows)

    return {
        'survivor_features' : all_survivor_features,
        'survivor_vims'     : all_survivor_vims,
        'initial_feat'      : init_feats,
        'initial_vim'       : init_vims,
        'initial_survivor'  : init_survivors,
        'initial_module'    : init_modules,
        'module_list'       : module_list,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Public API — Phase 2: v14 selection
# ─────────────────────────────────────────────────────────────────────────────

def run_select_rfe(X_surv_mat, y, feature_names_surv,
                   drop_fraction, number_selected,   # number_selected: ignored (no K trim)
                   mtry_factor, ntree_factor, min_ntree,
                   nodesize, classification, seed, n_jobs, max_modsize,
                   shap_model="full",
                   n_boots=50, pi_thr=0.60,
                   compute_interactions=True,
                   # ── shadow / threshold ────────────────────────────────────
                   shadow_mode="split",        # "split" or "within_module"
                   shadow_percentile=95,
                   pi_thr_indep=None,          # None = elbow; float = fixed
                   threshold_mode="pi_thr",    # "pi_thr" or "elbow"
                   # ── mtry ──────────────────────────────────────────────────
                   mtry_rule="auto", mtry_on_real=False,
                   # ── adaptive bootstrap stopping ───────────────────────────
                   early_stop_boots=False, early_stop_tol=0.01,
                   early_stop_check_every=10,
                   # ── DML ───────────────────────────────────────────────────
                   use_dml_residual=False,
                   dml_n_folds=5,
                   dml_ntree=300,
                   dml_scope="indep",
                   # ── group info ────────────────────────────────────────────
                   module_assignments_surv=None,
                   indep_modules=None,
                   verbose=True):
    """
    v18 selection: shadow stability -> stable union -> single-shot joint TreeSHAP.

    No pass-2 RFE.  No K trim.  Returns all stable features.

    shadow_mode = "split"
        Corr pool -> shadow stability -> stable_corr  (pi_thr threshold)
        Ind  pool -> shadow stability -> stable_ind   (elbow or pi_thr_indep)
        stable_union = stable_corr | stable_ind
        -> single RF on X[stable_union] + TreeSHAP -> return all |stable_union| features

    shadow_mode = "within_module"
        Per-module shadow stability -> per-module stable sets
        stable_union = union of all module stable sets
        -> single RF + TreeSHAP -> return all features

    use_dml_residual=True (optional):
        Before ind-pool shadow stability, compute y_resid = y - RF_cv(X_corr).
        Ind shadow stability runs on y_resid (boosted SNR for ind features).
        Final single-shot RF uses original y.

    Returns dict with:
        final_features  : stable features sorted by mean|SHAP| descending
        final_vims      : mean|SHAP| (importance) per feature
        stability_freq  : shadow boot frequency per feature
        shap_matrix     : signed n x |stable| SHAP values (rows=obs, cols=features)
        selection_list_feats / selection_list_vims : compat with R output structure
    """
    # ── coerce types ─────────────────────────────────────────────────────────
    X_surv_mat          = np.array(X_surv_mat, dtype=np.float64)
    y                   = np.array(y).ravel()
    orig_feature_names  = [str(f) for f in feature_names_surv]
    feature_names_surv  = orig_feature_names[:]
    drop_fraction       = float(drop_fraction)
    # number_selected intentionally not coerced -- not used in v18
    mtry_factor         = float(mtry_factor)
    ntree_factor        = float(ntree_factor)
    min_ntree           = int(min_ntree)
    nodesize            = int(nodesize)
    classification      = bool(classification)
    seed                = int(seed)
    n_jobs              = int(n_jobs)
    max_modsize         = int(max_modsize)
    use_shap            = (shap_model in ("full", "partial"))
    compute_interactions = bool(compute_interactions)
    n_boots             = int(n_boots)
    pi_thr              = float(pi_thr)
    shadow_mode         = str(shadow_mode).lower()
    shadow_percentile   = float(shadow_percentile)
    threshold_mode      = str(threshold_mode).lower()
    if threshold_mode not in ("pi_thr", "elbow"):
        raise ValueError(f"threshold_mode must be 'pi_thr' or 'elbow', got '{threshold_mode}'")
    use_dml_residual    = bool(use_dml_residual)
    dml_n_folds         = max(2, int(dml_n_folds))
    dml_ntree           = max(10, int(dml_ntree))
    dml_scope           = str(dml_scope).lower()
    if dml_scope not in ("indep", "all_modules"):
        raise ValueError(f"dml_scope must be 'indep' or 'all_modules', got '{dml_scope}'")
    # pi_thr_indep: None -> elbow; otherwise coerce to float
    if pi_thr_indep is not None:
        try:
            pi_thr_indep = float(pi_thr_indep)
        except (TypeError, ValueError):
            pi_thr_indep = None

    indep_set = set(str(m) for m in indep_modules) if indep_modules else set()
    if module_assignments_surv is not None:
        raw_mods    = [str(m) for m in module_assignments_surv]
        name_to_mod = dict(zip(orig_feature_names, raw_mods))
    else:
        name_to_mod = {}

    indep_strategy = ("elbow" if pi_thr_indep is None
                      else f"fixed pi={pi_thr_indep:.2f}")
    if verbose: print(f"shadow_mode={shadow_mode}  shadow_percentile={shadow_percentile}"
          f"  pi_thr={pi_thr}  pi_thr_indep={indep_strategy}  n_boots={n_boots}")
    if verbose: print(f"indep_modules={indep_modules}  "
          f"use_dml_residual={use_dml_residual}"
          + (f"  dml_scope={dml_scope}  dml_n_folds={dml_n_folds}  dml_ntree={dml_ntree}"
             if use_dml_residual else "")
          )

    if classification:
        y = _map_labels(y)
    else:
        y = y.astype(np.float64)

    # (Phase-2 blockwise pre-reduction removed by design - the full survivor
    #  pool goes straight into shadow stability.)

    # ── Phase 2a: DML corr residualization  (v17 new) ────────────────────────
    # Compute y_resid ONCE, before the mode branches, using all corr survivors.
    # This requires knowing corr_idx — find it now regardless of shadow_mode.
    # (The per-mode branches will re-derive it; that's fine, it's cheap.)
    if use_dml_residual and indep_set and name_to_mod:
        all_corr_idx = [i for i, f in enumerate(feature_names_surv)
                        if name_to_mod.get(f, "") not in indep_set]
        if len(all_corr_idx) == 0:
            if verbose: print("DML: no corr survivors - skipping residualization")
            y_for_ind = y
        else:
            X_corr_surv = X_surv_mat[:, all_corr_idx]
            if verbose: print(f"DML: computing cross-fitted residuals "
                  f"({dml_n_folds} folds, ntree={dml_ntree}, "
                  f"p_corr={len(all_corr_idx)})")
            y_for_ind = _cross_fit_residuals(
                X_corr   = X_corr_surv,
                y        = y,
                k_folds  = dml_n_folds,
                ntree    = dml_ntree,
                mtry_factor = float(mtry_factor),
                nodesize    = int(nodesize),
                seed        = int(seed) + 88317,
                n_jobs      = int(n_jobs),
            )
    else:
        y_for_ind = y   # no DML: ind pool uses original y (v14 behaviour)

    # ── helper: run shadow stability + threshold for one column subset ───────
    def _shadow_on_cols(col_idx, label, seed_offset=0, is_indep=False, y_override=None):
        """
        Returns (names_sort, freq_map, stable, threshold, freqs_sort, elbow_rank)
          names_sort  : features sorted by freq desc
          freq_map    : {feature: freq}
          stable      : selected features (above threshold)
          threshold   : cutoff frequency applied
          freqs_sort  : parallel frequency list to names_sort
          elbow_rank  : 1-based elbow rank (ind+elbow only), else None
        """
        if len(col_idx) == 0:
            return [], {}, [], 0.0, [], None
        X_sub   = X_surv_mat[:, col_idx]
        names_s = [feature_names_surv[i] for i in col_idx]
        # Priority: explicit y_override > y_for_ind (if is_indep) > y
        if y_override is not None:
            y_arr = y_override
        elif is_indep:
            y_arr = y_for_ind
        else:
            y_arr = y
        _using_resid = (y_override is not None) or (is_indep and use_dml_residual)
        # A DML-residualized target is CONTINUOUS even for a classification task,
        # so this pool's shadow-stability RF must be fit as a REGRESSOR (else the
        # residual is cast to int -> degenerate).
        pool_classification = classification and not _using_resid
        if verbose: print(f"shadow stability on {label}: {len(col_idx)} features"
              + (" [y=y_resid]" if _using_resid else ""))
        names_sort, freqs_sort = _boruta_shadow_stability(
            X_pool           = X_sub,
            feat_names       = names_s,
            y                = y_arr,
            n_boots          = n_boots,
            mtry_factor      = mtry_factor,
            ntree_factor     = ntree_factor,
            min_ntree        = min_ntree,
            nodesize         = nodesize,
            classification   = pool_classification,
            seed             = seed + seed_offset,
            n_jobs           = n_jobs,
            use_shap         = use_shap,
            shadow_percentile= shadow_percentile,
            mtry_rule        = mtry_rule,
            mtry_on_real     = mtry_on_real,
            early_stop_boots = early_stop_boots,
            early_stop_tol   = early_stop_tol,
            early_stop_check_every = early_stop_check_every,
        )
        freq_map = dict(zip(names_sort, freqs_sort))
        # ── threshold: corr -> pi_thr; ind -> elbow (or pi_thr_indep) ─────────
        if threshold_mode == "elbow":
            stable, threshold, elbow_rank = _apply_stable_threshold(
                names_sort, freqs_sort, pi_thr, use_elbow=True, label=label
            )
        elif is_indep:
            use_elbow = (pi_thr_indep is None)
            thr       = pi_thr_indep if not use_elbow else pi_thr  # pi_thr unused when elbow
            stable, threshold, elbow_rank = _apply_stable_threshold(
                names_sort, freqs_sort, thr,
                use_elbow=use_elbow, label=label
            )
        else:
            stable, threshold, elbow_rank = _apply_stable_threshold(
                names_sort, freqs_sort, pi_thr,
                use_elbow=False, label=label
            )
        return names_sort, freq_map, stable, threshold, freqs_sort, elbow_rank

    # ─────────────────────────────────────────────────────────────────────────
    # MODE A: "split" — corr pool vs ind pool
    # ─────────────────────────────────────────────────────────────────────────
    if shadow_mode == "split":
        if indep_set and name_to_mod:
            corr_idx  = [i for i, f in enumerate(feature_names_surv)
                         if name_to_mod.get(f, "") not in indep_set]
            indep_idx = [i for i, f in enumerate(feature_names_surv)
                         if name_to_mod.get(f, "") in indep_set]
        else:
            # no group info — treat everything as corr, skip ind
            corr_idx  = list(range(len(feature_names_surv)))
            indep_idx = []

        # ── split+mDML: residualize corr pool using ind features ─────────────
        y_for_corr = y  # default (no override)
        if use_dml_residual and dml_scope == "all_modules" and len(indep_idx) > 0:
            X_ind_surv = X_surv_mat[:, indep_idx]
            if verbose: print(f"split+mDML: residualizing y for corr pool using "
                  f"{len(indep_idx)} ind features "
                  f"({dml_n_folds} folds, ntree={dml_ntree})")
            y_for_corr = _cross_fit_residuals(
                X_corr      = X_ind_surv,
                y           = y,
                k_folds     = dml_n_folds,
                ntree       = dml_ntree,
                mtry_factor = float(mtry_factor),
                nodesize    = int(nodesize),
                seed        = int(seed) + 88317 + 2,
                n_jobs      = int(n_jobs),
            )

        corr_override = y_for_corr if (use_dml_residual and dml_scope == "all_modules") else None
        names_corr,  freq_corr,  stable_corr,  thr_corr,  freqs_corr,  elbow_corr  = \
            _shadow_on_cols(corr_idx,  "corr",  seed_offset=0, is_indep=False,
                            y_override=corr_override)
        names_indep, freq_indep, stable_indep, thr_indep, freqs_indep, elbow_indep = \
            _shadow_on_cols(indep_idx, "indep", seed_offset=1, is_indep=True)

        freq_map_all = {**freq_corr, **freq_indep}
        stable_union = list(set(stable_corr) | set(stable_indep))

        # ── pool stability data for plotting ─────────────────────────────────
        pool_stability = []
        if names_corr:
            stable_set_c = set(stable_corr)
            pool_stability.append({
                'pool'        : 'corr',
                'feature_names': names_corr,
                'freqs'       : freqs_corr,
                'threshold'   : float(thr_corr),
                'selected'    : [f in stable_set_c for f in names_corr],
                'elbow_rank'  : elbow_corr,          # None for corr (fixed thr)
            })
        if names_indep:
            stable_set_i = set(stable_indep)
            pool_stability.append({
                'pool'        : 'indep',
                'feature_names': names_indep,
                'freqs'       : freqs_indep,
                'threshold'   : float(thr_indep),
                'selected'    : [f in stable_set_i for f in names_indep],
                'elbow_rank'  : elbow_indep,         # int if elbow was used
            })

    # ─────────────────────────────────────────────────────────────────────────
    # MODE B: "within_module" — one shadow run per module
    # ─────────────────────────────────────────────────────────────────────────
    elif shadow_mode == "within_module":
        pool_stability = []
        if not name_to_mod:
            # fallback: no module info — run as single pool (treat as corr)
            if verbose: print("within_module: no module assignments found, running as single pool")
            names_all, freq_map_all, stable_union, thr_all, freqs_all, elbow_all = \
                _shadow_on_cols(list(range(len(feature_names_surv))), "all",
                                seed_offset=0, is_indep=False)
            if names_all:
                stable_set_all = set(stable_union)
                pool_stability.append({
                    'pool'        : 'all',
                    'feature_names': names_all,
                    'freqs'       : freqs_all,
                    'threshold'   : float(thr_all),
                    'selected'    : [f in stable_set_all for f in names_all],
                    'elbow_rank'  : elbow_all,
                })
        else:
            # group features by module
            mod_to_idx = {}
            for i, f in enumerate(feature_names_surv):
                mod = name_to_mod.get(f, "__unknown__")
                mod_to_idx.setdefault(mod, []).append(i)

            freq_map_all = {}
            stable_union = []

            for m_offset, (mod_label, col_idx) in enumerate(mod_to_idx.items()):
                is_indep_mod = (mod_label in indep_set)

                # ── per-module DML: residualize y by all OTHER modules ────────
                y_mod = None   # None → _shadow_on_cols uses default logic
                if use_dml_residual and dml_scope == "all_modules":
                    other_idx = [i for i in range(len(feature_names_surv))
                                 if i not in col_idx]
                    if len(other_idx) == 0:
                        if verbose: print(f"mdml: no other-module features for "
                              f"module:{mod_label} — using original y")
                        y_mod = y
                    else:
                        X_other = X_surv_mat[:, other_idx]
                        if verbose: print(f"mdml: residualizing y for module:{mod_label} "
                              f"using {len(other_idx)} other-module features "
                              f"({dml_n_folds} folds, ntree={dml_ntree})")
                        y_mod = _cross_fit_residuals(
                            X_corr      = X_other,
                            y           = y,
                            k_folds     = dml_n_folds,
                            ntree       = dml_ntree,
                            mtry_factor = float(mtry_factor),
                            nodesize    = int(nodesize),
                            seed        = int(seed) + 88317 + m_offset * 997,
                            n_jobs      = int(n_jobs),
                        )

                names_mod, freq_mod, stable_mod, thr_mod, freqs_mod, elbow_mod = \
                    _shadow_on_cols(col_idx, f"module:{mod_label}",
                                    seed_offset=m_offset, is_indep=is_indep_mod,
                                    y_override=y_mod)
                freq_map_all.update(freq_mod)
                stable_union.extend(stable_mod)

                if names_mod:
                    stable_set_mod = set(stable_mod)
                    pool_stability.append({
                        'pool'        : mod_label,
                        'feature_names': names_mod,
                        'freqs'       : freqs_mod,
                        'threshold'   : float(thr_mod),
                        'selected'    : [f in stable_set_mod for f in names_mod],
                        'elbow_rank'  : elbow_mod,
                    })

            stable_union = list(set(stable_union))

    else:
        raise ValueError(f"shadow_mode must be 'split' or 'within_module', "
                         f"got '{shadow_mode}'")

    n_stable = len(stable_union)
    if verbose: print(f"stable_union: {n_stable} features"
          + (" (empty — returning empty result)" if n_stable == 0 else
             " -> single-shot RF + TreeSHAP"))

    # ── all-features sorted by freq (for selection_list compat) ──────────────
    all_sorted = sorted(freq_map_all.keys(),
                        key=lambda f: freq_map_all[f], reverse=True)
    all_freqs  = [freq_map_all[f] for f in all_sorted]

    # ── empty stable union: return empty (no padding) ─────────────────────────
    if n_stable == 0:
        return {
            'final_features'       : [],
            'final_vims'           : [],
            'stability_freq'       : [],
            'shap_matrix'          : [],
            'interaction_matrix'   : [],
            'pool_stability'       : pool_stability,
            'selection_list_feats' : [all_sorted, []],
            'selection_list_vims'  : [all_freqs,  []],
        }

    # ─────────────────────────────────────────────────────────────────────────
    # v18 core: single-shot RF on all stable features, then TreeSHAP
    # No K trim. No shortfall padding. What survived stability is the output.
    # ─────────────────────────────────────────────────────────────────────────

    # Sort stable_union by stability freq (descending) for consistent ordering
    stable_sorted = sorted(stable_union,
                           key=lambda f: freq_map_all.get(f, 0.0), reverse=True)

    # Build sub-matrix for stable features
    feat_idx_map = {f: i for i, f in enumerate(feature_names_surv)}
    stable_cols  = [feat_idx_map[f] for f in stable_sorted]
    X_stable     = X_surv_mat[:, stable_cols]

    ntree_s = _compute_ntree(n_stable, ntree_factor, min_ntree)
    mtry_s  = _compute_mtry(n_stable, mtry_factor, classification, mtry_rule)
    if verbose: print(f"Fitting single RF: n_stable={n_stable}  ntree={ntree_s}  mtry={mtry_s}")

    rf_stable = _build_rf(ntree_s, mtry_s, nodesize,
                          (seed + 8888) % (2**31 - 1), n_jobs, classification)
    rf_stable.fit(X_stable, y)

    # ── fastest TreeSHAP: observational explainer, no background ─────────────
    import shap as _shap
    explainer = _shap.TreeExplainer(rf_stable)
    shap_raw  = explainer.shap_values(X_stable, check_additivity=False)

    # ── extract signed SHAP matrix (n x n_stable) ────────────────────────────
    if isinstance(shap_raw, list):
        # classification: list of (n x p) arrays, one per class
        if classification and len(rf_stable.classes_) == 2:
            sv_signed = np.array(shap_raw[1])          # positive class
        elif classification:
            sv_signed = np.mean([np.array(s) for s in shap_raw], axis=0)
        else:
            sv_signed = np.array(shap_raw[0])          # regression
    else:
        if shap_raw.ndim == 3:                          # newer shap: (n, p, n_classes)
            if classification and shap_raw.shape[2] == 2:
                sv_signed = shap_raw[:, :, 1]
            elif classification:
                sv_signed = np.mean(shap_raw, axis=2)
            else:
                sv_signed = shap_raw[:, :, 0]
        else:
            sv_signed = np.array(shap_raw)             # (n, p)

    # mean |SHAP| per feature for ranking
    shap_mean_abs = np.mean(np.abs(sv_signed), axis=0)  # (n_stable,)

    # sort by mean |SHAP| descending
    order          = np.argsort(-shap_mean_abs)
    final_features = [stable_sorted[i] for i in order]
    final_vims     = [float(shap_mean_abs[i]) for i in order]
    stab_freqs     = [float(freq_map_all.get(f, 0.0)) for f in final_features]

    # signed SHAP matrix reordered to match final_features column order
    shap_matrix_ordered = sv_signed[:, order]            # (n, n_stable)

    # ── TreeSHAP interaction values (optional, exact) ─────────────────────────
    interaction_matrix_ordered = None
    if compute_interactions:
        try:
            if verbose: print(f"Computing TreeSHAP interaction values "
                  f"({n_stable}×{n_stable} matrix) …")
            int_mat_raw = _shap_interaction_matrix(rf_stable, X_stable, classification)
            # reorder rows and cols to match final_features order
            interaction_matrix_ordered = int_mat_raw[np.ix_(order, order)]
            if verbose: print("Interaction matrix complete.")
        except Exception as e:
            if verbose: print(f"Interaction values failed ({type(e).__name__}: {e}); "
                  f"skipping interaction matrix.")
            interaction_matrix_ordered = None

    if verbose: print(f"Single-shot SHAP complete: {n_stable} features returned"
          + (f"  top3={final_features[:3]}" if n_stable >= 3 else "")
          + (" [DML residual used for stability]" if use_dml_residual else ""))

    return {
        'final_features'       : list(final_features),
        'final_vims'           : list(final_vims),           # mean|SHAP|
        'stability_freq'       : list(stab_freqs),           # shadow boot freq
        'shap_matrix'          : shap_matrix_ordered.tolist(), # signed n x |U|
        'interaction_matrix'   : (interaction_matrix_ordered.tolist()
                                   if interaction_matrix_ordered is not None else []),
        'pool_stability'       : pool_stability,             # per-pool selection data
        'selection_list_feats' : [all_sorted, list(stable_sorted)],
        'selection_list_vims'  : [all_freqs,
                                   [freq_map_all.get(f, 0.0) for f in stable_sorted]],
    }


# ─────────────────────────────────────────────────────────────────────────────
# Public API — Phase 3: Final RF + exact SHAP  (identical to v12)
# ─────────────────────────────────────────────────────────────────────────────

def run_final_rf_shap(X_final_mat, y, feature_names_final,
                      n_estimators, mtry, nodesize,
                      classification, seed, n_jobs):
    X_final_mat         = np.array(X_final_mat, dtype=np.float64)
    y                   = np.array(y).ravel()
    feature_names_final = [str(f) for f in feature_names_final]
    n_estimators        = int(n_estimators)
    mtry                = int(mtry)
    nodesize            = int(nodesize)
    classification      = bool(classification)
    seed                = int(seed)
    n_jobs              = int(n_jobs)

    if classification:
        y = _map_labels(y)
    else:
        y = y.astype(np.float64)

    rf = _build_rf(n_estimators, mtry, nodesize, seed, n_jobs, classification)
    rf.fit(X_final_mat, y.astype(int) if classification else y.astype(float))

    mean_abs  = _shap_mean_abs(rf, X_final_mat, classification)
    shap_full = _shap_full_signed(rf, X_final_mat, classification)

    return {
        'mean_abs_shap' : mean_abs.tolist(),
        'shap_matrix'   : shap_full.tolist(),
        'feature_names' : feature_names_final,
    }
