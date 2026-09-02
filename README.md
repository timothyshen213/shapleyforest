# `mossyforest`: Mossy Forest

Network-aware stable feature selection via Shapley values and random forests.

**UNDER DEVELOPMENT R Package.** `mossyforest` selects a small, stable set of
features from high-dimensional, correlated data by combining:

1. **Module-stratified screening** — recursive feature elimination (RFE) within
   each correlated module (or WGCNA-derived module).
2. **Shadow-stability selection** — Boruta-style bootstrap selection where each
   feature must beat its own permuted "shadow" across many resamples.
3. **A single TreeSHAP pass** on the stable set for importance, confidence
   intervals, and a pairwise interaction matrix.

A [Python implementation](https://github.com/timothyzshen/mossyforest) mirrors
this package and is bit-exact with the `backend = "python"` path (see
[Parity](#parity-with-the-python-package)); prefer it for very large data sets.

---

## Installation

```r
if (!require("devtools", quietly = TRUE))
       install.packages("devtools")
library(devtools)
devtools::install_github("timothyzshen/mossyforest")
```

The default `backend = "python"` uses `reticulate` and needs a Python
environment with `scikit-learn`, `shap`, and `numpy`. Run `mf_setup()` once to
provision it. The `backend = "R"` path is pure R (`ranger`) and needs no Python;
`r_shap = "treeshap"` or `"fastshap"` additionally needs the matching R package.

### A note on downloading `WGCNA`

`WGCNA` (used by `wmf()` for automatic module detection) installs from
Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install(c("GO.db", "preprocessCore", "impute"))  # dependencies
install.packages("WGCNA")
```

---

## Quick start

```r
library(mossyforest)

# X: data.frame of samples x features; y: numeric (regression) or factor (classification)
res <- mf(
  X, y,
  module_membership = c(gene1 = "M1", gene2 = "M1", gene3 = "M2", ...),
  screen_params     = screen_control(),   # Phase 1 defaults
  select_params     = select_control(),   # Phase 2 defaults
  num_processors    = 4,                   # parallel bootstraps / modules
  seed              = 42
)

res$final_SHAP        # stable features: importance, mean|SHAP|, 95% CIs, freq
res$stability_freq    # per-feature shadow-bootstrap selection frequency
detect_interaction(res, thresh = 0.01)    # ranked pairwise TreeSHAP interactions
```

### Automatic WGCNA module detection

```r
res <- wmf(X, y, WGCNA_params = WGCNA_control(power = 6, min_features = 20))
res$WGCNA_object      # raw module labels
```

---

## Backends

| `backend` | engine | SHAP for the final step | notes |
|---|---|---|---|
| `"python"` (default) | scikit-learn RF | TreeSHAP (`shap`) | bit-exact with the Python package; needs `reticulate` |
| `"R"` | `ranger` | set by `r_shap` | pure R, no Python |

For `backend = "R"`, `r_shap` selects how features are ranked in **screening,
stability selection, and the final step**:

- `"permutation"` (default) — ranger permutation importance; fastest.
- `"treeshap"` — TreeSHAP throughout (needs the `treeshap` package; regression
  and binary classification only). Slower but SHAP-consistent.
- `"fastshap"` — model-agnostic SHAP via the `fastshap` package.

The pure-R backend is *structurally* aligned with the Python package but cannot
be numerically identical to it (ranger vs scikit-learn, and a different SHAP
engine).

---

## Parameters

### `screen_control` — Phase 1 (screening)
| param | default | meaning |
|---|---|---|
| `drop_fraction` | 0.25 | fraction of features dropped per RFE step |
| `keep_fraction` | 0.50 | fraction of each module retained after screening |
| `mtry_factor` | 1.0 | `mtry` multiplier |
| `mtry_rule` | `"auto"` | `mtry` rule: `"auto"` (p/3 clf, √p reg), `"sqrt"`, `"p_over_3"` |
| `min_ntree` | 500 | minimum trees per forest |
| `ntree_factor` | 1 | trees = `max(p * factor, min_ntree)` |

### `select_control` — Phase 2 (shadow stability)
| param | default | meaning |
|---|---|---|
| `drop_fraction` | 0.10 | fraction dropped per survivor-pool RFE step |
| `mtry_factor` | 1.0 | `mtry` multiplier |
| `mtry_rule` | `"auto"` | see `screen_control` |
| `mtry_on_real` | `FALSE` | size `mtry` on real features only (not the shadows) |
| `min_ntree` / `ntree_factor` | 500 / 1 | tree-count controls |
| `shadow_mode` | `"split"` | `"split"` (corr + unassigned pools) or `"within_module"` |
| `n_boots` | 50 | bootstrap replicates |
| `pi_thr` | 0.60 | stability threshold for the correlated pool |
| `pi_thr_unassigned` | `NULL` | unassigned-pool threshold (`NULL` = elbow detection) |
| `threshold_mode` | `"pi_thr"` | `"pi_thr"` (R default) or `"elbow"` (both pools) |
| `shadow_percentile` | 95 | percentile of the shadow null used as the cutoff |
| `unassigned_modules` | `NULL` | module names treated as the unassigned pool |
| `use_cfres` | `FALSE` | cross-fitted cfRes residualization of `y` |
| `cfres_n_folds` / `cfres_ntree` | 5 / 300 | cfRes cross-fitting controls |
| `cfres_scope` | `"unassigned"` | `"unassigned"` or `"all_modules"` |
| `early_stop_boots` | `FALSE` | stop bootstraps early once frequencies converge |
| `early_stop_tol` / `early_stop_check_every` | 0.01 / 10 | early-stop controls |

### `mf()` execution
| param | default | meaning |
|---|---|---|
| `backend` | `"python"` | `"python"` or `"R"` (see [Backends](#backends)) |
| `r_shap` | `"permutation"` | SHAP method for `backend = "R"` |
| `num_processors` | 1 | parallel workers for bootstraps + modules |
| `compute_interactions` | `TRUE` | retain the RF so interactions compute lazily |
| `final_ntree` | — | trees in the final model |
| `seed` | — | RNG seed (matches the Python package's scheme) |

---

## Result object (class `mossy_forest`)

| element | description |
|---|---|
| `final_SHAP` | stable features: importance, mean\|SHAP\|, 95% CIs, stability freq |
| `shap_obj$shapley_values` | signed n×p SHAP matrix |
| `shap_obj$interaction_matrix` | p×p mean \|TreeSHAP interaction\| |
| `stability_freq` | per-feature shadow-bootstrap selection frequency |
| `stability_data` | per-pool stability curves (for `plot_stability_elbow`) |
| `survivor_list` | Phase-1 survivors per module |
| `final_rf` / `final_X` | the fitted final forest and its design matrix |
| `runtimes` | seconds for Screen / Selection / Final_RF |

---

## Plots

```r
plot_importance(res)              # beeswarm / bar SHAP importance
plot_signed_importance(res)
plot_waterfall(res); plot_force(res); plot_decisions(res)
plot_stability(res); plot_stability_elbow(res)     # per-pool selection curves
plot_potential_interactions(res)  # TreeSHAP interaction heatmap
plot_modules(res); plot_elbow(res)
```

---

## Parity with the Python package

The `backend = "python"` path and the Python `mossyforest` package are a
bit-exact match on identical data, seed, and stack: screening survivors and the
selected set are identical, stability frequencies agree to 0.0000, and the
signed SHAP matrix matches (max \|Δ\| = 0.000). Re-verify with
`benchmarks/parity_in_process.py` in the Python repository.

---

## License

GPL-3.0-or-later.
