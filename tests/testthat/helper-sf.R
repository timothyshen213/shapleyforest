# Shared fixtures and skip helpers for the mossyforest test suite.

# A small, fast, module-structured data set with known true signals:
#   module A  : A1..A4  correlated with latent z1
#   module B  : B1..B4  correlated with latent z2
#   independent: indSIG (real signal, no module partners) + noise1..noise4
# Regression target depends on z1 - z2 + indSIG, so every "true" feature below
# carries real signal and every noise* column is pure noise.
make_sf_data <- function(n = 160, sd_within = 0.30, sd_noise = 0.35, seed = 11) {
  set.seed(seed)
  z1  <- rnorm(n)
  z2  <- rnorm(n)
  ind <- rnorm(n)
  cols <- c(
    lapply(1:4, function(i) z1 + rnorm(n, sd = sd_within)),
    lapply(1:4, function(i) z2 + rnorm(n, sd = sd_within)),
    list(ind),
    lapply(1:4, function(i) rnorm(n))
  )
  nm <- c(paste0("A", 1:4), paste0("B", 1:4), "indSIG", paste0("noise", 1:4))
  X  <- as.data.frame(setNames(cols, nm))
  mm <- setNames(
    c(rep("A", 4), rep("B", 4), rep("unassigned", 5)),
    nm
  )
  y_reg <- as.numeric(z1 - z2 + 1.5 * ind + rnorm(n, sd = sd_noise))
  y_clf <- factor(as.integer(y_reg > median(y_reg)))
  list(
    X = X, mm = mm, y_reg = y_reg, y_clf = y_clf,
    true_signals = c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "indSIG"),
    noise = paste0("noise", 1:4)
  )
}

# Fast control objects so end-to-end tests stay well under a second each.
fast_screen <- function(...) screen_control(min_ntree = 60L, ...)
fast_select <- function(...) {
  select_control(n_boots = 10L, min_ntree = 60L, pi_thr = 0.30,
                 shadow_mode = "split", unassigned_modules = "unassigned", ...)
}

# Run a pure-R backend fit with the fast controls.
fit_R <- function(d, r_shap = "permutation", ...) {
  mf(d$X, d$y_reg, module_membership = d$mm,
     screen_params = fast_screen(), select_params = fast_select(),
     final_ntree = 60L, num_processors = 1L, verbose = 0L, seed = 3L,
     backend = "R", r_shap = r_shap, ...)
}

# Skip a test unless the Python backend (reticulate + sklearn + shap) is usable.
skip_if_no_py_backend <- function() {
  testthat::skip_if_not_installed("reticulate")
  ok <- tryCatch({
    reticulate::import("sklearn", delay_load = FALSE)
    reticulate::import("shap",    delay_load = FALSE)
    reticulate::import("numpy",   delay_load = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(isTRUE(ok), "Python backend (sklearn/shap/numpy) not available")
}
