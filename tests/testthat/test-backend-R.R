# End-to-end behaviour of the pure-R backend (ranger, no Python required).

test_that("pure-R backend (permutation) returns the documented result contract", {
  d <- make_sf_data()
  r <- fit_R(d, r_shap = "permutation")

  expect_s3_class(r, "mossy_forest")
  expect_setequal(
    names(r),
    c("final_rf", "final_X", "module_membership", "WGCNA_object",
      "survivor_list", "selection_list", "final_SHAP", "shap_obj",
      "runtimes", "feature_list", "stability_freq", "stability_data")
  )

  expect_s3_class(r$final_SHAP, "data.frame")
  expect_true(all(c("feature_name", "variable_importance", "module_membership")
                  %in% names(r$final_SHAP)))
  expect_gte(nrow(r$final_SHAP), 1L)

  stable <- r$final_SHAP$feature_name
  expect_true(all(stable %in% colnames(d$X)))          # never invents a feature
  expect_named(r$stability_freq)                       # named per-feature freq
})

test_that("screening keeps the true signals and can drop noise", {
  d <- make_sf_data()
  r <- fit_R(d, r_shap = "permutation")

  survivors <- unlist(lapply(r$survivor_list, function(df) {
    if (is.null(df) || nrow(df) == 0L) character(0) else as.character(df[[1L]])
  }), use.names = FALSE)

  # Strong signals should survive module-stratified screening.
  kept_signals <- intersect(survivors, d$true_signals)
  expect_gte(length(kept_signals), 5L)
})

test_that("the final ranger model can predict on new rows", {
  d <- make_sf_data()
  r <- fit_R(d, r_shap = "permutation")

  expect_false(is.null(r$final_rf))
  pr <- predict(r$final_rf, d$X[, r$final_SHAP$feature_name, drop = FALSE])$predictions
  expect_length(pr, nrow(d$X))
  expect_true(all(is.finite(pr)))
})

test_that("classification target runs end-to-end on the pure-R backend", {
  d <- make_sf_data()
  r <- mf(d$X, d$y_clf, module_membership = d$mm,
          screen_params = fast_screen(), select_params = fast_select(),
          final_ntree = 60L, num_processors = 1L, verbose = 0L, seed = 3L,
          backend = "R", r_shap = "permutation")
  expect_s3_class(r, "mossy_forest")
  expect_gte(nrow(r$final_SHAP), 1L)
})

test_that("r_shap = 'treeshap' drives screening + stability + final without error", {
  skip_if_not_installed("treeshap")
  d <- make_sf_data()
  r <- fit_R(d, r_shap = "treeshap")

  expect_s3_class(r, "mossy_forest")
  expect_gte(nrow(r$final_SHAP), 1L)
  expect_true(all(r$final_SHAP$feature_name %in% colnames(d$X)))
})

test_that("cfRes residualization on the unassigned pool runs (pure-R)", {
  d <- make_sf_data()
  r <- mf(d$X, d$y_reg, module_membership = d$mm,
          screen_params = fast_screen(),
          select_params = fast_select(use_cfres = TRUE,
                                       cfres_n_folds = 3L, cfres_ntree = 80L),
          final_ntree = 60L, num_processors = 1L, verbose = 0L, seed = 3L,
          backend = "R", r_shap = "permutation")
  expect_s3_class(r, "mossy_forest")
})

# ── directional (signed) importance on the pure-R backend ───────────────────
# Direction needs a per-observation SHAP matrix, so it is only available for
# r_shap = "fastshap" / "treeshap" — never "permutation" (a scalar VIM has
# nothing to correlate against). See plot_signed_importance().

direction_cols <- c("direction", "dir_corr", "signed_importance",
                     "direction_ambiguous")

test_that("permutation VIM never produces direction columns", {
  d <- make_sf_data()
  r <- fit_R(d, r_shap = "permutation")
  expect_false(any(direction_cols %in% names(r$final_SHAP)))
})

test_that("fastshap produces direction columns for regression", {
  skip_if_not_installed("fastshap")
  d <- make_sf_data()
  r <- fit_R(d, r_shap = "fastshap")

  expect_true(all(direction_cols %in% names(r$final_SHAP)))
  fs <- r$final_SHAP
  expect_true(all(fs$direction %in% c(-1L, 0L, 1L)))
  expect_true(all(abs(fs$dir_corr) <= 1))
  expect_equal(fs$direction_ambiguous, abs(fs$dir_corr) < 0.1)
  expect_equal(fs$signed_importance, fs$direction * fs$mean_abs_shap,
               tolerance = 1e-6)
})

test_that("treeshap produces direction columns for regression", {
  skip_if_not_installed("treeshap")
  d <- make_sf_data()
  r <- fit_R(d, r_shap = "treeshap")

  expect_true(all(direction_cols %in% names(r$final_SHAP)))
  fs <- r$final_SHAP
  expect_true(all(fs$direction %in% c(-1L, 0L, 1L)))
  expect_equal(fs$signed_importance, fs$direction * fs$mean_abs_shap,
               tolerance = 1e-6)
})

test_that("fastshap produces direction columns for classification too", {
  skip_if_not_installed("fastshap")
  d <- make_sf_data()
  r <- mf(d$X, d$y_clf, module_membership = d$mm,
          screen_params = fast_screen(), select_params = fast_select(),
          final_ntree = 60L, num_processors = 1L, verbose = 0L, seed = 3L,
          backend = "R", r_shap = "fastshap")
  expect_true(all(direction_cols %in% names(r$final_SHAP)))
})

test_that("treeshap on classification fails upstream (known treeshap::ranger.unify limitation)", {
  # treeshap::ranger.unify() does not support ranger probability forests
  # (probability = TRUE, which mf() always uses for classification): it looks
  # for a "prediction" column that classification trees don't have. This is a
  # limitation of the treeshap package, not of direction detection here. If a
  # future treeshap release adds support, this test will start failing loudly
  # -- that's the intended signal to revisit r_shap = "treeshap" + classification.
  skip_if_not_installed("treeshap")
  d <- make_sf_data()
  expect_error(
    mf(d$X, d$y_clf, module_membership = d$mm,
       screen_params = fast_screen(), select_params = fast_select(),
       final_ntree = 60L, num_processors = 1L, verbose = 0L, seed = 3L,
       backend = "R", r_shap = "treeshap")
  )
})
