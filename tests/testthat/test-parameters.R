# Guards the R control objects against drift from the Python reference package
# (mossyforest_py). If someone reorders or renames a field on one side, these
# fail loudly instead of the divergence going unnoticed. The expected orders
# below are copied verbatim from src/mossyforest/parameters.py.

test_that("screen_control field order matches the Python reference", {
  expected <- c("drop_fraction", "keep_fraction", "mtry_factor",
                "mtry_rule", "min_ntree", "ntree_factor")
  expect_identical(names(formals(screen_control)), expected)
  expect_identical(names(screen_control()), expected)
})

test_that("select_control field order matches the Python reference", {
  expected <- c("drop_fraction", "mtry_factor", "mtry_rule", "mtry_on_real",
                "min_ntree", "ntree_factor", "shadow_mode", "n_boots",
                "pi_thr", "pi_thr_unassigned", "threshold_mode", "shadow_percentile",
                "unassigned_modules", "use_cfres", "cfres_n_folds", "cfres_ntree",
                "cfres_scope", "early_stop_boots", "early_stop_tol",
                "early_stop_check_every")
  expect_identical(names(formals(select_control)), expected)
  expect_identical(names(select_control()), expected)
})

test_that("control defaults match the Python reference", {
  sc <- screen_control()
  expect_equal(sc$drop_fraction, 0.25)
  expect_equal(sc$keep_fraction, 0.50)
  expect_equal(sc$mtry_rule, "auto")
  expect_equal(sc$min_ntree, 500L)

  sl <- select_control()
  expect_equal(sl$drop_fraction, 0.10)
  expect_equal(sl$pi_thr, 0.60)
  expect_equal(sl$shadow_mode, "split")
  expect_equal(sl$threshold_mode, "pi_thr")
  expect_equal(sl$cfres_scope, "unassigned")
  expect_false(sl$mtry_on_real)
  expect_false(sl$use_cfres)
  expect_false(sl$early_stop_boots)
  expect_null(sl$pi_thr_unassigned)
  expect_null(sl$unassigned_modules)
})

test_that("control objects carry the right S3 class", {
  expect_s3_class(screen_control(), "screen_control")
  expect_s3_class(select_control(), "select_control")
})

test_that("invalid enum arguments are rejected", {
  expect_error(screen_control(mtry_rule = "nope"))
  expect_error(select_control(mtry_rule = "nope"))
  expect_error(select_control(shadow_mode = "nope"))
  expect_error(select_control(threshold_mode = "nope"))
  expect_error(select_control(cfres_scope = "nope"))
  expect_error(select_control(ntree_factor = 0))
})
