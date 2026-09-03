# mf() argument surface: input validation and the compute_interactions plumbing.
# The Python-backend behaviour is only exercised where that backend is available.

test_that("mf() validates its core inputs", {
  d <- make_sf_data(n = 60)
  expect_error(
    mf(as.matrix(d$X), d$y_reg, module_membership = d$mm, backend = "R"),
    "X must be a data.frame"
  )
  expect_error(
    mf(d$X, d$X, module_membership = d$mm, backend = "R"),
    "y must be"
  )
  expect_error(mf(d$X, d$y_reg, module_membership = d$mm, backend = "nope"))
  expect_error(mf(d$X, d$y_reg, module_membership = d$mm,
                  backend = "R", r_shap = "nope"))
})

test_that("compute_interactions toggles cleanly on the Python backend", {
  skip_if_no_py_backend()
  d <- make_sf_data()

  run <- function(ci) {
    mf(d$X, d$y_reg, module_membership = d$mm,
       screen_params = fast_screen(), select_params = fast_select(),
       final_ntree = 60L, num_processors = 1L, verbose = 0L, seed = 3L,
       backend = "python", compute_interactions = ci)
  }

  r_on  <- run(TRUE)
  r_off <- run(FALSE)
  expect_s3_class(r_on,  "mossy_forest")
  expect_s3_class(r_off, "mossy_forest")

  # Same selection regardless of the interaction toggle.
  expect_setequal(r_on$final_SHAP$feature_name, r_off$final_SHAP$feature_name)
})

test_that("Python and pure-R backends agree on the screening survivors", {
  skip_if_no_py_backend()
  d <- make_sf_data()

  common <- list(
    module_membership = d$mm, screen_params = fast_screen(),
    select_params = fast_select(), final_ntree = 60L,
    num_processors = 1L, verbose = 0L, seed = 3L
  )
  r_py <- do.call(mf, c(list(d$X, d$y_reg), common, backend = "python"))
  r_r  <- do.call(mf, c(list(d$X, d$y_reg), common, backend = "R",
                        r_shap = "permutation"))

  # Different SHAP engines pick different final sets, but both must keep the
  # strong signals available after screening.
  surv <- function(r) unlist(lapply(r$survivor_list, function(df)
    if (is.null(df) || nrow(df) == 0L) character(0) else as.character(df[[1L]])),
    use.names = FALSE)
  expect_gte(length(intersect(surv(r_py), d$true_signals)), 5L)
  expect_gte(length(intersect(surv(r_r),  d$true_signals)), 5L)
})

test_that("classification direction is w.r.t. levels(y)[1], not factor position (regression guard)", {
  # Guards a real bug: as.integer(y) on an R factor returns each level's 1-based
  # POSITION (1, 2, ...), not a 0/1 label. That sent sklearn classes_ = [1, 2]
  # instead of [0, 1], and every hardcoded classes_[1] pick in
  # mf_python_backend.py then explained levels(y)[2] instead of levels(y)[1] --
  # silently flipping the sign of direction / dir_corr / signed_importance for
  # every classification result. This factor is built exactly like the bug
  # report: levels = c(1, 0), labels = c("case", "control"), so "case" is
  # level 1 -- as.integer(y) would give 1 for case / 2 for control, the
  # data-independent shape that used to trigger the flip.
  skip_if_no_py_backend()
  set.seed(42)
  n <- 400
  f1 <- rnorm(n)                        # true driver: higher f1 -> more likely "case"
  f2 <- rnorm(n)                        # noise
  casecontrol <- rbinom(n, 1, plogis(1.5 * f1))
  y <- factor(casecontrol, levels = c(1, 0), labels = c("case", "control"))
  expect_identical(levels(y), c("case", "control"))
  expect_identical(as.integer(y)[1:3] %in% 1:2, rep(TRUE, 3))  # position-coded, not 0/1

  X  <- data.frame(f1 = f1, f2 = f2)
  mm <- c(f1 = "M", f2 = "M")

  r <- mf(X, y, module_membership = mm,
          screen_params = screen_control(min_ntree = 100L),
          select_params = select_control(n_boots = 10L, min_ntree = 100L, pi_thr = 0.3),
          final_ntree = 200L, verbose = 0L, seed = 1L, backend = "python")

  f1_row <- r$final_SHAP[r$final_SHAP$feature_name == "f1", ]
  expect_equal(nrow(f1_row), 1L)
  # f1 increases P(case) = levels(y)[1] by construction -> direction must be +1.
  expect_equal(f1_row$direction, 1L)
  expect_gt(f1_row$dir_corr, 0)
})
