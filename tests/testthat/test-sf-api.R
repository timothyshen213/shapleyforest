# bf() argument surface: input validation and the compute_interactions plumbing.
# The Python-backend behaviour is only exercised where that backend is available.

test_that("bf() validates its core inputs", {
  d <- make_sf_data(n = 60)
  expect_error(
    bf(as.matrix(d$X), d$y_reg, module_membership = d$mm, backend = "R"),
    "X must be a data.frame"
  )
  expect_error(
    bf(d$X, d$X, module_membership = d$mm, backend = "R"),
    "y must be"
  )
  expect_error(bf(d$X, d$y_reg, module_membership = d$mm, backend = "nope"))
  expect_error(bf(d$X, d$y_reg, module_membership = d$mm,
                  backend = "R", r_shap = "nope"))
})

test_that("compute_interactions toggles cleanly on the Python backend", {
  skip_if_no_py_backend()
  d <- make_sf_data()

  run <- function(ci) {
    bf(d$X, d$y_reg, module_membership = d$mm,
       screen_params = fast_screen(), select_params = fast_select(),
       final_ntree = 60L, num_processors = 1L, verbose = 0L, seed = 3L,
       backend = "python", compute_interactions = ci)
  }

  r_on  <- run(TRUE)
  r_off <- run(FALSE)
  expect_s3_class(r_on,  "bonsai_forest")
  expect_s3_class(r_off, "bonsai_forest")

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
  r_py <- do.call(bf, c(list(d$X, d$y_reg), common, backend = "python"))
  r_r  <- do.call(bf, c(list(d$X, d$y_reg), common, backend = "R",
                        r_shap = "permutation"))

  # Different SHAP engines pick different final sets, but both must keep the
  # strong signals available after screening.
  surv <- function(r) unlist(lapply(r$survivor_list, function(df)
    if (is.null(df) || nrow(df) == 0L) character(0) else as.character(df[[1L]])),
    use.names = FALSE)
  expect_gte(length(intersect(surv(r_py), d$true_signals)), 5L)
  expect_gte(length(intersect(surv(r_r),  d$true_signals)), 5L)
})
