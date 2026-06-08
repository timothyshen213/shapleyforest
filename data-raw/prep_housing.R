# data-raw/prep_housing.R
#
# Run this script ONCE to create the bundled housing dataset.
# Requires: install.packages("modeldata")
#
# A non-biological, high-dimensional, highly-correlated tabular dataset —
# useful for demonstrating Shapley Forest on data whose features (square
# footage, room counts, garage size, ...) are intuitive to a general audience.
#
# Saves a single list `housing` to data/, with elements:
#   housing$X  — 2930 × 33 numeric matrix (Ames property attributes)
#   housing$y  — numeric vector length 2930 (sale price, USD)
if (!requireNamespace("modeldata", quietly = TRUE))
  install.packages("modeldata")

data(ames, package = "modeldata")

# ── Outcome ────────────────────────────────────────────────────────────────────
y <- as.numeric(ames$Sale_Price)
names(y) <- paste0("house_", seq_along(y))

# ── Feature matrix ─────────────────────────────────────────────────────────────
# Keep the numeric property attributes (drop Sale_Price itself — that's `y`)
is_num <- vapply(ames, is.numeric, logical(1))
is_num[["Sale_Price"]] <- FALSE

X <- as.matrix(ames[, is_num])
mode(X) <- "double"
rownames(X) <- names(y)

housing <- list(X = X, y = y)

cat("housing$X :", nrow(housing$X), "x", ncol(housing$X), "\n")
cat("housing$y : range $", paste(round(range(housing$y)), collapse = " - $"), "\n")

# ── Save ───────────────────────────────────────────────────────────────────────
usethis::use_data(housing, overwrite = TRUE)
