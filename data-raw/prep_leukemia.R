# data-raw/prep_leukemia.R
#
# Run this script ONCE to create the bundled leukemia dataset.
# Requires: BiocManager::install(c("ALL", "Biobase"))
#
# Saves a single list `leukemia` to data/, with elements:
#   leukemia$X  — 128 × 500 numeric matrix (top-500 variable probe-sets, log2)
#   leukemia$y  — integer vector length 128 (1 = B-cell, 0 = T-cell)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALL")
BiocManager::install("Biobase")

library(ALL)
library(Biobase)

data(ALL)

# ── Outcome ────────────────────────────────────────────────────────────────────
y <- as.integer(grepl("^B", ALL$BT))
names(y) <- sampleNames(ALL)

# ── Feature matrix ─────────────────────────────────────────────────────────────
expr     <- t(exprs(ALL))                         # 128 × 12625, already log2
gene_var <- apply(expr, 2, var)
top500   <- order(gene_var, decreasing = TRUE)[seq_len(500L)]
X        <- expr[, top500]

leukemia <- list(X = X, y = y)

cat("leukemia$X :", nrow(leukemia$X), "x", ncol(leukemia$X), "\n")
cat("leukemia$y :", table(leukemia$y), "\n")

# ── Save ───────────────────────────────────────────────────────────────────────
usethis::use_data(leukemia, overwrite = TRUE)
