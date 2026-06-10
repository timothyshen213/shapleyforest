# data-raw/prep_leukemia.R
#
# Run this script ONCE to create the bundled leukemia dataset.
# Requires: BiocManager::install(c("ALL", "Biobase"))
#
# Saves a single list `leukemia` to data/, with elements:
#   leukemia$X  — n × 500 numeric matrix (top-500 variable probe-sets, log2)
#   leukemia$y  — integer vector length n (1 = BCR/ABL, 0 = NEG)
#
# Contrast: among B-lineage ALL samples, the BCR/ABL (Philadelphia-chromosome,
# Ph+) molecular subtype vs cytogenetically normal (NEG).  This is the standard
# *balanced* ALL benchmark used throughout the Bioconductor "Case Studies"
# (~37 BCR/ABL vs ~42 NEG), replacing the heavily imbalanced B-vs-T lineage
# split (95 / 33).
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in c("ALL", "Biobase"))
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)

library(ALL)
library(Biobase)

data(ALL)

# ── Sample selection ───────────────────────────────────────────────────────────
# B-lineage only, and only the two balanced molecular subtypes of interest.
mol      <- as.character(ALL$mol.biol)
is_bcell <- grepl("^B", ALL$BT)
keep     <- is_bcell & mol %in% c("BCR/ABL", "NEG")

ALLsub <- ALL[, keep]

# ── Outcome ────────────────────────────────────────────────────────────────────
y <- as.integer(as.character(ALLsub$mol.biol) == "BCR/ABL")   # 1 = BCR/ABL, 0 = NEG
names(y) <- sampleNames(ALLsub)

# ── Feature matrix ─────────────────────────────────────────────────────────────
expr     <- t(exprs(ALLsub))                      # n × 12625, already log2
gene_var <- apply(expr, 2, var)
top500   <- order(gene_var, decreasing = TRUE)[seq_len(500L)]
X        <- expr[, top500]

leukemia <- list(X = X, y = y)

cat("leukemia$X :", nrow(leukemia$X), "x", ncol(leukemia$X), "\n")
cat("leukemia$y :", paste(names(table(leukemia$y)), table(leukemia$y),
                          sep = "=", collapse = "  "), "\n")

# ── Save ───────────────────────────────────────────────────────────────────────
usethis::use_data(leukemia, overwrite = TRUE)
