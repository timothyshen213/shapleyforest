# data-raw/prep_wdi.R
#
# Run this script ONCE to create the bundled `wdi` dataset.
# Requires: install.packages(c("WDI", "data.table"))
#
# A non-biological, genuinely high-dimensional (p > n), highly-correlated
# tabular dataset — useful for demonstrating Mossy Forest on data whose
# features (GDP per capita, literacy, CO2 emissions, ...) are intuitive to a
# general audience while posing the same challenge as genomic panels.
#
# Saves a single list `wdi` to data/, with elements:
#   wdi$X      — 181 × 501 numeric matrix (World Bank indicators, 2020)
#   wdi$y      — numeric vector length 181 (life expectancy at birth, years)
#   wdi$labels — named chr: indicator code -> human-readable indicator name
if (!requireNamespace("WDI", quietly = TRUE))        install.packages("WDI")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

YEAR      <- "2020"
IND_COVER <- 0.80   # keep indicators present for >= 80% of countries
CTY_COVER <- 0.80   # keep countries with >= 80% of the retained indicators
TARGET    <- "SP.DYN.LE00.IN"   # life expectancy at birth, total (years)

# ── Download the World Bank WDI bulk CSVs ───────────────────────────────────────
tmp <- tempfile(fileext = ".zip")
# HTTP/1.1 avoids intermittent HTTP/2 stream errors on this large (~270 MB) file
utils::download.file(
  "https://databankfiles.worldbank.org/public/ddpext_download/WDI_CSV.zip",
  tmp, mode = "wb", method = "curl", extra = "--http1.1 --retry 3")
ex <- tempfile(); dir.create(ex)
utils::unzip(tmp, files = c("WDICSV.csv", "WDICountry.csv"), exdir = ex)

dat <- fread(file.path(ex, "WDICSV.csv"),     header = TRUE)
cty <- fread(file.path(ex, "WDICountry.csv"), header = TRUE)

# ── Reshape one year to a country x indicator matrix ────────────────────────────
# Real countries only: World Bank aggregates ("World", "High income", ...) have
# a blank Region field.
real   <- cty[Region != "" & !is.na(Region), `Country Code`]
labels <- unique(dat[, .(code = `Indicator Code`, name = `Indicator Name`)])
cname  <- unique(cty[, .(code = `Country Code`, nm = `Short Name`)])

dyr <- dat[`Country Code` %in% real,
           .(CountryCode = `Country Code`, Indicator = `Indicator Code`,
             val = get(YEAR))]
w   <- dcast(dyr, CountryCode ~ Indicator, value.var = "val")
mat <- as.matrix(w[, -1]); rownames(mat) <- w$CountryCode

# ── Outcome ─────────────────────────────────────────────────────────────────────
y_all <- mat[, TARGET]

# ── Feature filter ──────────────────────────────────────────────────────────────
keep_ind <- colMeans(!is.na(mat)) >= IND_COVER
# Drop tautological life-table-derived series (life expectancy, mortality,
# survival, deaths). These are mechanically part of the outcome, not
# independent predictors — matched on the human-readable indicator name.
lab_name  <- setNames(labels$name, labels$code)
taut_rx   <- "life expectancy|mortalit|surviv|death|dying|stillbirth|perinatal|health pillar"
taut_code <- names(lab_name)[grepl(taut_rx, lab_name, ignore.case = TRUE)]
keep_ind[colnames(mat) %in% taut_code] <- FALSE
X <- mat[, keep_ind, drop = FALSE]

# ── Country filter: has target + enough retained indicators ─────────────────────
rok <- !is.na(y_all) & (rowMeans(!is.na(X)) >= CTY_COVER)
X   <- X[rok, , drop = FALSE]
y   <- y_all[rok]

# ── Median-impute residual gaps; drop near-constant indicators ──────────────────
for (j in seq_len(ncol(X))) {
  v <- X[, j]; if (anyNA(v)) X[is.na(v), j] <- median(v, na.rm = TRUE)
}
nzv <- apply(X, 2, function(v) length(unique(v)) > 1 && sd(v) > 0)
X   <- X[, nzv, drop = FALSE]

# ── Readable names ──────────────────────────────────────────────────────────────
rownames(X) <- make.unique(cname[match(rownames(X), code), nm]); names(y) <- rownames(X)
lab <- setNames(labels[match(colnames(X), code), name], colnames(X))

wdi <- list(X = X, y = y, labels = lab)

cat("wdi$X :", nrow(wdi$X), "x", ncol(wdi$X), "\n")
cat("wdi$y : life expectancy", paste(round(range(wdi$y), 1), collapse = " - "), "years\n")

# ── Save ────────────────────────────────────────────────────────────────────────
usethis::use_data(wdi, overwrite = TRUE)
