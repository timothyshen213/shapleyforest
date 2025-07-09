library(fuzzyforest)
library(WGCNA)
library(future)
library(ranger)
options(stringsAsFactors = FALSE)

# Set seed for reproducibility
set.seed(89)

# ---- Step 1: Simulate gene expression data ----
n_samples <- 100   # number of samples
n_genes <- 1000    # number of genes

# Simulate correlated gene modules
simulate_module <- function(n_samples, n_genes, rho = 0.8) {
  cov_mat <- matrix(rho, nrow = n_genes, ncol = n_genes)
  diag(cov_mat) <- 1
  t(MASS::mvrnorm(n = n_samples, mu = rep(0, n_genes), Sigma = cov_mat))
}

# Create 5 modules of 200 genes each
modules <- lapply(1:5, function(i) simulate_module(n_samples, 200, rho = 0.7))
datExpr <- do.call(rbind, modules)  # genes x samples
datExpr <- t(datExpr)               # transpose to samples x genes
colnames(datExpr) <- paste0("Gene", 1:ncol(datExpr))

# ---- Step 2: Create a classification Y variable ----
# Let class be determined by expression of first module
mod1_activity <- rowMeans(datExpr[, 1:50])
Y <- ifelse(mod1_activity > median(mod1_activity), "High", "Low")
Y <- factor(Y)

# Check for good genes/samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodGenes, gsg$goodSamples]
}

# Pick soft-threshold power
powers <- c(1:10)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 3

# Run WGCNA blockwiseModules to get modules
net <- blockwiseModules(datExpr,
                        power = softPower,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)

moduleLabels <- net$colors
table(moduleLabels)

module_membership <- factor(moduleLabels)

# Transpose back to samples x genes for fuzzy forest
X_for_fuzzy <-as.data.frame((datExpr))

# Set fuzzy forest parameters
screen_params <- screen_control(
  drop_fraction = 0.5,
  keep_fraction = 0.25,
  mtry_factor = 0.25
)

select_params <- select_control(
  number_selected = 10,
  drop_fraction = 0.5,
  mtry_factor = 0.25
)

n_cores <- parallel::detectCores() - 1

set.seed(89)

Y <- as.factor(Y)

start <- proc.time()

sff <- shapff(X_for_fuzzy, Y, module_membership = moduleLabels,
              select_params = select_params,
              shap_model= "full", 
              screen_params = screen_params, 
              auto_initial = 2,
              nodesize = 1, 
              debug = 1, 
              verbose = 1, 
              initial = TRUE, 
              num_processors = n_cores, 
              parallel = 1,
              min_features = 10)

end <- proc.time()
par_run <- end - start 

library(shapleyforest)

start <- proc.time()

sff <- shapleyforest::shapff(X_for_fuzzy, Y, module_membership = moduleLabels,
              select_params = select_params,
              shap_model= "full", 
              screen_params = screen_params, 
              auto_initial = 2,
              nodesize = 1, 
              debug = 1, 
              verbose = 1, 
              initial = TRUE, 
              num_processors = 1, 
              min_features = 10)

end <- proc.time()
end - start 

