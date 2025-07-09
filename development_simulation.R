.libPaths("/u/home/t/tzshen/R")

library(WGCNA)
library(randomForest)
library(mvtnorm)
library(fuzzyforest)

set.seed(2002)

keep_frac <- c(0.1)
drop_frac <- c(0.25)
mtry_factor <- c(1)
p <- c(1000)
n <- c(100)

current_sim_params <- c(n, p, mtry_factor, drop_frac, keep_frac)
names(current_sim_params) <- c("n", "p", "mtry_factor", "drop_fraction", "keep_fraction")

sim_mod <- function(n, p, corr) {
  sigma <- matrix(corr, nrow = p, ncol = p)
  diag(sigma) <- 1
  X <- rmvnorm(n, sigma = sigma)
  return(X)
}

n <- as.numeric(current_sim_params[1])
p <- as.numeric(current_sim_params[2])
mtry_factor <- as.numeric(current_sim_params[3])
keep_fraction <- as.numeric(current_sim_params[4])
drop_fraction <- as.numeric(current_sim_params[5])

corr <- 0.8

if (p == 1000) {
  number_of_groups <- 10
  number_of_mods <- number_of_groups - 1
  p_per_group <- p/number_of_groups
  vim_list <- c(1:3, 901:903)
  vim_interest <- c(1:1000)
  beta_list <- rep(c(5, 5, 2), 2)
}

# simulation
all_modules <- lapply(1:number_of_mods, function(j) sim_mod(n, p_per_group, corr))
all_modules[[number_of_groups]] <- matrix(rnorm(p_per_group * n), nrow = n, ncol = p_per_group)
X <- do.call(cbind, all_modules)
beta <- rep(0, p_per_group * (number_of_mods + 1))
beta[vim_list] <- beta_list
y <- X %*% beta + rnorm(n, sd = 0.1)
X <- as.data.frame(X)
names(X) <- paste("V", 1:p, sep = "")

# wgcna
powers <- c(1:20)
sft <- pickSoftThreshold(X, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6

net <- blockwiseModules(X,
                        power = softPower,
                        TOMType = "signed",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)


moduleLabels <- net$colors
table(moduleLabels)
module_membership <- factor(moduleLabels)

# shapleyforest (via ranger & parallel)
runtime <- system.time({
  mtry_factor <- 1
  screen_params <- screen_control(drop_fraction = drop_fraction, keep_fraction = keep_fraction,
                                  mtry_factor = mtry_factor)
  select_params <- select_control(number_selected = 10, drop_fraction = drop_fraction,
                                  mtry_factor = mtry_factor)
  y <- as.numeric(y)
  sff <- shapff(X, y, module_membership = moduleLabels,
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
})

