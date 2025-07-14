.libPaths("/u/home/t/tzshen/R")

id <- as.integer(Sys.getenv("SGE_TASK_ID"))
id <- id
if (is.na(id)) id <- 1 # if run manually only do for id = 1

library(WGCNA)
#library(randomForest)
library(mvtnorm)
library(fuzzyforest)
library(ranger)
library(future)
library(future.apply)
library(doFuture)
library(tidyverse)

source(here("./test/test.R"))

set.seed(id)

rep_num <- 100
keep_frac <- c(0.01, 0.05, 0.1, 0.15, 0.25)
drop_frac <- c(0.05, 0.1, 0.25, 0.5)
mtry_factor <- c(0.5, 1, 2)
p <- c(100, 1000)
n <- c(100)

param_list <- list(keep_frac, drop_frac, mtry_factor, p, n)
param_settings <- expand.grid(param_list)
param_settings <- param_settings[, 5:1]
names(param_settings) <- c("n", "p", "mtry_factor", "drop_fraction", "keep_fraction")

param_settings
current_sim_params <- param_settings[ceiling((id)/rep_num), ]

sim_number <- 10
sim_results <- list()
sim_results_1 <- list()
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
if (p == 100) {
  number_of_groups <- 4
  number_of_mods <- number_of_groups - 1
  p_per_group <- p/number_of_groups
  vim_list <- c(1:3, 76:78)
  vim_interest <- c(1:4, 76:79)
  beta_list <- rep(c(5, 5, 2), 2)
}
if (p == 1000) {
  number_of_groups <- 10
  number_of_mods <- number_of_groups - 1
  p_per_group <- p/number_of_groups
  vim_list <- c(1:3, 901:903)
  vim_interest <- c(1:4, 901:904)
  beta_list <- rep(c(5, 5, 2), 2)
}

registerDoSEQ()
set.seed(id)

all_modules <- lapply(1:number_of_mods, function(j) sim_mod(n, p_per_group, corr))
all_modules[[number_of_groups]] <- matrix(rnorm(p_per_group * n), nrow = n, ncol = p_per_group)
X <- do.call(cbind, all_modules)
beta <- rep(0, p_per_group * (number_of_mods + 1))
beta[vim_list] <- beta_list
y <- X %*% beta + rnorm(n, sd = 0.1)
X <- as.data.frame(X)
names(X) <- paste("V", 1:p, sep = "")
#mtry_factor <- 1
screen_params <- screen_control(drop_fraction = drop_fraction, keep_fraction = keep_fraction, 
                                mtry_factor = mtry_factor)
select_params <- select_control(number_selected = 10, drop_fraction = drop_fraction, 
                                mtry_factor = mtry_factor)
y <- as.numeric(y)
powers <- c(1:20)
sft <- pickSoftThreshold(X, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 3
net <- blockwiseModules(X,
                        power = softPower,
                        TOMType = "signed",
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)
moduleLabels <- net$colors
module_membership <- factor(moduleLabels)

for (l in 1:sim_number) {
  ff_1 <- shapff(X, y, module_membership = module_membership,
               select_params = select_params,
               shap_model = "full",
               screen_params = screen_params,
               auto_initial = 4,
               nodesize = 1,
               debug = 1,
               verbose = 1,
               initial = TRUE,
               num_processors = 1,
               min_features = 10,
               seed = id)
  shap_feature_1 <- data.frame(index = ff_1$final_SHAP[1], 
                             ff_1$final_SHAP[2], row.names = NULL)
  vim_1 <- shap_feature_1
  sim_results_1[[l]] <- vim_1
  
  ff <- shapff(X, y, module_membership = module_membership,
                        select_params = select_params,
                        shap_model = "full",
                        screen_params = screen_params,
                        auto_initial = 4,
                        nodesize = 1,
                        debug = 1,
                        verbose = 1,
                        initial = TRUE,
                        num_processors = 8,
                        min_features = 10,
                        seed = id)
  
  shap_feature <- data.frame(index = ff$final_SHAP[1], 
                             ff$final_SHAP[2], row.names = NULL)
  vim <- shap_feature
  sim_results[[l]] <- vim
}

summarize_sim <- function(sim_results, vim_interest) {
  num_interesting_feat <- length(vim_interest)
  vims <- matrix(0, nrow = num_interesting_feat, ncol = length(sim_results))
  selected <- matrix(0, nrow = num_interesting_feat, ncol = length(sim_results))
  
  for (r in 1:num_interesting_feat) {
    current_feature <- paste("V", vim_interest[r], sep = "")
    for (t in 1:length(sim_results)) {
      current_vims <- sim_results[[t]]
      if (current_feature %in% current_vims[, 1]) {
        current_rank <- which(current_vims[, 1] == current_feature)
        vims[r, t] <- current_vims[current_rank, 2]
        selected[r, t] <- 1
      } else {
        vims[r, t] <- 0
      }
    }
  }
  
  mean_vims <- apply(vims, 1, mean)
  selected_props <- apply(selected, 1, mean)
  
  sim_out <- data.frame(
    feature = paste0("V", vim_interest),
    mean_vim = mean_vims,
    selected_prop = selected_props,
    row.names = NULL
  )
  
  return(sim_out)
}



out <- list(summarize_sim(sim_results, vim_interest), 
            vim_interest, 
            summarize_sim(sim_results_1, vim_interest), 
            current_sim_params, 
            runtime = list(ff$runtimes, ff_1$runtimes))

df_single <- out[[1]] %>%
  filter(gsub("V", "", feature) %in% as.character(out[[2]])) %>%
  mutate(feature_num = as.numeric(gsub("V", "", feature)),
         Label = "Single")

df_par <- out[[3]] %>%
  filter(gsub("V", "", feature) %in% as.character(out[[2]])) %>%
  mutate(feature_num = as.numeric(gsub("V", "", feature)),
         Label = "Parallel")

# Combine
plot_data <- bind_rows(df_single, df_par) %>%
  arrange(Label, feature_num) %>%
  mutate(feature = factor(as.character(feature, levels = unique(feature_num))))

# Plot
ggplot(plot_data, aes(x = feature, y = selected_prop, color = Label, group = Label)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(title = "Selected Proportion Comparison: Single vs Non Parallel",
       x = "Feature",
       y = "Selected Proportion",
       color = "Run Type") +
  theme_minimal()



dir.create("out", showWarnings = FALSE)
out_filename <- paste("parallel/out", id, sep = "")
save(out, file = out_filename)