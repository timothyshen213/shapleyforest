library(devtools)
library(roxygen2)

setwd("C:/Users/timot/OneDrive/Documents/1 - UCLA/Research/shapff")
detach("package:shapff", unload = TRUE)
remove.packages("shapff")

library(randomForest)
library(fastshap)
library(shapviz)
library(dplyr)
library(purrr)
library(reshape2)
library(ggplot2)

devtools::install_github("timothyshen213/shapff")

library(WGCNA)
library(fuzzyforest)
library(mvtnorm)
library(treeshap)
library(shapff)
options(warn = 1)
shap_ff <- function(X, y, Z=NULL, shap_model = 1, shap_type = "shapley", module_membership,
                   screen_params = fuzzyforest:::screen_control(min_ntree=5000),
                   select_params = fuzzyforest:::select_control(min_ntree=5000),
                   final_ntree = 5000,
                   num_processors = 1, parallel = 1, nodesize, 
                   test_features=NULL, test_y=NULL, nsim = 1, 
                   final_nsim = 100) {
  if(!is.null(Z)) {
    if (!is.data.frame(Z)) {
      stop("Z must be a data.frame.",
           call. = FALSE)
    }
  }
  if (!(is.vector(y) || is.factor(y))) {
    stop("y must be vector or factor")
  }
  if(!is.data.frame(X)) {
    stop("X must be a data.frame.", call. = FALSE)
  }
  CLASSIFICATION <- is.factor(y)
  if(CLASSIFICATION == TRUE) {
    if(missing(nodesize)){
      nodesize <- 1
    }
  }
  if(CLASSIFICATION == FALSE) {
    if(missing(nodesize)){
      nodesize <- 5
    }
  }
  
  screen_control <- screen_params
  select_control <-  select_params
  module_list <- unique(module_membership)
  
  #parallelization
  if(parallel == 1 && num_processors > 1) {
    cl = parallel::makeCluster(num_processors)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }
  if(parallel == 2 && num_processors > 1) {
    cl <- snow::makeCluster(num_processors, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    on.exit(snow::stopCluster(cl), add = TRUE)
  }
  survivors <- vector('list', length(module_list))
  drop_fraction <- screen_control$drop_fraction
  mtry_factor <- screen_control$mtry_factor
  ntree_factor <- screen_control$ntree_factor
  min_ntree <- screen_control$min_ntree
  keep_fraction <- screen_control$keep_fraction
  if(ncol(X)*keep_fraction < select_control$number_selected){
    warning(c("ncol(X)*keep_fraction < number_selected", "\n",
              "number_selected will be set to floor(ncol(X)*keep_fraction)"))
    select_control$number_selected <- max(floor(ncol(X)*keep_fraction), 1)
  }
  
  total_iterations <- length(module_list)
  progress_bar <- txtProgressBar(min = 0, max = total_iterations, style = 3)
  
  for (i in 1:length(module_list)) {
    setTxtProgressBar(progress_bar, i)
    module <- X[, which(module_membership == module_list[i]), drop=FALSE]
    num_features <- ncol(module)
    #TUNING PARAMETER mtry_factor
    if(CLASSIFICATION == TRUE) {
      mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
      if(missing(nodesize)){
        nodesize <- 1
      }
    }
    if(CLASSIFICATION == FALSE) {
      mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
      if(missing(nodesize)){
        nodesize <- 5
      }
    }
    #TUNING PARAMETER ntree_factor
    ntree <- max(num_features*ntree_factor, min_ntree)
    #TUNING PARAMETER keep_fraction
    target = ceiling(num_features * keep_fraction)
    
    while (num_features >= target){
      if(num_processors > 1) {
        rf = `%dopar%`(foreach(ntree = rep(ntree/num_processors, num_processors)
                               , .combine = randomforest::combine, .packages = 'randomForest'),
                       #second argument to '%dopar%'
                       randomForest(module , y, ntree = ntree, mtry = mtry,
                                    importance = TRUE, scale = FALSE, nodesize=nodesize))
        
      }
      if(num_processors == 1) {
        rf <- randomForest(module, y, ntree = ntree, mtry = mtry,
                           importance = TRUE, scale = FALSE,
                           nodesize = nodesize)
      }
      
      #Feature selection via fastshap
      if (shap_model == 1){
        shap <- suppressMessages(fastshap::explain(rf, X = module, nsim = nsim, pred_wrapper = predict))
        var_importance <- colMeans(abs(shap))
        var_importance <- sort(var_importance, decreasing = TRUE)
        var_importance <- data.frame(Feature = var_importance)
      }
      #Feature selection via permutation vims
      if (shap_model == 0){
        var_importance <- importance(rf, type=1, scale=FALSE)
        var_importance <- var_importance[order(var_importance[, 1],
                                               decreasing=TRUE), ,drop=FALSE]
      }
      reduction <- ceiling(num_features*drop_fraction)
      
      if(num_features - reduction > target) {
        trimmed_varlist <- var_importance[1:(num_features - reduction), ,drop=FALSE]
        features <- row.names(trimmed_varlist)
        module <- module[, which(names(module) %in% features)]
        num_features <- length(features)
        if(CLASSIFICATION == TRUE) {
          mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
        }
        if(CLASSIFICATION == FALSE) {
          mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
        }
        ntree <- max(num_features*ntree_factor, min_ntree)
      }
      else {
        num_features <- target - 1
        mod_varlist <- var_importance[, 1][1:target]
        features <- row.names(var_importance)[1:target]
        survivors[[i]] <- cbind(features, mod_varlist)
        row.names(survivors[[i]]) <- NULL
        survivors[[i]] <- as.data.frame(survivors[[i]])
        survivors[[i]][, 1] <- as.character(survivors[[i]][, 1])
        survivors[[i]][, 2] <- as.numeric(as.character(survivors[[i]][, 2]))
      }
    }
  }
  
  cat("\nSelection Step ...")
  
  survivor_list <- survivors
  names(survivor_list) <- module_list
  survivors <- do.call('rbind', survivors)
  survivors <- as.data.frame(survivors, stringsAsFactors = FALSE)
  survivors[, 2] <- as.numeric(survivors[, 2])
  names(survivors) <- c("featureID", "Permutation VIM")
  X_surv <- X[, names(X) %in% survivors[,1]]
  if(!is.null(Z)) {
    X_surv <- cbind(X_surv, Z, stringsAsFactors=FALSE)
  }
  select_args <- list(X_surv, y, num_processors, nodesize)
  select_args <- c(select_args, select_control)
  names(select_args)[1:4] <- c("X", "y", "num_processors", "nodesize")
  
  #Selection step via SHAP values
  if (shap_model == 1){
    select_results <- shapselect_RF(select_args$X, select_args$y, select_args$drop_fraction, 
                                    select_args$number_selected, select_args$mtry_factors,
                                    select_args$ntree_factor, select_args$min_ntree,
                                    select_args$num_processors, select_args$nodesize, select_args$cl,
                                    nsim = nsim)
  }
  #Selection step via permutation vims
  if (shap_model == 0){
    select_results <- do.call(fuzzyforest::select_RF, select_args)
  }
  
  final_list <- select_results[[1]]
  selection_list <- select_results[[2]]
  final_list[, 2] <- round(as.numeric(final_list[, 2]), 4)
  row.names(final_list) <- NULL
  colnames(final_list) <- c("feature_name", "variable_importance")
  final_list <- as.data.frame(final_list, stringsAsFactors=FALSE)
  final_list[, 2] <- as.numeric(final_list[, 2])
  final_list <- cbind(final_list, rep(".", dim(final_list)[1]),
                      stringsAsFactors=FALSE)
  names(final_list)[3] <- c("module_membership")
  select_X <- names(X)[which(names(X) %in% final_list[, 1])]
  select_mods <- module_membership[which(names(X) %in% final_list[,1])]
  select_order <- final_list[, 1][which(final_list[,1] %in% names(X))]
  select_mods <- select_mods[match(select_order, select_X)]
  final_list[, 3][final_list[, 1] %in% names(X)] <- select_mods
  final_X <- X[, names(X) %in% final_list[, 1], drop=FALSE]
  if(!is.null(Z)) {
    final_X <- cbind(final_X, Z[, names(Z) %in% final_list[, 1], drop=FALSE],
                     stringsAsFactors=FALSE)
  }
  current_p <- dim(final_X)[2]
  if(CLASSIFICATION == TRUE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*current_p/3),
                      current_p)
  }
  if(CLASSIFICATION == FALSE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*current_p),
                      current_p)
  }
  if(!is.null(test_features)) {
    test_features <- test_features[, which(names(test_features) %in%
                                             names(final_X))]
  }
  final_rf <- randomForest(x=final_X, y=y, mtry=final_mtry, ntree=final_ntree,
                           importance=TRUE, nodesize=nodesize,
                           xtest=test_features, ytest=test_y)
  final_module_membership <- as.data.frame(cbind(names(X), module_membership),
                                           stringsAsFactors=FALSE)
  names(final_module_membership) <- c("feature_name", "module")
  
  #Final SHAP value calculation
  if (shap_type == "shapley"){
    shap_final_obj <- fastshap::explain(final_rf, X = final_X, nsim = final_nsim, 
                                        pred_wrapper = predict, shap_only = FALSE)
    shap_final <- shap_final_obj$shapley_values
    var_importance_final <- colMeans(abs(shap_final))
    var_importance_final <- sort(var_importance_final, decreasing = TRUE)
    var_importance_final <- data.frame(vim = var_importance_final)
  }
  if (shap_type == "tree"){
    unified_model <- randomForest.unify(final_rf, final_X)
    shap_final_obj <- treeshap::treeshap(unified_model, x = final_X, 
                                         interactions = TRUE, verbose = FALSE)
    shap_final <- shap_final_obj$shaps
    var_importance_final <- colMeans(abs(shap_final))
    var_importance_final <- sort(var_importance_final, decreasing = TRUE)
    var_importance_final <- data.frame(vim = var_importance_final)
  }
  shap_final_list <- data.frame(feature_name = rownames(var_importance_final),
                                variable_importance = var_importance_final[,1])
  shap_final_list[, 2] <- round(as.numeric(shap_final_list[, 2]), 4)
  row.names(shap_final_list) <- NULL
  shap_final_list <- as.data.frame(shap_final_list, stringsAsFactors=FALSE)
  shap_final_list[, 2] <- as.numeric(shap_final_list[, 2])
  shap_final_list <- cbind(shap_final_list, rep(".", dim(shap_final_list)[1]),
                           stringsAsFactors=FALSE)
  names(shap_final_list)[3] <- c("module_membership")
  
  select_X <- names(X)[which(names(X) %in% shap_final_list[, 1])]
  select_mods <- module_membership[which(names(X) %in% shap_final_list[,1])]
  select_order <- shap_final_list[, 1][which(shap_final_list[,1] %in% names(X))]
  select_mods <- select_mods[match(select_order, select_X)]
  shap_final_list[, 3][shap_final_list[, 1] %in% names(final_X)] <- select_mods
  
  out <- shap_fuzzy_forest(final_rf, final_module_membership,
                           survivor_list=survivor_list,
                           selection_list=selection_list, 
                           final_shap = shap_final_list,
                           shap_obj = shap_final_obj,
                           final_X = final_X)
  
  cat("Done \n")
  
  return(out)
}

shap_wff <- function(X, y, Z=NULL, shap_model = 1, shap_type = "tree", WGCNA_params=WGCNA_control(p=6),
                    screen_params=fuzzyforest:::screen_control(min_ntree=5000),
                    select_params=fuzzyforest:::select_control(min_ntree=5000),
                    final_ntree=500, num_processors, parallel=1, nodesize,
                    test_features=NULL, test_y=NULL, nsim=1, final_nsim=100) {
  #browser()
  if ( !("package:WGCNA" %in% search()) ) {
    stop("WGCNA must be loaded and attached. Type library(WGCNA) to do so.",
         call. = FALSE)
  }
  if (!(is.vector(y) || is.factor(y))) {
    stop("y must be vector or factor")
  }
  #WGCNA yields errors if X is an integer, here we convert integers to numeric.
  integer_test <- sapply(X, is.integer)
  if( sum(integer_test) > 0 ) {
    ints <- which(integer_test == TRUE)
    X[, ints] <- apply(X[, ints, drop=FALSE], 2, as.numeric)
  }
  numeric_test <- sapply(X, is.numeric)
  if (sum(numeric_test) != dim(X)[2]) {
    stop("The columns of X must be numeric.")
  }
  CLASSIFICATION <- is.factor(y)
  if(CLASSIFICATION == TRUE) {
    if(missing(nodesize)){
      nodesize <- 1
    }
  }
  if(CLASSIFICATION == FALSE) {
    if(missing(nodesize)){
      nodesize <- 5
    }
  }
  if (shap_model != 0 && shap_model != 1) {
    stop("shap must be 0 or 1. Type help(shapff) or help(shapwff) for details.")
  }
  
  WGCNA_control <- WGCNA_params
  screen_control <- screen_params
  select_control <-  select_params
  WGCNA_args <- list(X,WGCNA_control$power)
  WGCNA_args <- c(WGCNA_args, WGCNA_control$extra_args)
  names(WGCNA_args) <- c("datExpr", "power", names(WGCNA_control$extra_args))
  bwise <- do.call("blockwiseModules", WGCNA_args)
  module_membership <- bwise$colors
  screen_drop_fraction <- screen_control$drop_fraction
  screen_keep_fraction <- screen_control$keep_fraction
  screen_mtry_factor <- screen_control$mtry_factor
  screen_ntree_factor <- screen_control$ntree_factor
  screen_min_ntree <- screen_control$min_ntree
  cat("Screening Step ... \n")
  out <- shap_ff(X, y, Z, shap_model, shap_type, module_membership,
                screen_control, select_control, final_ntree,
                num_processors, nodesize=nodesize,
                test_features=test_features, test_y=test_y)
  out$WGCNA_object <- bwise
  return(out)
}


id <- 1005
set.seed(1)

rep_num <- 100
keep_frac <- c(0.01, 0.05, 0.1, 0.15, 0.25)
drop_frac <- c(0.05, 0.1, 0.25, 0.5)
mtry_factor <- c(0.5, 1, 2)
p <- c(100,1000)
n <- c(100)

param_list <- list(keep_frac, drop_frac, mtry_factor, p, n)
param_settings <- expand.grid(param_list)
param_settings <- param_settings[, 5:1]
names(param_settings) <- c("n", "p", "mtry_factor", "drop_fraction", "keep_fraction")

current_sim_params <- param_settings[ceiling((id)/rep_num), ]

## LINEAR TEST ##
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


all_modules <- lapply(1:number_of_mods, function(j) sim_mod(n, p_per_group, corr))
all_modules[[number_of_groups]] <- matrix(rnorm(p_per_group * n), nrow = n, ncol = p_per_group)
X <- do.call(cbind, all_modules)
beta <- rep(0, p_per_group * (number_of_mods + 1))
beta[vim_list] <- beta_list
y <- X %*% beta + rnorm(n, sd = 0.1)
X <- as.data.frame(X)
names(X) <- paste("V", 1:p, sep = "")
mtry_factor <- 1
screen_params <- screen_control(drop_fraction = drop_fraction, keep_fraction = keep_fraction, 
                                    mtry_factor = mtry_factor)
select_params <- select_control(number_selected = 10, drop_fraction = drop_fraction, 
                                    mtry_factor = mtry_factor)
y <- as.numeric(y)
ff <- shapwff(X, y, shap_model = "full", screen_params = screen_params, select_params = select_params, 
                   num_processors = 1, nodesize = 1, min_features = 30, verbose = 1, debug = 2, final_nsim = 1)

# via treeshap
plot_importance(ff)
plot_interactions(ff, kind = "beeswarm", max_display = 10) # plots interactions, diag = main effect
plot_interactions(ff, kind = "matrix", max_display = 10) # interaction matrix
plot_dependence(ff, features = "V1", color_var = names(ff$final_X), interaction = FALSE)

plot_decisions(ff)
