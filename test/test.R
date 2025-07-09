.libPaths("/u/home/t/tzshen/R")

id <- as.integer(Sys.getenv("SGE_TASK_ID"))
id <- id
if (is.na(id)) id <- 1 # if run manually only do for id = 1

num_cores <- 8

library(WGCNA)
#library(randomForest)
library(mvtnorm)
library(fuzzyforest)
library(ranger)
library(future)
library(future.apply)
library(doFuture)
library(tidyverse)

shapff <- function(X, y, Z=NULL, shap_model = "full", module_membership,
                   min_features = 20, verbose = 1, debug = 2, 
                   initial = TRUE, auto_initial = NULL, 
                   screen_params = fuzzyforest:::screen_control(min_ntree=5000),
                   select_params = fuzzyforest:::select_control(min_ntree=5000),
                   final_ntree = 5000,
                   num_processors = n_cores, nodesize, 
                   test_features=NULL, test_y=NULL, nsim = 1, 
                   final_nsim = 100, seed = 1234) {
  
  options(doFuture.rng.onMisuse = "ignore")
  set.seed(seed)
  
  ## validating prerequisites
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
  if (verbose == 0){
    options(warn = -1)
  }
  if (num_processors > 1){
    parallel = TRUE
  } else {
    parallel = FALSE
  }
  
  # sets parameters for each step
  screen_control <- screen_params
  select_control <- select_params
  
  # obtains module membership (ie. from WGCNA)
  module_list <- unique(module_membership)
  
  # adjusts keep_fraction if below minimum threshold
  if(ncol(X)*screen_control$keep_fraction < select_control$number_selected){
    if (verbose != 0){
      warning(c("\n\n ncol(X)*keep_fraction < number_selected", "\n",
                "number_selected will be set to floor(ncol(X)*keep_fraction)"))
    }
    select_control$number_selected <- max(floor(ncol(X)*keep_fraction), 1)
  }
  
  ## Screening Step
  screen_result <- shapscreen_RF(
    X = X,
    y = y,
    module_list = module_list,
    module_membership = module_membership,
    screen_control = screen_control,
    select_control = select_control,
    shap_model = shap_model,
    CLASSIFICATION = CLASSIFICATION,
    min_features = min_features,
    nsim = nsim,
    nodesize = nodesize,
    num_processors = num_processors,
    verbose = verbose,
    debug = debug,
    parallel = parallel,
    seed = seed
  )
  
  assign("screen_result", screen_result, envir = .GlobalEnv)
  
  survivors <- screen_result$survivor_results
  initial_screen <- screen_result$initial_screen
  
  ## Initial Screening Output Procedure
  
  if (!is.null(auto_initial)){
    if (auto_initial %in% c("1", "2")) { # save output
      write.csv(initial_screen, "initial_screen.csv", row.names = FALSE)
      assign("initial_screen", initial_screen, envir = .GlobalEnv)
      if (verbose != 0){
        cat("\n\n Dataframe saved as 'initial_screen.csv'.\n")
      }
      if (auto_initial == "1"){ # stops running
        options(warn = 1)
        stop("\n\n Execution stopped as per user choice.\n")
      }
    } else { # skips all
      if (auto_initial == "3"){
        options(warn = 1)
        stop("\n\n Execution stopped as per user choice.\n")
      }
    }
    initial = FALSE
  }
  
  if (initial == TRUE){
    cat("\n\nDisplaying Initial Screen (First Step) Variable Importance ... ")
    print(knitr::kable(initial_screen))
    
    # User given an option to save initial screening
    save_prompt <- readline(prompt = "Choose an option on output (1, 2, 3, or 4): 
    1. Save and Stop 
    2. Save and Continue 
    3. Don't Save and Stop
    4. Don't Save and Continue")  
    
    # if encounters invalid response to prompt, asks again
    while (!save_prompt %in% c("1", "2", "3", "4")) {
      save_prompt <- readline(prompt = "Invalid choice. Please enter 1, 2, 3 or 4: ")
    }
    
    # follows the user defined save procedure
    if (tolower(save_prompt) %in% c("1", "2")) { # save output
      write.csv(initial_screen, "initial_screen.csv", row.names = FALSE)
      assign("initial_screen", initial_screen, envir = .GlobalEnv)
      cat("\n\nDataframe saved as 'initial_screen.csv'.\n")
      if (tolower(save_prompt) == "1"){ # stops running
        options(warn = 1)
        stop("\n\nExecution stopped as per user choice.\n")
      }
    } else { # skips all
      if (tolower(save_prompt) == "3"){
        options(warn = 1)
        stop("\n\nExecution stopped as per user choice.\n")
      }
    }
  }
  
  
  # verbose UI
  if (verbose != 0){cat("\nSelection Step ...")}
  
  (survivor_list <- survivors)
  (names(survivor_list) <- module_list)
  survivors <- do.call('rbind', survivors)
  survivors <- as.data.frame(survivors, stringsAsFactors = FALSE)
  survivors[, 2] <- as.numeric(survivors[, 2])
  names(survivors) <- c("featureID", "Permutation VIM")
  X_surv <- X[, names(X) %in% survivors[,1]]
  if(!is.null(Z)) {
    X_surv <- cbind(X_surv, Z, stringsAsFactors=FALSE)
  }
  
  ## Selection Step
  
  # sets up selection step parameters
  select_args <- list(X_surv, y, num_processors, nodesize)
  select_args <- c(select_args, select_control)
  names(select_args)[1:4] <- c("X", "y", "num_processors", "nodesize")
  
  # RFE via fastshap
  if (shap_model == "full"){
    select_results <- shapselect_RF(select_args$X, select_args$y, select_args$drop_fraction, 
                                    select_args$number_selected, CLASSIFICATION,
                                    select_args$mtry_factors,
                                    select_args$ntree_factor, select_args$min_ntree,
                                    select_args$num_processors, select_args$nodesize, select_args$cl,
                                    nsim = nsim)
  }
  
  # RFE via permutation VIMs
  if (shap_model == "after"){
    select_results <- do.call(fuzzyforest::select_RF, select_args)
  }
  # gathers final surviving list
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
  
  # selection step results for X, module membership
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
  
  # sets Random Forest Variable for final run
  current_p <- dim(final_X)[2]
  if(CLASSIFICATION == TRUE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*current_p/3),
                      current_p)
  }
  if(CLASSIFICATION == FALSE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*current_p),
                      current_p)
  }
  
  # gets test predictors from final_X for Random Forest
  if(!is.null(test_features)) {
    test_features <- test_features[, which(names(test_features) %in%
                                             names(final_X))]
  }
  
  plan(sequential)
  
  # final Random Forest
  final_rf <- ranger(
    dependent.variable.name = "y",
    data = data.frame(y = y, final_X),
    mtry = final_mtry,
    num.trees = final_ntree,
    importance = "permutation",
    min.node.size = nodesize,
    probability = CLASSIFICATION,        
    classification = TRUE,
    verbose = FALSE,
    num.threads = 1
  )
  
  # extracts module membership of final survivors
  final_module_membership <- as.data.frame(cbind(names(X), module_membership),
                                           stringsAsFactors=FALSE)
  names(final_module_membership) <- c("feature_name", "module")
  
  # final variable importance measure via fastshap
  predict_function <- if (CLASSIFICATION) {
    if (num_classes == 2) {
      function(object, newdata) {
        preds <- predict(object, data = newdata)$predictions
        return(preds[, 1])  # ranger returns a matrix with one column: prob of class 1
      }
    } else if (num_classes > 2) {
      function(object, newdata) {
        preds <- predict(object, data = newdata)$predictions
        return(preds)  # returns matrix with class probabilities
      }
    } else {
      stop("Invalid or single-class data in y")
    }
  } else {
    function(object, newdata) {
      preds <- predict(object, data = newdata)$predictions
      return(preds)  # regression: numeric vector
    }
  }
  
  plan(multisession, workers = num_processors)
  
  shap_final_obj <- suppressMessages(fastshap::explain(
    object = final_rf,
    X = final_X,
    pred_wrapper = predict_function,
    nsim = final_nsim,
    adjust = TRUE,
    parallel = TRUE,
    shap_only = FALSE, 
    .packages = "ranger"
  ))
  
  plan(sequential)
  
  shap_final <- shap_final_obj$shapley_values
  var_importance_final <- colMeans(abs(shap_final))
  var_importance_final <- sort(var_importance_final, decreasing = TRUE)
  var_importance_final <- data.frame(vim = var_importance_final)
  
  # extracts final shapley values for each survivor
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
  # output function
  out <- shapley_forest(final_rf, final_module_membership,
                        survivor_list=survivor_list,
                        selection_list=selection_list, 
                        final_shap = shap_final_list,
                        shap_obj = shap_final_obj,
                        final_X = final_X)
  
  # verbose message
  if (verbose != 0){cat("Done \n")}
  
  return(out)
}

shapscreen_RF <- function(X, y, module_list, module_membership, screen_control, 
                          select_control, shap_model = "full", 
                          CLASSIFICATION, 
                          min_features, nsim, nodesize, 
                          num_processors, verbose, 
                          debug, parallel, 
                          seed) {
  ## Screening Step
  survivors <- vector('list', length(module_list)) # initiates survivor list
  
  # parameters for screening
  drop_fraction <- screen_control$drop_fraction 
  mtry_factor <- screen_control$mtry_factor
  ntree_factor <- screen_control$ntree_factor
  min_ntree <- screen_control$min_ntree
  keep_fraction <- screen_control$keep_fraction
  
  # begins RFE
  total_iterations <- length(module_list)
  
  # sets UI Progress Bar
  if (verbose !=0){
    cat("\nScreening Step ...")
    if (!parallel){
      progress_bar <- txtProgressBar(min = 0, max = total_iterations, style = 3)
    }
    else {
      cat("\nParallelizing with", num_processors, "processors\n")
    }
  }
  
  # initiates initial_screen module list
  initial_screen <- data.frame("Module Color" = character(0), "Feature Name" = character(0),
                               "VIM" = numeric(0), "Survivor" = logical(0))
  
  # Recursive Feature Elimination
  if (!parallel){
    # runs RFE for given module
    while (num_features >= target){
      
      # debug warning message
      if (debug != -1 || debug != 1){
        # if module has low features, only keep non zero important features
        if (num_features <= min_features){
          if (verbose != 0){
            warning(sprintf("\n\n Module %s has fewer than %d features! 
                            All non-zero important features will be kept during screening.", 
                            module_list[i], min_features))
          }
          keep = TRUE
        } else{
          keep = FALSE
        }
      }
      
      # runs Random Forest
      rf <- randomForest(module, y, ntree = ntree, mtry = mtry,
                         importance = TRUE, scale = FALSE,
                         nodesize = nodesize)
      
      ## run Feature Selection
      
      # via fastshap
      if (shap_model == "full"){
        # classification case
        if (CLASSIFICATION == TRUE){
          num_classes <- nlevels(y)
          if (num_classes == 2){
            prediction <- function(object, newdata) {
              prob <- predict(object, newdata = newdata, type = "prob")
              return(prob[,2])
            }
          }
          else if (num_classes > 2) {
            prediction <- function(object, newdata) {
              prob <- predict(object, newdata = newdata, type = "prob")
              return(prob)
            }
          } else {
            options(warn = 1)
            stop("Invalid or single-class data in y")
          }
          shap <- suppressMessages(
            fastshap::explain(rf, 
                              X = module, 
                              nsim = nsim, 
                              pred_wrapper = prediction)
          )
        }
        
        # regression case
        if (CLASSIFICATION == FALSE){
          shap <- suppressMessages(
            fastshap::explain(rf, 
                              X = module, 
                              nsim = nsim, 
                              pred_wrapper = predict)
          )
        }
        var_importance <- colMeans(abs(shap))
        var_importance <- sort(var_importance, decreasing = TRUE)
        var_importance <- data.frame(Feature = var_importance)
      }
      
      # via permutation VIMs
      if (shap_model == "after"){
        var_importance <- importance(rf, type=1, scale=FALSE)
        var_importance <- var_importance[order(var_importance[, 1],
                                               decreasing=TRUE), ,drop=FALSE]
      }
      
      # sets reduction value for feature elimination
      reduction <- ceiling(num_features*drop_fraction)
      
      # debug feature
      if (debug != -1 || debug != 1){
        if(keep == TRUE){
          # keeps only non-zero important values
          trimmed_varlist <- var_importance[var_importance > 0, , drop = FALSE]
          features <- row.names(trimmed_varlist)
          module <- module[, which(names(module) %in% features)]
          target = num_features - reduction
          num_features <- length(features)
        }
      }
      
      # if target not reached, continue RFE
      if(num_features - reduction > target) {
        # keep top variable important features
        trimmed_varlist <- var_importance[1:(num_features - reduction), ,
                                          drop=FALSE]
        features <- row.names(trimmed_varlist)
        module <- module[, which(names(module) %in% features)]
        num_features <- length(features)
        
        # adjust Random Forest variables for next iteration
        if(CLASSIFICATION == TRUE) {
          mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
        }
        if(CLASSIFICATION == FALSE) {
          mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
        }
        ntree <- max(num_features*ntree_factor, min_ntree)
      }
      
      # if target is reached, stop RFE
      else {
        num_features <- target - 1 # stops while loop
        mod_varlist <- var_importance[, 1][1:target] # module membership for survivors
        features <- row.names(var_importance)[1:target] # feature names for survivors
        
        # extracts survivor list for given module
        survivors[[i]] <- cbind(features, mod_varlist)
        row.names(survivors[[i]]) <- NULL
        survivors[[i]] <- as.data.frame(survivors[[i]])
        survivors[[i]][, 1] <- as.character(survivors[[i]][, 1])
        survivors[[i]][, 2] <- as.numeric(as.character(survivors[[i]][, 2]))
      }
      
      # prints initial screening output for each module
      if (is_initial == TRUE){
        n <- nrow(var_importance)
        mod_name <- rep(module_name, times = n)
        feature_name <- rownames(var_importance)
        trimmed_features <- rownames(trimmed_varlist)
        survived <- feature_name %in% trimmed_features
        temp <- data.frame("Module Color" = mod_name,
                           "Feature Name" = feature_name, 
                           "VIM" = var_importance[,1], 
                           "Survivor" = survived)
        initial_screen <- rbind(initial_screen, temp)
        is_initial = FALSE
      }
    }
    survivor_results <- survivors
  }
  else {
    # progress bar
    # if (verbose !=0){
    #   handlers(global = TRUE)
    #   handlers("progress") 
    # }
    
    # starts parallelization
    registerDoFuture() 
    plan(multisession, workers = num_processors)
    
    # parallelizes RFE at each WGCNA module
    results <- future_lapply(seq_along(module_list), function(i) {
      
      # gets specific WGCNA module
      module_name <- module_list[i]
      module <- X[, which(module_membership == module_name), drop = FALSE]
      num_features <- ncol(module)
      
      # classification or regression set up for parameters
      CLASSIFICATION <- is.factor(y)
      mtry <- if (CLASSIFICATION) min(ceiling(mtry_factor * num_features / 3), num_features)
      else min(ceiling(mtry_factor * sqrt(num_features)), num_features)
      ntree <- max(num_features * ntree_factor, min_ntree)
      nodesize_i <- if (CLASSIFICATION) 1 else 5
      target <- ceiling(num_features * keep_fraction)
      is_initial <- TRUE
      
      # begins RFE for specific module
      initial_screen_mod <- NULL
      survivors_i <- NULL
      
      while (num_features >= target) {
        
        # debug warning message
        if (debug != -1 || debug != 1){
          # if module has low features, only keep non zero important features
          if (num_features <= min_features){
            if (verbose != 0){
              warning(sprintf("\n\n Module %s has fewer than %d features! 
                                All non-zero important features will be kept during screening.", 
                              module_list[i], min_features))
            }
            keep = TRUE
          } else{
            keep = FALSE
          }
        }
        # random forest
        rf <- ranger(
          formula = y ~ .,
          data = data.frame(y = y, module),
          num.trees = ntree,
          mtry = mtry,
          importance = "permutation",
          min.node.size = nodesize_i,
          keep.inbag = TRUE,
          probability = CLASSIFICATION,
          num.threads = 1
        )
        
        # full shapleyforest
        if (shap_model == "full"){
          # sets ups prediction function for fastshap
          predict_function <- if (CLASSIFICATION) {
            if (num_classes == 2) {
              function(object, newdata) {
                preds <- predict(object, data = newdata)$predictions
                return(preds[, 1])  # ranger returns a matrix with one column: prob of class 1
              }
            } else if (num_classes > 2) {
              function(object, newdata) {
                preds <- predict(object, data = newdata)$predictions
                return(preds)  # returns matrix with class probabilities
              }
            } else {
              stop("Invalid or single-class data in y")
            }
          } else {
            function(object, newdata) {
              preds <- predict(object, data = newdata)$predictions
              return(preds)  # regression: numeric vector
            }
          }
          
          # fastshap
          shap <- suppressMessages(fastshap::explain(
            object = rf,
            X = module,
            pred_wrapper = predict_function,
            nsim = 1,
            adjust = FALSE,
            parallel = FALSE,
            .packages = "ranger"
          ))
          
          # gathers absolute shap values
          var_importance <- colMeans(abs(shap))
          var_importance <- sort(var_importance, decreasing = TRUE)
          var_importance <- data.frame(Feature = var_importance)
        }
        
        # after shapleyforest
        if (shap_model == "after"){
          # stores permutation VIM
          var_importance <- rf$variable.importance
          var_importance <- sort(var_importance, decreasing = TRUE)
          var_importance <- data.frame(Variable = names(var_importance),
                                       Importance = as.vector(var_importance))
        }
        
        # keeps top features and updates new module list
        reduction <- ceiling(num_features * drop_fraction)
        
        # if still above target, trim features
        if (num_features - reduction > target) {
          trimmed_varlist <- var_importance[1:(num_features - reduction), , drop = FALSE]
          features <- rownames(trimmed_varlist)
          module <- module[, names(module) %in% features, drop = FALSE]
          num_features <- length(features)
          
          # update mtry/ntree
          mtry <- if (CLASSIFICATION) min(ceiling(mtry_factor * num_features / 3), num_features)
          else min(ceiling(mtry_factor * sqrt(num_features)), num_features)
          ntree <- max(num_features * ntree_factor, min_ntree)
        } else {
          # RFE stops, record survivors
          num_features <- target - 1
          mod_varlist <- var_importance[, 1][1:target]
          features <- rownames(var_importance)[1:target]
          survivors_i <- data.frame(features = features, mod_varlist = mod_varlist)
        }
        
        
        # keeps initial screening
        if (is_initial) {
          feature_name <- rownames(var_importance)
          survived <- feature_name %in% rownames(trimmed_varlist)
          initial_screen_mod <- data.frame(
            "Module Color" = module_name,
            "Feature Name" = feature_name,
            "VIM" = var_importance[, 1],
            "Survivor" = survived
          )
          is_initial <- FALSE
        }
      }
      
      # maintains output for each cluster during parallelizing
      list(
        survivor = survivors_i,
        screen_df = initial_screen_mod
      )
    }, future.seed = seed)
    
    # parallelizing stops
    plan(sequential)
    
    assign("result", results, envir = .GlobalEnv)
    
    # combines all survivor lists from all modules
    survivor_results <- lapply(results, `[[`, "survivor")
    initial_screen <- do.call(rbind, lapply(results, `[[`, "screen_df"))
    names(survivors) <- module_list
  }
  
  return(list(
    survivor_results = survivor_results,
    initial_screen = initial_screen
  ))
}

shapselect_RF <- function(X, y, drop_fraction, number_selected, CLASSIFICATION, mtry_factor,
                          ntree_factor, min_ntree,
                          num_processors, nodesize, cl, nsim) {
  # initialize lists
  selection_list <- list()
  feature_list <- NULL
  
  # initializes parallelization
  if (num_processors > 1) {
    plan(multisession, workers = num_processors)
    parallel = TRUE
  } else { parallel = FALSE}
  
  # sets parameters for RFE
  num_features <- ncol(X)
  mtry <- min(ceiling(mtry_factor*sqrt(num_features)), dim(X)[2])
  ntree <- max(num_features*ntree_factor, min_ntree)
  target <- number_selected
  current_X <- X
  
  cat(num_features, target)
  
  ## begins Recursive Feature Elimination
  i <- 1
  while (num_features >= target){ # runs until it hits target
    # runs random forest
    rf <- ranger(
      formula = y ~ .,
      data = data.frame(y = y, current_X),
      num.trees = ntree,
      mtry = mtry,
      importance = "permutation",
      min.node.size = nodesize,
      keep.inbag = TRUE,
      probability = CLASSIFICATION,
      num.threads = 1
    )
    ## calculates SHAP values
    # sets ups prediction function for fastshap
    predict_function <- if (CLASSIFICATION) {
      if (num_classes == 2) {
        function(object, newdata) {
          preds <- predict(object, data = newdata)$predictions
          return(preds[, 1])  # ranger returns a matrix with one column: prob of class 1
        }
      } else if (num_classes > 2) {
        function(object, newdata) {
          preds <- predict(object, data = newdata)$predictions
          return(preds)  # returns matrix with class probabilities
        }
      } else {
        stop("Invalid or single-class data in y")
      }
    } else {
      function(object, newdata) {
        preds <- predict(object, data = newdata)$predictions
        return(preds)  # regression: numeric vector
      }
    }
    
    # fastshap
    shap <- suppressMessages(fastshap::explain(
      object = rf,
      X = current_X,
      pred_wrapper = predict_function,
      nsim = 1,
      adjust = FALSE,
      parallel = parallel,
      .packages = "ranger"
    ))
    var_importance <- colMeans(abs(shap))
    var_importance <- sort(var_importance, decreasing = TRUE)
    var_importance <- data.frame(Feature = var_importance)
    
    # adds to selection list for iteration i
    selection_list[[i]] <- data.frame(row.names(var_importance),
                                      round(var_importance[, 1], 4),
                                      stringsAsFactors=FALSE)
    
    names(selection_list[[i]]) <- c("feature_name", "variable_importance")
    i <- i + 1
    
    # sets reduction value for feature elimination
    reduction <- ceiling(num_features*drop_fraction)
    
    # if target not reached, continue RFE
    if(num_features - reduction > target) {
      # keep top variable important features
      row_names_var_importance <- row.names(var_importance)
      trimmed_varlist <- var_importance[1:(num_features - reduction), , drop = FALSE]
      features <- row_names_var_importance[1:(num_features - reduction)]
      row.names(trimmed_varlist) <- features
      current_X <- current_X[, which(names(current_X) %in% features)]
      num_features <- length(features)
      # adjust Random Forest variables for next iteration
      mtry <- min(ceiling(mtry_factor*sqrt(num_features)), dim(current_X)[2])
      ntree <- max(num_features*ntree_factor, min_ntree)
    }
    # if target is reached, stop RFE
    else {
      num_features <- target - 1 # stops while loop
      mod_varlist <- var_importance[, 1][1:target] # module membership for survivors
      features <- row.names(var_importance)[1:target] # feature names for survivors
      
      # extracts survivor list for given module
      feature_list <- cbind(features, mod_varlist)
      selection_list[[i]] <- as.data.frame(cbind(features, round(mod_varlist, 4)),
                                           stringsAsFactors=FALSE)
      names(selection_list[[i]]) <- c("feature_name", "variable_importance")
    }
  }
  
  # stops parallezation
  if(num_processors > 1) {
    plan(sequential)
    registerDoSEQ()  
  }
  
  # output
  out <- list(feature_list, selection_list)
  return(out)
}

# library(compiler)
# shapselect_RF <- cmpfun(shapselect_RF)

shapley_forest <- function(final_rf, final_X, module_membership,
                           WGCNA_object=NULL, survivor_list, selection_list, 
                           final_shap, shap_obj) {
  # define output list
  out <- list()
  
  # outputs
  out[[1]] <- final_rf
  out[[2]] <- final_X
  out[[3]] <- module_membership
  out[[4]] <- WGCNA_object
  out[[5]] <- survivor_list
  out[[6]] <- selection_list
  out[[7]] <- final_shap
  out[[8]] <- shap_obj
  
  # define column name for each output
  names(out) <- c("final_rf", "final_X", "module_membership",
                  "WGCNA_object", "survivor_list", "selection_list",
                  "final_SHAP", "shap_obj")
  
  # defines class of object type
  class(out) <- "shapley_forest"
  
  return(out)
}

set.seed(id)

rep_num <- 100
keep_frac <- c(0.01, 0.05, 0.1, 0.15, 0.25)
#keep_frac <- c(0.01, 0.1, 0.25)
drop_frac <- c(0.05, 0.1, 0.25, 0.5)
#drop_frac <- c(0.05, 0.25, 0.5)
mtry_factor <- c(0.5, 1, 2)
#mtry_factor <- c(1, 2)
p <- c(100, 1000)
n <- c(100)

param_list <- list(keep_frac, drop_frac, mtry_factor, p, n)
param_settings <- expand.grid(param_list)
param_settings <- param_settings[, 5:1]
names(param_settings) <- c("n", "p", "mtry_factor", "drop_fraction", "keep_fraction")

param_settings
current_sim_params <- param_settings[ceiling((id)/rep_num), ]

sim_number <- 5
sim_results <- list()
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

runtime <- system.time({
  for (l in 1:sim_number) {
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
    ff <- shapff(X, y, module_membership = module_membership,
                 select_params = select_params,
                 shap_model = "full",
                 screen_params = screen_params,
                 auto_initial = 2,
                 nodesize = 1,
                 debug = 1,
                 verbose = 2,
                 initial = TRUE,
                 num_processors = 8,
                 seed = 2002)
    shap_feature <- data.frame(index = ff$feature_list[1], 
                               ff$feature_list[2], row.names = NULL)
    vim <- shap_feature
    sim_results[[l]] <- vim
  }
})

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
  sim_out <- cbind(mean_vims, selected_props)
  return(sim_out)
}


out <- list(summarize_sim(sim_results, vim_interest), vim_interest, current_sim_params, runtime = runtime)

dir.create("out", showWarnings = FALSE)
out_filename <- paste("parallel/out", id, sep = "")
save(out, file = out_filename)