library(WGCNA)
#library(randomForest)
library(mvtnorm)
library(fuzzyforest)
library(ranger)
library(future)
library(future.apply)
library(doFuture)

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
  
  ## Screening Step
  survivors <- vector('list', length(module_list)) # initiates survivor list
  
  # parameters for screening
  drop_fraction <- screen_control$drop_fraction 
  mtry_factor <- screen_control$mtry_factor
  ntree_factor <- screen_control$ntree_factor
  min_ntree <- screen_control$min_ntree
  keep_fraction <- screen_control$keep_fraction
  
  # adjusts keep_fraction if below minimum threshold
  if(ncol(X)*keep_fraction < select_control$number_selected){
    if (verbose != 0){
      warning(c("\n\n ncol(X)*keep_fraction < number_selected", "\n",
                "number_selected will be set to floor(ncol(X)*keep_fraction)"))
    }
    select_control$number_selected <- max(floor(ncol(X)*keep_fraction), 1)
  }
  
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
    
  }
  else {
    # starts parallelization
    registerDoFuture() 
    plan(multisession, workers = num_processors)
    
    # parallelizes RFE at each WGCNA module
    survivor_results <- future_lapply(seq_along(module_list), function(i) {
      
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
          probability = CLASSIFICATION
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
        trimmed_varlist <- var_importance[1:(num_features - reduction), , drop = FALSE]
        features <- row.names(trimmed_varlist)
        module <- module[, which(names(module) %in% features)]
        num_features <- length(features)
        
        # if RFE reaches target amount of features remaining, stores survivors
        if (num_features <= target) {
          mod_varlist <- var_importance[, 1][1:target]
          features <- row.names(var_importance)[1:target]
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
  }
  
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
  
  # combines all survivor lists from all modules
  survivors <- lapply(survivor_results, `[[`, "survivor")
  initial_screen <- do.call(rbind, lapply(survivor_results, `[[`, "screen_df"))
  names(survivors) <- module_list
  
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
  
  ## Selection Step
  
  # sets up selection step parameters
  select_args <- list(X_surv, y, num_processors, nodesize)
  select_args <- c(select_args, select_control)
  names(select_args)[1:4] <- c("X", "y", "num_processors", "nodesize")
  
  # RFE via fastshap
  if (shap_model == "full"){
    select_results <- shapselect_RF(select_args$X, select_args$y, select_args$drop_fraction, 
                                    select_args$number_selected, CLASSIFICATION, select_args$mtry_factors,
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
    verbose = FALSE
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

shapselect_RF <- function(X, y, drop_fraction, number_selected, CLASSIFICATION, mtry_factor,
                          ntree_factor, min_ntree,
                          num_processors, nodesize, cl, nsim) {
  # initialize lists
  selection_list <- list()
  feature_list <- NULL
  
  # initializes parallelization
  if (num_processors > 1) {
    plan(multisession, workers = num_processors)
  }
  
  # sets parameters for RFE
  num_features <- ncol(X)
  mtry <- min(ceiling(mtry_factor*sqrt(num_features)), dim(X)[2])
  ntree <- max(num_features*ntree_factor, min_ntree)
  target <- number_selected
  current_X <- X
  
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
      probability = CLASSIFICATION
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
      parallel = FALSE,
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


#########################
#---------data----------#
#########################

visit_a_x <- qread("/u/home/t/tzshen/hiv-visit-ag/ewas/data/visit_a_x_meanvar.qs")
visit_a_y <- qread("/u/home/t/tzshen/hiv-visit-ag/ewas/data/Y_visit_a_ewas.qs")

wgcna_module_info <- read_csv("/u/home/t/tzshen/hiv-visit-ag/ewas/data/module_membership_signed_kME_mean.csv", 
                              show_col_types = FALSE)

dynamicColors <- wgcna_module_info$PrimaryModule

#########################
#---------param---------#
#########################

screen_params <- screen_control(
  drop_fraction = 0.5,         # Drop half the features at each screening step
  keep_fraction = 0.25,        # Only keep 25% of features to carry forward
  mtry_factor = 0.25           # Fewer variables considered at each tree split
)

select_params <- select_control(
  number_selected = 10,        # Small number of final selected features
  drop_fraction = 0.5,         # Drop aggressively during selection
  mtry_factor = 0.25           # Keep tree building light
)


#########################
#-----shapleyforest-----#
#########################

n_cores <- parallel::detectCores() - 1

set.seed(2002)

X <- visit_a_x

y <- visit_a_y %>% select(casecontrol)

y <- as.factor(y$casecontrol)

sff <- shapff(X, y, module_membership = dynamicColors,
              select_params = select_params,
              shap_model = "full",
              screen_params = screen_params,
              auto_initial = 2,
              nodesize = 1,
              debug = 1,
              verbose = 1,
              initial = TRUE,
              num_processors = n_cores,
              min_features = 10,
              seed = 2002)

saveRDS(sff, "/u/home/t/tzshen/hiv-visit-ag/ewas/data/shapforest_initial.rds")


