#' Screening step of shapley forest through SHAP values.
#'
#' Runs the screening step of shapley forest algorithm.
#' Returns data.frame with the final surviving SHAP values and features.
#' 
#' @param X                 A data.frame. With columns denoting a feature vector.
#'                          Could include additional covariates not a part of
#'                          the original modules.
#' @param y                 A Response vector.
#' @param module_list       A list containing the modules. 
#' @param module_membership A vector that specifies the module membership for each
#'                          each feature. See \code{wsf} for possible method.
#' @param screen_control    Defines the parameter settings for the screening step.
#' @param select_control    Defines the parameter setting for the selection step.
#' @param shap_model        Binary indicator for \code{sf} model. If `full`, \code{sf}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `after`, \code{sf} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `full` is default.
#' @param CLASSIFICATION    Binary indicator if the model is a classification model.
#' @param min_features      Defines minimum feature allowed for each module. 
#' @param nsim              Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the screening step. Default is `1`. Increasing
#'                          \code{nsim} leads to more accurate results, but at the cost
#'                          of computational cost.
#' @param nodesize          Minimum `nodesize`
#' @param num_processors    Number of processors used to fit random forests.
#' @param verbose           Defines the warning message protocol. See \code{sf} for information.
#' @param debug             Sets the debugging procedures. See \code{sf} for more information
#' @param parallel          Binary indicator if parallel is called. 
#' @param seed              Seed to be set for reproducibility. Note: for a given seed,
#'                          parallel and non-parallel versions will differ.
#' 
#' @import fuzzyforest
#' @import ranger
#' @import fastshap
#' 
#' @importFrom randomForest combine
#' 
#' @return A data.frame with final surviving features.
#' 
screen_RFE <- function(X, y, module_list, module_membership, screen_control, 
                       select_control, shap_model = "full", 
                       CLASSIFICATION, 
                       min_features, nsim, nodesize, 
                       num_processors, verbose, 
                       debug, parallel,
                       seed) {
  
  # sets seeds for each module
  module_seeds <- seed + seq_along(module_list)
  
  ## Screening Step
  survivors <- vector('list', length(module_list)) # initiates survivor list
  screen_df <- vector('list', length(module_list)) # initiates survivor list
  
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
    for (i in 1:length(module_list)) {
      
      # set seed
      seed <- module_seeds[i]
      set.seed(seed)
      
      # progress bar
      if (verbose != 0){
        setTxtProgressBar(progress_bar, i)
      }
      
      # extracts a given module
      module_name <- module_list[i]
      module <- X[, which(module_membership == module_name), drop=FALSE]
      num_features <- ncol(module)
      
      # sets mtry and nodesize for Random Forest
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
      
      # sets ntree for Random Forest
      ntree <- max(num_features*ntree_factor, min_ntree)
      
      # sets target variable to stop RFE
      target = ceiling(num_features * keep_fraction)
      
      # flag for initial screening step
      is_initial = TRUE
      
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
        rf <- ranger(
          formula = y ~ .,
          data = data.frame(y = y, module),
          num.trees = ntree,
          mtry = mtry,
          importance = "permutation",
          min.node.size = nodesize,
          keep.inbag = TRUE,
          probability = CLASSIFICATION,
          num.threads = 1,
          seed = seed
        )
        
        # full shapleyforest
        if (shap_model == "full"){
          # sets ups prediction function for fastshap
          predict_function <- if (CLASSIFICATION) {
            num_classes <- if (CLASSIFICATION) length(levels(y)) else NULL
            if (num_classes == 2) {
              function(object, newdata) {
                preds <- predict(object, data = newdata, seed = seed)$predictions
                return(preds[, 1])  # ranger returns a matrix with one column: prob of class 1
              }
            } else if (num_classes > 2) {
              function(object, newdata) {
                preds <- predict(object, data = newdata, seed = seed)$predictions
                return(preds)  # returns matrix with class probabilities
              }
            } else {
              stop("Invalid or single-class data in y")
            }
          } else {
            function(object, newdata) {
              preds <- predict(object, data = newdata, seed = seed)$predictions
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
          var_importance <- data.frame(Feature = var_importance)
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
          screen_df[[i]] <- data.frame("Module Color" = mod_name,
                                       "Feature Name" = feature_name, 
                                       "VIM" = var_importance[,1], 
                                       "Survivor" = survived)
          is_initial = FALSE
        }
      }
    }
    # Convert non-parallel outputs into same structure as parallel
    survivor_results <- vector("list", length(survivors))
    for (i in seq_along(survivors)) {
      module_name <- module_list[i]
      survivor_results[[i]] <- list(
        survivor = survivors[[i]],
        screen_df = screen_df[[i]],
        seed = seed
      )
    }
  } else {
    # starts parallelization
    registerDoFuture() 
    plan(multisession, workers = num_processors)
    
    # parallelizes RFE at each WGCNA module
    survivor_results <- future_lapply(seq_along(module_list), function(i) {
      seed <- module_seeds[i]
      set.seed(seed)
      
      # gets specific WGCNA module
      module_name <- module_list[i]
      module <- X[, which(module_membership == module_name), drop = FALSE]
      num_features <- ncol(module)
      
      # classification or regression set up for parameters
      CLASSIFICATION <- is.factor(y)
      mtry <- if (CLASSIFICATION) min(ceiling(mtry_factor * num_features / 3), num_features)
      else min(ceiling(mtry_factor * sqrt(num_features)), num_features)
      nodesize_i <- if (CLASSIFICATION) 1 else 5
      
      # sets ntree for Random Forest
      ntree <- max(num_features * ntree_factor, min_ntree)
      
      # sets target variable to stop RFE
      target <- ceiling(num_features * keep_fraction)
      
      # flag for initial screening step
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
          num.threads = 1,
          seed = seed
        )
        
        # full shapleyforest
        if (shap_model == "full"){
          # sets ups prediction function for fastshap
          predict_function <- if (CLASSIFICATION) {
            num_classes <- if (CLASSIFICATION) length(levels(y)) else NULL
            if (num_classes == 2) {
              function(object, newdata) {
                preds <- predict(object, data = newdata, seed = seed)$predictions
                return(preds[, 1])  # ranger returns a matrix with one column: prob of class 1
              }
            } else if (num_classes > 2) {
              function(object, newdata) {
                preds <- predict(object, data = newdata, seed = seed)$predictions
                return(preds)  # returns matrix with class probabilities
              }
            } else {
              stop("Invalid or single-class data in y")
            }
          } else {
            function(object, newdata) {
              preds <- predict(object, data = newdata, seed = seed)$predictions
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
          var_importance <- data.frame(Feature = var_importance)
        }
        
        # keeps top features and updates new module list
        reduction <- ceiling(num_features * drop_fraction)
        
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
        if (num_features - reduction > target) {
          # keep top variable important features
          trimmed_varlist <- var_importance[1:(num_features - reduction), , drop = FALSE]
          features <- rownames(trimmed_varlist)
          module <- module[, names(module) %in% features, drop = FALSE]
          num_features <- length(features)
          
          # adjust Random Forest variables for next iteration
          mtry <- if (CLASSIFICATION) min(ceiling(mtry_factor * num_features / 3), num_features)
          else min(ceiling(mtry_factor * sqrt(num_features)), num_features)
          
          ntree <- max(num_features * ntree_factor, min_ntree)
          
          # if target is reached, stop RFE
        } else {
          num_features <- target - 1 # stops while loop
          mod_varlist <- var_importance[, 1][1:target] # module membership for survivors
          features <- rownames(var_importance)[1:target] # feature name for survivors
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
        screen_df = initial_screen_mod,
        seed = seed
      )
    }, future.seed = NULL)
    
    # parallelizing stops
    plan(sequential)
  }
  
  initial_screen_list <- lapply(survivor_results, function(x) x$screen_df)
  initial_screen <- do.call(rbind, initial_screen_list)
  
  return(list(
    survivor_results = survivor_results,
    initial_screen = initial_screen
  ))
}

#' Selection step of shapley forest through SHAP values.
#'
#' Runs the selection step of shapley forest algorithm.
#' Returns data.frame with the final surviving SHAP values and features.
#' @param X                 A data.frame. With columns denoting a feature vector.
#'                          Could include additional covariates not a part of
#'                          the original modules.
#' @param y                 A Response vector.
#' @param drop_fraction     A value between 0 and 1 representing the percentage of features
#'                          to be dropped at each iteration.
#' @param shap_model        Binary indicator for \code{sf} model. If `full`, \code{sf}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `after`, \code{sf} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `full` is default.
#' @param number_selected   Number of features selected at the end of shapley forest.
#' @param CLASSIFICATION    Binary indicator if the model is a classification model.
#' @param mtry_factor       Number to adjust \code{mtry} for random forest.
#'                          If regression, \code{mtry} is set to
#'                          \code{ceiling}(\eqn{\sqrt(p)}*\code{mtry_factor}).
#'                          If classication, \code{mtry} is set to
#'                          \code{ceiling}((p/3)*\code{mtry_factor}).  If either
#'                          of these numbers is greater than p, \code{mtry} is
#'                          set to p.
#' @param ntree_factor      A number greater than 1 to adjust \code{ntree}.  \code{ntree} for each
#'                          random is \code{ntree_factor} times the number
#'                          of features.  For each random forest, \code{ntree}
#'                          is set to \code{max}(\code{min_ntree},
#'                          \code{ntree_factor}*\code{p}).
#' @param min_ntree         Minimum number of trees grown in each random forest.
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum `nodesize`
#' @param nsim              Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the screening step. Default is `1`. Increasing
#'                          \code{nsim} leads to more accurate results, but at the cost
#'                          of computational cost.
#' @param seed              Seed to be set for reproducibility. Note: for a given seed,
#'                          parallel and non-parallel versions will differ.
#' @import fuzzyforest
#' @import ranger
#' @import fastshap
#' 
#' @importFrom randomForest combine
#' 
#' @return A data.frame with final surviving features.

select_RFE <- function(X, y, drop_fraction, shap_model = "full",
                       number_selected, CLASSIFICATION, mtry_factor,
                       ntree_factor, min_ntree,
                       num_processors, nodesize, nsim, seed) {
  # initialize lists
  selection_list <- list()
  feature_list <- NULL
  
  # initializes parallelization
  if (num_processors > 1) {
    plan(multisession, workers = num_processors)
    parallel = TRUE
  } else { 
    parallel = FALSE
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
      probability = CLASSIFICATION,
      num.threads = 1,
      seed = seed
    )
    
    if (shap_model == "full"){
      ## calculates SHAP values
      # sets ups prediction function for fastshap
      predict_function <- if (CLASSIFICATION) {
        num_classes <- nlevels(y)
        if (num_classes == 2) {
          function(object, newdata) {
            preds <- predict(object, data = newdata, seed = seed)$predictions
            return(preds[, 1])  # ranger returns a matrix with one column: prob of class 1
          }
        } else if (num_classes > 2) {
          function(object, newdata) {
            preds <- predict(object, data = newdata, seed = seed)$predictions
            return(preds)  # returns matrix with class probabilities
          }
        } else {
          stop("Invalid or single-class data in y")
        }
      } else {
        function(object, newdata) {
          preds <- predict(object, data = newdata, seed = seed)$predictions
          return(preds) 
        }
      }
      
      set.seed(seed)
      
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
    }
    
    if (shap_model == "after"){
      var_importance <- rf$variable.importance
      var_importance <- sort(var_importance, decreasing = TRUE)
      var_importance <- data.frame(Feature = var_importance)
    }
    
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