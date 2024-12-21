#' Runs shapley forest algorithm
#'
#' Runs shapley forest algorithm for feature importance through
#' the use of SHAPley values.
#' 
#' @export
#' @param X                 A data.frame. where each column represents a
#'                          feature vector.
#' @param y                 Response vector. If performing classification, `y` should 
#'                          be a factor. If performing regression, `y`
#'                          should be numeric vector.
#' @param Z                 A data.frame of additional features that will bypass
#'                          the screening step.
#' @param shap_model        Binary indicator for \code{shapff} model. If `full`, \code{shapff}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `after`, \code{shapff} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `1` is default.
#' @param module_membership A vector that specifies the module membership for each
#'                          each feature. See \code{shapwff} for possible method.
#' @param min_features      Defines minimum feature allowed for each module. If `debug` is 
#'                          not `-1`, modules below `min_features` will only keep non-zero
#'                          important features during each Recursive Feature Elimination
#'                          iteration
#' @param verbose           Defines the warning message protocol. If `0`, no warning or UI
#'                          will be displayed. If `1`, warnings and UI progress bar will
#'                          be displayed.
#' @param debug             Sets the debugging procedures. If `-1`, all debugging functions
#'                          will be bypassed. If `0`, debugging at the WGCNA will be bypassed.
#'                          Note for \code{shapff}, `0` has no effect. If `1`, debugging during
#'                          Recursive Feature Elimination at both screening and selection step 
#'                          will be bypassed. If `2`, all debugging functions will be ran. Below
#'                          are the debugging features. Debugging at WGCNA detects if each module
#'                          is below the \code{min_features}. Debugging at RFE will keep only 
#'                          non zero important feature at each elimination step for modules below
#'                          \code{min_features}.
#' @param initial           Binary indicator to print out initial screening step results (ie the
#'                          results from the first Recursive Feature Elimination at the screening
#'                          step for each module). If `True`, \code{shapff} will pause after RFE
#'                          allowing users to select output method for initial screening. If `False`,
#'                          it will bypass all initial screening procedure.
#' @param auto_initial      Bypass readline prompt for \code{initial}. If `1`, `initial_screening.csv`
#'                          will be saved in directory and stops. If `2`, `initial_screening.csv`
#'                          will be saved in directory and proceeds. If `3`, nothing saved and stops.
#'                          If `4`, nothing saved and proceeds. Default is `NULL`. Note if \code{initial},
#'                          is set to `TRUE`, `auto_initial` will automically be set to `NULL`.
#' @param screen_params     Defines the parameter settings for the screening step
#'                          of \link[fuzzyforest]{fuzzyforest}.
#'                          See \code{\link[fuzzyforest]{screen_control}} for
#'                          details. \code{screen_params} is an object of type
#'                          \code{screen_control}.
#' @param select_params     Defines the parameter setting for the selection step
#'                          of \link[fuzzyforest]{fuzzyforest}.
#'                          See \code{\link[fuzzyforest]{select_control}} for details.
#'                          \code{select_params} is an object of type
#'                          \code{select_control}.
#' @param final_ntree       The number of trees grown in the final random forest in
#'                          the selection step. This random forest contains all
#'                          the surviving features.
#' @param num_processors    Number of processors used to run random forests.
#' @param parallel          Type of parellization to be used. `1` if 
#'                          \code{\link[doParallel]{doParallel}}. `2` if 
#'                          \code{\link[doSNOW]{doSNOW}}. `1` is the default.
#' @param nodesize          Minimum terminal nodesize. 1 if classification.
#'                          5 if regression.  If the sample size is very large,
#'                          the trees will be grown extremely deep.
#'                          This may lead to issues with memory usage and may
#'                          lead to significant increases in the time it takes
#'                          the algorithm to run. In this case,
#'                          it may be useful to increase \code{nodesize}.
#' @param test_features     A data.frame containing features from a test set.
#'                          The data.frame should contain the features in both
#'                          X and Z. Used during final Random Forest call (after
#'                          screening and selection step).
#' @param test_y            The responses for the test set. Used during final Random 
#'                          Forest call (after screening and selection step).
#' @param nsim              Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the screening step. Default is `1`. Increasing
#'                          \code{nsim} leads to more accurate results, but at the cost
#'                          of computational cost.
#' @param final_nsim        Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the selection step. Default is `1`. \code{final_nsim}
#'                          should be as large as feasibly possible.
#' @return Returns an object of type `shapley_forest`, which is a list containing the essential 
#' output of shapley forests, including a data.frame of selected features and the random forest 
#' model fitted using those features.
#' 
#' @import fuzzyforest
#' @import randomForest
#' @import fastshap
#' 
#' @importFrom randomForest combine


shapff <- function(X, y, Z=NULL, shap_model = "full", module_membership,
                   min_features = 20, verbose = 1, debug = 2, 
                   initial = TRUE, auto_initial = NULL, 
                   screen_params = fuzzyforest:::screen_control(min_ntree=5000),
                   select_params = fuzzyforest:::select_control(min_ntree=5000),
                   final_ntree = 5000,
                   num_processors = 1, parallel = 1, nodesize, 
                   test_features=NULL, test_y=NULL, nsim = 1, 
                   final_nsim = 100) {
  
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
    options("warn" = -1)
  }
  
  # sets parameters for each step
  screen_control <- screen_params
  select_control <- select_params
  
  # obtains module membership (ie. from WGCNA)
  module_list <- unique(module_membership)
  
  # sets up paralleling via DoParallel
  if(parallel == 1 && num_processors > 1) {
    cl = parallel::makeCluster(num_processors)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }
  # sets up paralleling via DoSNOW
  if(parallel == 2 && num_processors > 1) {
    cl <- snow::makeCluster(num_processors, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    on.exit(snow::stopCluster(cl), add = TRUE)
  }
  
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
    warning(c("ncol(X)*keep_fraction < number_selected", "\n",
              "number_selected will be set to floor(ncol(X)*keep_fraction)"))
    select_control$number_selected <- max(floor(ncol(X)*keep_fraction), 1)
  }
  
  # begins RFE
  total_iterations <- length(module_list)
  
  # sets UI Progress Bar
  if (verbose !=0){
    progress_bar <- txtProgressBar(min = 0, max = total_iterations, style = 3)
  }
  
  # initiates initial_screen module list
  initial_screen <- data.frame("Module Color" = character(0), "Feature Name" = character(0),
                    "VIM" = numeric(0), "Survivor" = logical(0))
  
  # Recursive Feature Elimination
  for (i in 1:length(module_list)) {
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
          warning(sprintf("Module %s has fewer than %d features! All non-zero important features will be kept during screening.", 
                          module_list[i], min_features))
          keep = TRUE
        } else{
          keep = FALSE
        }
      }
      
      # runs Random Forest
      if(num_processors > 1) {
        rf = `%dopar%`(foreach(ntree = rep(ntree/num_processors, num_processors)
                               , .combine = randomForest::combine, .packages = 'randomForest'),
                       #second argument to '%dopar%'
                       randomForest(module , y, ntree = ntree, mtry = mtry,
                                    importance = TRUE, scale = FALSE, nodesize=nodesize))
      }
      if(num_processors == 1) {
        rf <- randomForest(module, y, ntree = ntree, mtry = mtry,
                           importance = TRUE, scale = FALSE,
                           nodesize = nodesize)
      }
      
      # run Feature Selection
      
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
          shap <- suppressMessages(fastshap::explain(rf, X = module, nsim = nsim, pred_wrapper = prediction))
        }
        
        # regression case
        if (CLASSIFICATION == FALSE){
          shap <- suppressMessages(fastshap::explain(rf, X = module, nsim = nsim, pred_wrapper = predict))
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
        trimmed_varlist <- var_importance[1:(num_features - reduction), ,drop=FALSE]
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
        temp <- data.frame("Module Color" = mod_name,"Feature Name" = feature_name, "VIM" = var_importance[,1], 
                           "Survivor" = survived)
        initial_screen <- rbind(initial_screen, temp)
        is_initial = FALSE
      }
    }

  }
  
  
  ## Initial Screening Output Procedure
  if (initial == TRUE){
    cat("\nDisplaying Initial Screen (First Step) Variable Importance ... ")
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
      cat("Dataframe saved as 'initial_screen.csv'.\n")
      if (tolower(save_prompt) == "1"){ # stops running
        options(warn = 1)
        stop("Execution stopped as per user choice.\n")
      }
    } else { # skips all
      if (tolower(save_prompt) == "3"){
        options(warn = 1)
        stop("Execution stopped as per user choice.\n")
      }
    }
  }
  
  if (auto_initial != NULL){
    if (auto_initial %in% c("1", "2")) { # save output
      write.csv(initial_screen, "initial_screen.csv", row.names = FALSE)
      assign("initial_screen", initial_screen, envir = .GlobalEnv)
      cat("Dataframe saved as 'initial_screen.csv'.\n")
      if (auto_initial == "1"){ # stops running
        options(warn = 1)
        stop("Execution stopped as per user choice.\n")
      }
    } else { # skips all
      if (auto_initial == "3"){
        options(warn = 1)
        stop("Execution stopped as per user choice.\n")
      }
    }
  }
  
  # verbose UI
  if (verbose != 0){cat("\nSelection Step ...")}
  
  # combines all survivor lists from all modules
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
  
  # gathers final surivoring list
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
  
  # final Random Forest
  final_rf <- randomForest(x=final_X, y=y, mtry=final_mtry, ntree=final_ntree,
                           importance=TRUE, nodesize=nodesize,
                           xtest=test_features, ytest=test_y)
  
  # extracts module membership of final survivors
  final_module_membership <- as.data.frame(cbind(names(X), module_membership),
                                           stringsAsFactors=FALSE)
  names(final_module_membership) <- c("feature_name", "module")
  
  # final variable importance measure via fastshap
  shap_type = "shapley"
  if (shap_type =="shapley"){
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
        shap_final_obj <- fastshap::explain(final_rf, X = final_X, nsim = final_nsim, 
                                            pred_wrapper = prediction, shap_only = FALSE)
    }
    
    # regression case
    if (CLASSIFICATION == FALSE){
      shap_final_obj <- fastshap::explain(final_rf, X = final_X, nsim = final_nsim, 
                                          pred_wrapper = predict, shap_only = FALSE)
    }
    shap_final <- shap_final_obj$shapley_values
    var_importance_final <- colMeans(abs(shap_final))
    var_importance_final <- sort(var_importance_final, decreasing = TRUE)
    var_importance_final <- data.frame(vim = var_importance_final)
  }
  
  #if (shap_type == "tree"){
  #  unified_model <- randomForest.unify(final_rf, final_X)
  #  shap_final_obj <- treeshap::treeshap(unified_model, x = final_X, 
  #                                       interactions = TRUE, verbose = FALSE)
  #  shap_final <- shap_final_obj$shaps
  #  var_importance_final <- colMeans(abs(shap_final))
  #  var_importance_final <- sort(var_importance_final, decreasing = TRUE)
  #  var_importance_final <- data.frame(vim = var_importance_final)
  #}
  
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
  out <- shap_fuzzy_forest(final_rf, final_module_membership,
                           survivor_list=survivor_list,
                           selection_list=selection_list, 
                           final_shap = shap_final_list,
                           shap_obj = shap_final_obj,
                           final_X = final_X)
  
  # verbose message
  if (verbose != 0){cat("Done \n")}
  
  return(out)
}

#' Runs shapley forest algorithm along with WGCNA
#'
#' Runs a weighted gene correlated network analysis. Then runs shapley forest 
#' algorithm for feature importance through the use of SHAPley values.
#' 
#' @export
#' @param X                 A data.frame. where each column represents a
#'                          feature vector.
#' @param y                 Response vector. If performing classification, `y` should 
#'                          be a factor. If performing regression, `y`
#'                          should be numeric vector.
#' @param Z                 A data.frame of additional features that will bypass
#'                          the screening step.
#' @param WGCNA_params      WGCNA parameters.
#'                          See \code{\link[WGCNA]{blockwiseModules}} and
#'                          \code{\link[fuzzyforest]{WGCNA_control}} for details.
#'                          \code{WGCNA_params} is an object of type
#'                          \code{WGCNA_control}.
#' @param shap_model        Binary indicator for \code{shapff} model. If `full`, \code{shapff}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `after`, \code{shapff} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `1` is default.
#' @param module_membership A vector that specifies the module membership for each
#'                          each feature. See \code{shapwff} for possible method.
#' @param min_features      Defines minimum feature allowed for each module. If `debug` is 
#'                          not `-1`, modules below `min_features` will only keep non-zero
#'                          important features during each Recursive Feature Elimination
#'                          iteration
#' @param verbose           Defines the warning message protocol. If `0`, no warning or UI
#'                          will be displayed. If `1`, warnings and UI progress bar will
#'                          be displayed.
#' @param debug             Sets the debugging procedures. If `-1`, all debugging functions
#'                          will be bypassed. If `0`, debugging at the WGCNA will be bypassed.
#'                          Note for \code{shapff}, `0` has no effect. If `1`, debugging during
#'                          Recursive Feature Elimination at both screening and selection step 
#'                          will be bypassed. If `2`, all debugging functions will be ran. Below
#'                          are the debugging features. Debugging at WGCNA detects if each module
#'                          is below the \code{min_features}. Debugging at RFE will keep only 
#'                          non zero important feature at each elimination step for modules below
#'                          \code{min_features}.
#' @param initial           Binary indicator to print out initial screening step results (ie the
#'                          results from the first Recursive Feature Elimination at the screening
#'                          step for each module). If `True`, \code{shapff} will pause after RFE
#'                          allowing users to select output method for initial screening. If `False`,
#'                          it will bypass all initial screening procedure.
#' @param auto_initial      Bypass readline prompt for \code{initial}. If `1`, `initial_screening.csv`
#'                          will be saved in directory and stops. If `2`, `initial_screening.csv`
#'                          will be saved in directory and proceeds. If `3`, nothing saved and stops.
#'                          If `4`, nothing saved and proceeds. Default is `NULL`. Note if \code{initial},
#'                          is set to `TRUE`, `auto_initial` will automically be set to `NULL`.
#' @param screen_params     Defines the parameter settings for the screening step
#'                          of \link[fuzzyforest]{fuzzyforest}.
#'                          See \code{\link[fuzzyforest]{screen_control}} for
#'                          details. \code{screen_params} is an object of type
#'                          \code{screen_control}.
#' @param select_params     Defines the parameter setting for the selection step
#'                          of \link[fuzzyforest]{fuzzyforest}.
#'                          See \code{\link[fuzzyforest]{select_control}} for details.
#'                          \code{select_params} is an object of type
#'                          \code{select_control}.
#' @param final_ntree       The number of trees grown in the final random forest in
#'                          the selection step. This random forest contains all
#'                          the surviving features.
#' @param num_processors    Number of processors used to run random forests.
#' @param parallel          Type of parellization to be used. `1` if 
#'                          \code{\link[doParallel]{doParallel}}. `2` if 
#'                          \code{\link[doSNOW]{doSNOW}}. `1` is the default.
#' @param nodesize          Minimum terminal nodesize. 1 if classification.
#'                          5 if regression.  If the sample size is very large,
#'                          the trees will be grown extremely deep.
#'                          This may lead to issues with memory usage and may
#'                          lead to significant increases in the time it takes
#'                          the algorithm to run. In this case,
#'                          it may be useful to increase \code{nodesize}.
#' @param test_features     A data.frame containing features from a test set.
#'                          The data.frame should contain the features in both
#'                          X and Z. Used during final Random Forest call (after
#'                          screening and selection step).
#' @param test_y            The responses for the test set. Used during final Random 
#'                          Forest call (after screening and selection step).
#' @param nsim              Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the screening step. Default is `1`. Increasing
#'                          \code{nsim} leads to more accurate results, but at the cost
#'                          of computational cost.
#' @param final_nsim        Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the selection step. Default is `1`. \code{final_nsim}
#'                          should be as large as feasibly possible.
#' @return Returns an object of type fuzzy_forest, which is a list containing the essential 
#' output of fuzzy forests, including a data.frame of selected features and the random forest 
#' model fitted using those features.
#' 
#' @import fuzzyforest
#' @import randomForest
#' @import fastshap
#' 
# @importFrom randomForest combine
#' 
#' @references
#' Leo Breiman (2001). Random Forests. Machine Learning, 45(1), 5-32.
#'
#' Lundberg, S. M., & Lee, S. I. (2017). A unified approach to interpreting model predictions. Advances in neural information processing systems, 30.
#'
#' Daniel Conn, Tuck Ngun, Christina M. Ramirez (2015). Fuzzy Forests: a New
#' WGCNA Based Random Forest Algorithm for Correlated, High-Dimensional Data,
#' Journal of Statistical Software, Manuscript in progress.
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#' @examples
#' TODO
#' 

shapwff <- function(X, y, Z=NULL, shap_model = "full",
                    WGCNA_params=WGCNA_control(p=6),
                    min_features=20, verbose = 1, debug = 2, 
                    initial = TRUE, auto_initial = NULL,
                    screen_params=fuzzyforest:::screen_control(min_ntree=5000),
                    select_params=fuzzyforest:::select_control(min_ntree=5000),
                    final_ntree=500, num_processors, parallel=1, nodesize,
                    test_features=NULL, test_y=NULL, nsim=1, final_nsim=100) {
  
  if ( !("package:WGCNA" %in% search()) ) {
    stop("WGCNA must be loaded and attached. Type library(WGCNA) to do so.",
         call. = FALSE)
  }
  if (!(is.vector(y) || is.factor(y))) {
    stop("y must be vector or factor")
  }
  # WGCNA yields errors if X is an integer, here we convert integers to numeric.
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
  if (shap_model != "full" && shap_model != "after") {
    stop("shap_model must be `full` or `after`. Type help(shapff) or help(shapwff) for details.")
  }
  
  if (!verbose %in% c(0, 1)) {
    stop("verbose must be 0 or 1")
  }
  
  if (!debug %in% c(-1, 0, 1, 2)) {
    stop("debug must be -1, 0, 1, or 2")
  }
  
  if (verbose == 0 && !debug %in% c(-1, 0)){
    stop("if debug is 1 or 2, verbose must be 1")
  }
  
  if (!is.logical(initial) || length(initial) != 1) {
    stop("initial must be boolean.")
  }
  
  if (!is.null(auto_initial) && !auto_initial %in% c(1, 2, 3, 4) ){
    stop("auto_initial must be NULL, 1, 2, 3, or 4")
  }
  
  if (initial == FALSE){
    auto_initial = NULL
  }
  
  if (verbose == 0){
    options("warn" = -1)
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
  min_features <- min_features
  
  if (!debug %in% c(-1, 0)){
    feature_counts <- table(module_membership)
    low_frequency_modules <- feature_counts[feature_counts <= min_features]
    
    if (length(low_frequency_modules) > 0) {
      warning(sprintf("WGCNA - Some modules contain fewer than % s features.", min_features))
      response <- readline(prompt = "Do you wish to continue? (yes/no): ")
      if (tolower(response) != "yes") {
        stop(cat(sprintf("Process terminated by the user. Low Frequency Modules:\n%s", 
                         paste(capture.output(print(low_frequency_modules)), 
                               collapse = "\n")), "\n"))
      }
    }
  }
  
  if (verbose != "0"){ cat("Screening Step ... \n")}
  out <- shapff(X, y, Z, shap_model, module_membership,
                min_features, verbose, debug, initial, screen_control, 
                select_control, final_ntree,
                num_processors, nodesize=nodesize,
                test_features=test_features, test_y=test_y)
  out$WGCNA_object <- bwise
  
  options(warn = 1)
  return(out)
}
