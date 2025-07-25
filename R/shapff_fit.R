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
#' @param shap_model        Binary indicator for \code{sf} model. If `full`, \code{sf}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `after`, \code{sf} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `full` is default.
#' @param module_membership A vector that specifies the module membership for each
#'                          each feature. See \code{wsf} for possible method.
#' @param min_features      Defines minimum feature allowed for each module. If `debug` is 
#'                          not `-1`, modules below `min_features` will only keep non-zero
#'                          important features during each Recursive Feature Elimination
#'                          iteration
#' @param verbose           Defines the warning message protocol. If `0`, no warning or UI
#'                          will be displayed. If `1`, warnings and UI progress bar will
#'                          be displayed.
#' @param debug             Sets the debugging procedures. If `-1`, all debugging functions
#'                          will be bypassed. If `0`, debugging at the WGCNA will be bypassed.
#'                          Note for \code{sf}, `0` has no effect. If `1`, debugging during
#'                          Recursive Feature Elimination at both screening and selection step 
#'                          will be bypassed. If `2`, all debugging functions will be ran. Below
#'                          are the debugging features. Debugging at WGCNA detects if each module
#'                          is below the \code{min_features}. Debugging at RFE will keep only 
#'                          non zero important feature at each elimination step for modules below
#'                          \code{min_features}.
#' @param initial           Binary indicator to print out initial screening step results (ie the
#'                          results from the first Recursive Feature Elimination at the screening
#'                          step for each module). If `True`, \code{sf} will pause after RFE
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
#' @param num_processors    Number of processors used to run shapley forests. If 
#'                          `num_processors` is greater than 1, shapley forest will
#'                          parallelize at each module in the screening step. Furthermore,
#'                          it will call `fastshap` and `ranger` parallel 
#'                          functions at the selection step. See their respective 
#'                          documentation.
#'                          NOTE: Setting `num_processors > 1` may cause runtimes
#'                          to slow with small datasets (ie n,p < 100).
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
#' @param seed              RNG seed for reproducibility. Default is based on current
#'                          system time (`as.integer(Sys.time()))`).
#' @return Returns an object of type `shapley_forest`, which is a list containing the essential 
#' output of shapley forests, including a data.frame of selected features and the random forest 
#' model fitted using those features. See \code{shapley_forest} for more details.
#' 
#' @import fuzzyforest
#' @import ranger
#' @import fastshap
#' 
#' 
#' @references
#' TO DO
sf <- function(X, y, Z=NULL, shap_model = "full", module_membership,
               min_features = 20, verbose = 1, debug = 2, 
               initial = TRUE, auto_initial = NULL, 
               screen_params = fuzzyforest:::screen_control(min_ntree=5000),
               select_params = fuzzyforest:::select_control(min_ntree=5000),
               final_ntree = 5000,
               num_processors = 1, nodesize, 
               test_features=NULL, test_y=NULL, nsim = 1, 
               final_nsim = 100, 
               seed = set.seed(as.integer(Sys.time()))) {
  
  ## set RNG seed
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
  
  # initialize runtime
  runtime <- list(Screen = NA, Selection = NA, Final_RF = NA)
  
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
  
  # begin runtime for screening step
  start_time <- Sys.time()
  
  ## Screening Step
  screen_result <- screen_RFE(
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
  
  # extract survivor results and initial screen
  survivor_results <- screen_result$survivor_results
  initial_screen <- screen_result$initial_screen
  
  # end runtime for screening step
  end_time <- Sys.time()
  runtime$Screen <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
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
  
  # begin runtime for selection step
  start_time <- Sys.time()
  
  # sets up selection step parameters
  select_args <- list(X_surv, y, num_processors, nodesize)
  select_args <- c(select_args, select_control)
  names(select_args)[1:4] <- c("X", "y", "num_processors", "nodesize")
  
  # RFE via fastshap
  if (shap_model == "full"){
    select_results <- select_RFE(select_args$X, select_args$y, 
                                 select_args$drop_fraction, shap_model,
                                 select_args$number_selected, CLASSIFICATION,
                                 select_args$mtry_factors,
                                 select_args$ntree_factor, select_args$min_ntree,
                                 select_args$num_processors, select_args$nodesize,
                                 nsim = nsim, seed = seed)
  }
  
  # RFE via permutation VIMs
  if (shap_model == "after"){
    select_results <- select_RFE(select_args$X, select_args$y, 
                                 select_args$drop_fraction, shap_model,
                                 select_args$number_selected, CLASSIFICATION,
                                 select_args$mtry_factors,
                                 select_args$ntree_factor, select_args$min_ntree,
                                 select_args$num_processors, select_args$nodesize,
                                 nsim = nsim, seed = seed)
  }
  
  cat("Done. \n")
  
  # end runtime for selection step
  end_time <- Sys.time()
  runtime$Selection <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("Running Final Random Forest...")
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
  
  # begin runtime for final RF
  start_time <- Sys.time()
  
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
    num.threads = num_processors,
    seed = seed
  )
  
  # extracts module membership of final survivors
  final_module_membership <- as.data.frame(cbind(names(X), module_membership),
                                           stringsAsFactors=FALSE)
  names(final_module_membership) <- c("feature_name", "module")
  
  # final variable importance measure via fastshap
  predict_function <- if (CLASSIFICATION) {
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
  
  if(num_processors > 1 ){
    plan(multisession, workers = num_processors)
  }
  
  set.seed(seed)
  
  shap_final_obj <- suppressMessages(fastshap::explain(
    object = final_rf,
    X = final_X,
    pred_wrapper = predict_function,
    nsim = final_nsim,
    adjust = TRUE,
    parallel = parallel,
    shap_only = FALSE, 
    .packages = "ranger"
  ))
  
  if(num_processors > 1 ){
    plan(sequential)
  }
  
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
  
  # end runtime for final RF
  end_time <- Sys.time()
  runtime$Final_RF <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # output function
  out <- shapley_forest(final_rf, final_module_membership,
                        survivor_list=survivor_list,
                        selection_list=selection_list, 
                        final_shap = shap_final_list,
                        shap_obj = shap_final_obj,
                        final_X = final_X,
                        runtime = runtime)
  
  # verbose message
  if (verbose != 0){cat("Done \n")}
  
  # make sure make sure parallel is off now
  if (num_processors > 1) {
    registerDoSEQ()
  }
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
#' @param shap_model        Binary indicator for \code{sf} model. If `full`, \code{sf}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `after`, \code{sf} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `full` is default.
#' @param min_features      Defines minimum feature allowed for each module. If `debug` is `1`
#'                          or `2`, modules below `min_features` after WGCNA call will be
#'                          prompted to the user. See \code{debug} for more info. During
#'                          screening and selection step, if `debug` is 
#'                          not `-1`, modules below `min_features` will only keep non-zero
#'                          important features during each Recursive Feature Elimination
#'                          iteration.
#' @param verbose           Defines the warning message protocol. If `0`, no warning or UI
#'                          will be displayed. If `1`, warnings and UI progress bar will
#'                          be displayed.
#' @param debug             Sets the debugging procedures. If `-1`, all debugging functions
#'                          will be bypassed. If `0`, debugging at the WGCNA will be bypassed.
#'                          If `1`, debugging during Recursive Feature Elimination at both screening 
#'                          and selection step will be bypassed. 
#'                          If `2`, all debugging functions will be ran. Below
#'                          are the debugging features. Debugging at WGCNA detects if each module
#'                          is below the \code{min_features}. Debugging at RFE will keep only 
#'                          non zero important feature at each elimination step for modules below
#'                          \code{min_features}.
#' @param initial           Binary indicator to print out initial screening step results (ie the
#'                          results from the first Recursive Feature Elimination at the screening
#'                          step for each module). If `True`, \code{sf} will pause after RFE
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
#' @param num_processors    Number of processors used to run shapley forests.
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
#' @param seed              RNG seed for reproducibility. Default is based on current
#'                          system time (`as.integer(Sys.time()))`).
#' @return Returns an object of type `shapley_forest`, which is a list containing the essential 
#' output of shapley forests, including a data.frame of selected features and the random forest 
#' model fitted using those features. See \code{shapley_forest} for more details.
#' 
#' @import fuzzyforest
#' @import ranger
#' @import fastshap
#' 
#' 
#' @references
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#' 
#' Daniel Conn, Tuck Ngun, Christina M. Ramirez (2015). Fuzzy Forests: a New
#' WGCNA Based Random Forest Algorithm for Correlated, High-Dimensional Data,
#' Journal of Statistical Software, Manuscript in progress.
#' 
#' Leo Breiman (2001). Random Forests. Machine Learning, 45(1), 5-32.
#'
#' Lundberg, S. M., & Lee, S. I. (2017). A unified approach to interpreting model 
#' predictions. Advances in neural information processing systems, 30.
#'
#' @examples
wsf <- function(X, y, Z=NULL, shap_model = "full",
                    WGCNA_params=WGCNA_control(p=6),
                    min_features=20, verbose = 1, debug = 2, 
                    initial = TRUE, auto_initial = NULL,
                    screen_params=fuzzyforest:::screen_control(min_ntree=5000),
                    select_params=fuzzyforest:::select_control(min_ntree=5000),
                    final_ntree=500, num_processors, nodesize,
                    test_features=NULL, test_y=NULL, nsim=1, final_nsim=100) {
  
  ## validating prerequisites
  if ( !("package:WGCNA" %in% search()) ) {
    stop("WGCNA must be loaded and attached. Type library(WGCNA) to do so.",
         call. = FALSE)
  }
  if (!(is.vector(y) || is.factor(y))) {
    stop("y must be vector or factor")
  }
  
  # Convert X to numeric for WGCNA
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
    stop("shap_model must be `full` or `after`. Type help(sf) or help(wsf) for details.")
  }
  
  if (!verbose %in% c(0, 1)) {
    stop("verbose must be 0 or 1")
  }
  
  if (!debug %in% c(-1, 0, 1, 2)) {
    stop("debug must be -1, 0, 1, or 2")
  }
  
  if (verbose == 0 && !debug %in% c(-1, 1)){
    stop("if debug is 0 or 2, verbose must be 1")
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
    options(warn = -1)
  }
  
  # Set Up for Parameters at each step
  WGCNA_control <- WGCNA_params
  screen_control <- screen_params
  select_control <-  select_params
  
  # Define WGCNA arguments
  WGCNA_args <- list(X,WGCNA_control$power)
  WGCNA_args <- c(WGCNA_args, WGCNA_control$extra_args)
  names(WGCNA_args) <- c("datExpr", "power", names(WGCNA_control$extra_args))
  
  # Run WGCNA dependent on verbose setting
  if (verbose == 0){ # supresses messages
    invisible(capture.output({
      suppressWarnings(suppressMessages({
        bwise <- do.call("blockwiseModules", WGCNA_args)
      }))
    }))
  } else { # allows for WGCNA messages
    bwise <- do.call("blockwiseModules", WGCNA_args)
  }
  
  # Gathers module membership for sf
  module_membership <- bwise$colors
  
  # Debug: if low frequency modules exists, users are warned
  
  min_features <- min_features # defined low frequency threshold
  
  if (!debug %in% c(-1, 0)){ # debug only runs if desired by user
    # finds low frequency modules
    feature_counts <- table(module_membership) 
    low_frequency_modules <- feature_counts[feature_counts <= min_features]
    
    # prompts user if low frequency modules exist
    if (length(low_frequency_modules) > 0) {
      warning(sprintf("\n\n WGCNA - Some modules contain fewer than % s features.", min_features))
      response <- readline(prompt = "Do you wish to continue? (yes/no): ")
      if (tolower(response) != "yes") {
        stop(cat(sprintf("Process terminated by the user. Low Frequency Modules:\n%s", 
                         paste(capture.output(print(low_frequency_modules)), 
                               collapse = "\n")), "\n")) # stops and prints module freq counts
      }
    }
  }
  
  # Sets up for screening parameters
  screen_drop_fraction <- screen_control$drop_fraction
  screen_keep_fraction <- screen_control$keep_fraction
  screen_mtry_factor <- screen_control$mtry_factor
  screen_ntree_factor <- screen_control$ntree_factor
  screen_min_ntree <- screen_control$min_ntree
  
  # Begins Screening Step, calls sf
  #if (verbose != "0"){ cat("Screening Step ... \n")}
  out <- sf(X, y, Z, shap_model, module_membership,
                min_features, verbose, debug, initial, 
                auto_initial, screen_control, 
                select_control, final_ntree,
                num_processors, nodesize=nodesize,
                test_features=test_features, test_y=test_y,
                seed = seed)
  
  # adds WGCNA object at output
  out$WGCNA_object <- bwise
  
  # resets warn function
  options(warn = 1)
  
  return(out)
}
