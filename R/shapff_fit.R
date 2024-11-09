#' Fits fuzzy forest algorithm with SHAP values.
#'
#' Fits fuzzy forest algorithm with SHAP values for feature importance.
#' @export
#' @param X                 A data.frame. where each column represents a
#'                          feature vector.
#' @param y                 Response vector. If performing classification, `y` should 
#'                          be a factor or a character. If performing regression, `y`
#'                          should be numeric.
#' @param Z                 A data.frame of additional features that will bypass
#'                          the screening step.
#' @param shap_model        Binary indicator for \code{shapff} model. If `1`, \code{shapff}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `0`, \code{shapff} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `1` is default.
#' @param shap_type         Final SHAP calculation method. If `shapley`, \link[fastshap]{fastshap}
#'                          will be ran. If `tree`, \link[treeshap]{treeshap} will be ran. 
#'                          Running treeshap will have longer runtime, but allows for interaction values
#'                          and plots to be used.
#' @param module_membership A vector that specifies the module membership for each
#'                          each feature.
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
#' @param num_processors    Number of processors used to fit random forests.
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
#'                          X and Z.
#' @param test_y            The responses for the test set.
#' @param nsim              Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the screening step. Default is `1`. Increasing
#'                          \code{nsim} leads to more accurate results, but at the cost
#'                          of computational cost.
#' @param final_nsim        Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the selection step. Default is `1`. \codt{final_nsim}
#'                          should be as large as feasibly possible.
#' @return Returns an object of type fuzzy_forest, which is a list containing the essential 
#' output of fuzzy forests, including a data.frame of selected features and the random forest 
#' model fitted using those features.
#' 
#' @import fuzzyforest
#' @import randomForest
#' @import fastshap
#' 
#' @importFrom randomForest margin
#' @importFrom fastshap explain
#' @importFrom randomForest combine


shapff <- function(X, y, Z=NULL, shap_model = 1, shap_type = "shapley", module_membership,
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
              stop("Invalid or single-class data in y")
          }
          shap <- suppressMessages(fastshap::explain(rf, X = module, nsim = nsim, pred_wrapper = prediction))
        }
        if (CLASSIFICATION == FALSE){
          shap <- suppressMessages(fastshap::explain(rf, X = module, nsim = nsim, pred_wrapper = predict))
        }
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
                                    select_args$number_selected, CLASSIFICATION, select_args$mtry_factors,
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
  if (shap_type =="shapley"){
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
          stop("Invalid or single-class data in y")
        }
        shap_final_obj <- fastshap::explain(final_rf, X = final_X, nsim = final_nsim, 
                                            pred_wrapper = prediction, shap_only = FALSE)
    }
    if (CLASSIFICATION == FALSE){
      shap_final_obj <- fastshap::explain(final_rf, X = final_X, nsim = final_nsim, 
                                          pred_wrapper = predict, shap_only = FALSE)
    }
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
                           final_X = final_X, shap_type = shap_type)
  
  cat("Done \n")
  
  return(out)
}

#' Fits WGCNA fuzzy forest algorithm with SHAP values.
#'
#' Fits fuzzy forest algorithm with WGCNA and SHAP values for feature importance.
#' @export
#' @param X                 A data.frame. where each column represents a
#'                          feature vector. For clustering, WGCNA is employed on `X`.
#'                          Thus, `X` should be numeric. Non-numeric features may
#'                          be input in `Z`.
#' @param y                 Response vector. If performing classification, `y` should 
#'                          be a factor or a character. If performing regression, `y`
#'                          should be numeric.
#' @param Z                 A data.frame of additional features that will bypass
#'                          the screening step, including WGCNA.
#' @param shap_model        Binary indicator for \code{shapff} model. If `1`, \code{shapff}
#'                          runs SHAPley values at both screening and selection step.
#'                          If `0`, \code{shapff} model runs SHAPley values at the end
#'                          of final model and keeps permutation VIMs usage at other steps.
#'                          `1` is default.
#' @param shap_type         Final SHAP calculation method. If `shapley`, \link[fastshap]{fastshap}
#'                          will be ran. If `tree`, \link[treeshap]{treeshap} will be ran. 
#'                          Running treeshap will have longer runtime, but allows for interaction values
#'                          and plots to be used.
#' @param WGCNA_params      WGCNA parameters.
#'                          See \code{\link[WGCNA]{blockwiseModules}} and
#'                          \code{\link[fuzzyforest]{WGCNA_control}} for details.
#'                          \code{WGCNA_params} is an object of type
#'                          \code{WGCNA_control}.
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
#' @param num_processors    Number of processors used to fit random forests.
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
#'                          X and Z.
#' @param test_y            The responses for the test set.
#' @param nsim              Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the screening step. Default is `1`. Increasing
#'                          \code{nsim} leads to more accurate results, but at the cost
#'                          of computational cost.
#' @param final_nsim        Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the selection step. Default is `1`. \codt{final_nsim}
#'                          should be as large as feasibly possible.
#' @return Returns an object of type fuzzy_forest, which is a list containing the essential 
#' output of fuzzy forests, including a data.frame of selected features and the random forest 
#' model fitted using those features.
#' 
#' @import fuzzyforest
#' @import randomForest
#' @import fastshap
#' 
#' @importFrom randomForest margin
#' @importFrom fastshap explain
#' @importFrom randomForest combine
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
#' library(WGCNA)
#' library(randomForest)
#' library(fuzzyforest)
#' data(ctg)
#' y <- ctg$NSP
#' X <- ctg[, 2:22]
#' WGCNA_params <- WGCNA_control(p = 6, minModuleSize = 1, nThreads = 1)
#' mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
#' screen_params <- screen_control(drop_fraction = drop_fraction,
#'                                 keep_fraction = .25, min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' select_params <- select_control(drop_fraction = drop_fraction,
#'                                 number_selected = 5,
#'                                 min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' \donttest{
#' wff_fit <- wff(X, y, WGCNA_params = WGCNA_params,
#'                 screen_params = screen_params,
#'                 select_params = select_params,
#'                 final_ntree = 500)
#'
#' #extract variable importance rankings
#' vims <- wff_fit$feature_list
#'
#' #plot results
#' modplot(wff_fit)
#' }
#' 

shapwff <- function(X, y, Z=NULL, shap_model = 1, shap_type = "shapley", WGCNA_params=WGCNA_control(p=6),
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
  out <- shapff(X, y, Z, shap_model, shap_type, module_membership,
                screen_control, select_control, final_ntree,
                num_processors, nodesize=nodesize,
                test_features=test_features, test_y=test_y)
  out$WGCNA_object <- bwise
  return(out)
}
