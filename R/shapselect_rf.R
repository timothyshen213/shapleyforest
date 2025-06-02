#' Selection step of shapley forest through SHAP values.
#'
#' Runs the selection step of shapley forest algorithm if `shap_type` = "full".
#' Returns data.frame with the final surviving SHAP values and features.
#' @param X                 A data.frame. With columns denoting a feature vector.
#'                          Could include additional covariates not a part of
#'                          the original modules.
#' @param y                 A Response vector.
#' @param drop_fraction     A value between 0 and 1 representing the percentage of features
#'                          to be dropped at each iteration.
#' @param number_selected   Number of features selected at the end of shapley forest.
#' @param mtry_factor       Number to adjust \code{mtry} for random forest.
#'                          If regression, \code{mtry} is set to
#'                          \code{ceiling}(\eqn{\sqrt(p)}*\code{mtry_factor}).
#'                          If classication, \code{mtry} is set to
#'                          \code{ceiling}((p/3)*\code{mtry_factor}).  If either
#'                          of these numbers is greater than p, \code{mtry} is
#'                          set to p.
#' @param min_ntree         Minimum number of trees grown in each random forest.
#' @param ntree_factor      A number greater than 1 to adjust \code{ntree}.  \code{ntree} for each
#'                          random is \code{ntree_factor} times the number
#'                          of features.  For each random forest, \code{ntree}
#'                          is set to \code{max}(\code{min_ntree},
#'                          \code{ntree_factor}*\code{p}).
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum `nodesize`
#' @param cl                Initialized cluster for parallelization. Set \code{cl} 
#'                          to \code{NULL} if no parallelization is being ran.
#' @param nsim              Number of Monte Carlo repetitions for estimating SHAP
#'                          values in the screening step. Default is `1`. Increasing
#'                          \code{nsim} leads to more accurate results, but at the cost
#'                          of computational cost.
#' @import fuzzyforest
#' @import randomForest
#' @import fastshap
#' 
#' @importFrom randomForest combine
#' 
#' @return A data.frame with final surviving features.

shapselect_RF <- function(X, y, drop_fraction, number_selected, CLASSIFICATION, mtry_factor,
                          ntree_factor, min_ntree,
                          num_processors, nodesize, cl, nsim) {
  # initialize lists
  selection_list <- list()
  feature_list <- NULL
  
  # initializes parallelization
  if(num_processors > 1) {
    snow::stopCluster(cl)
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
    if(num_processors > 1) {
      if (parallel == 1){
        rf = `%dopar%`(foreach(ntree = rep(ntree/num_processors, num_processors)
                               , .combine = randomForest::combine, .packages = 'randomForest'),
                       #second argument to '%dopar%'
                       randomForest(module , y, ntree = ntree, mtry = mtry,
                                    importance = TRUE, scale = FALSE, nodesize=nodesize))
      }
      if (parallel == 1){
        rf <- foreach(ntree = rep(ntree / num_processors, num_processors),
                      .combine = randomForest::combine,
                      .packages = 'randomForest') %dopar% {
                        randomForest(module , y, ntree = ntree, mtry = mtry,
                                     importance = TRUE, scale = FALSE, nodesize = nodesize)
                      }
      }
    }
    if(num_processors == 1) {
      rf <- randomForest(current_X, y, ntree = ntree, mtry = mtry,
                         importance = TRUE, scale = FALSE,
                         nodesize = nodesize)
    }
    ## calculates SHAP values
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
        stop("Invalid or single-class data in y")
      }
      shap <- fastshap::explain(rf, X = current_X, nsim = nsim, 
                                          pred_wrapper = prediction, shap_only = FALSE)
      var_importance <- colMeans(abs(shap$shapley_values))
    }
    # regression case
    if (CLASSIFICATION == FALSE){
      shap <- fastshap::explain(rf, X = current_X, nsim = nsim, 
                                          pred_wrapper = predict, shap_only = FALSE)
      var_importance <- colMeans(abs(shap$shapley_values))
    }
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
    snow::stopCluster(cl)
  }
  
  # output
  out <- list(feature_list, selection_list)
  return(out)
}
