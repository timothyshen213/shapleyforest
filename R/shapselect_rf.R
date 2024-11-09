#' The selection step of fuzzy forest algorithm
#'
#' Carries out selection step of fuzzyforest algorithm using SHAP values.
#' Returns data.frame with variable importances and top rated features.
#' @param X                 A data.frame.
#'                          Each column corresponds to a feature vectors.
#'                          Could include additional covariates not a part of
#'                          the original modules.
#' @param y                 Response vector.
#' @param drop_fraction     A number between 0 and 1.  Percentage of features
#'                          dropped at each iteration.
#' @param number_selected   Number of features selected by fuzzyforest.
#' @param mtry_factor       In the case of regression, \code{mtry} is set to
#'                          \code{ceiling}(\eqn{\sqrt(p)}*\code{mtry_factor}).
#'                          In the case of classification, \code{mtry} is set to
#'                          \code{ceiling}((p/3)*\code{mtry_factor}).  If either
#'                          of these numbers is greater than p, \code{mtry} is
#'                          set to p.
#' @param min_ntree         Minimum number of trees grown in each random forest.
#' @param ntree_factor      A number greater than 1.  \code{ntree} for each
#'                          random is \code{ntree_factor} times the number
#'                          of features.  For each random forest, \code{ntree}
#'                          is set to \code{max}(\code{min_ntree},
#'                          \code{ntree_factor}*\code{p}).
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum nodesize
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
#' @importFrom randomForest margin
#' @importFrom fastshap explain
#' @importFrom randomForest combine
#' @return A data.frame with the top ranked features.

shapselect_RF <- function(X, y, drop_fraction, number_selected, CLASSIFICATION, mtry_factor,
                          ntree_factor, min_ntree,
                          num_processors, nodesize, cl, nsim) {
  selection_list <- list()
  feature_list <- NULL
  
  if(num_processors > 1) {
    snow::stopCluster(cl)
  }
  num_features <- ncol(X)
  mtry <- min(ceiling(mtry_factor*sqrt(num_features)), dim(X)[2])
  ntree <- max(num_features*ntree_factor, min_ntree)
  target <- number_selected
  current_X <- X
  i <- 1
  while (num_features >= target){
    if(num_processors > 1) {
      rf = `%dopar%`(foreach(ntree = rep(ntree/num_processors, num_processors)
                             , .combine = randomforest::combine, .packages = 'randomForest'),
                     #second argument to '%dopar%'
                     randomForest(current_X , y, ntree = ntree, mtry = mtry,
                                  importance = TRUE, scale = FALSE,
                                  nodesize=nodesize))
    }
    if(num_processors == 1) {
      rf <- randomForest(current_X, y, ntree = ntree, mtry = mtry,
                         importance = TRUE, scale = FALSE,
                         nodesize = nodesize)
    }
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
    if (CLASSIFICATION == FALSE){
      shap <- fastshap::explain(rf, X = current_X, nsim = nsim, 
                                          pred_wrapper = predict, shap_only = FALSE)
      var_importance <- colMeans(abs(shap$shapley_values))
    }
    var_importance <- sort(var_importance, decreasing = TRUE)
    var_importance <- data.frame(Feature = var_importance)
    
    selection_list[[i]] <- data.frame(row.names(var_importance),
                                      round(var_importance[, 1], 4),
                                      stringsAsFactors=FALSE)
    
    names(selection_list[[i]]) <- c("feature_name", "variable_importance")
    i <- i + 1
    reduction <- ceiling(num_features*drop_fraction)
    if(num_features - reduction > target) {
      row_names_var_importance <- row.names(var_importance)
      trimmed_varlist <- var_importance[1:(num_features - reduction), , drop = FALSE]
      features <- row_names_var_importance[1:(num_features - reduction)]
      row.names(trimmed_varlist) <- features
      current_X <- current_X[, which(names(current_X) %in% features)]
      num_features <- length(features)
      mtry <- min(ceiling(mtry_factor*sqrt(num_features)), dim(current_X)[2])
      ntree <- max(num_features*ntree_factor, min_ntree)
    }
    else {
      num_features <- target - 1
      mod_varlist <- var_importance[, 1][1:target]
      features <- row.names(var_importance)[1:target]
      feature_list <- cbind(features, mod_varlist)
      selection_list[[i]] <- as.data.frame(cbind(features, round(mod_varlist, 4)),
                                           stringsAsFactors=FALSE)
      names(selection_list[[i]]) <- c("feature_name", "variable_importance")
    }
  }
  if(num_processors > 1) {
    snow::stopCluster(cl)
  }
  out <- list(feature_list, selection_list)
  return(out)
}
