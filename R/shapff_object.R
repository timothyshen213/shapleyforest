#' SHAP Fuzzy Forest Object
#'
#' Fuzzy forests with SHAPley values returns an object of type.
#' fuzzyforest.
#' @export
#' @param final_rf          A final random forest fit using the features
#'                          selected by fuzzy forests.
#' @param module_membership Module membership of each feature.
#' @param WGCNA_object      If applicable, output of WGCNA analysis.
#' @param survivor_list     List of features that have survived screening step.
#' @param selection_list    List of features retained at each iteration of
#'                          selection step.
#' @param final_shap        List of survived feature along with its SHAP 
#'                          values
#' @param shap_obj          Final shap object from the final surviving features
#' @return An object of type fuzzy_forest.
shap_fuzzy_forest <- function(final_rf, final_X, module_membership,
                              WGCNA_object=NULL, survivor_list, selection_list, 
                              final_shap, shap_obj, shap_type) {
  out <- list()
  out[[1]] <- final_rf
  out[[2]] <- final_X
  out[[3]] <- module_membership
  out[[4]] <- WGCNA_object
  out[[5]] <- survivor_list
  out[[6]] <- selection_list
  out[[7]] <- final_shap
  out[[8]] <- shap_obj
  out[[9]] <- shap_type
  names(out) <- c("final_rf", "final_X", "module_membership",
                  "WGCNA_object", "survivor_list", "selection_list",
                  "final_SHAP", "shap_obj", "shap_type")
  class(out) <- "shap_fuzzy_forest"
  return(out)
}



#' Prints output from fuzzy forests algorithm.
#' @export
#' @param x   A fuzzy_forest object.
#' @param ... Additional arguments not in use.
#' @return data.frame with list of selected features and variable
#'          importance measures.
#' @note This work was partially funded by NSF IIS 1251151.
print.shap_fuzzy_forest <- function(x, ...) {
  print(x$final_SHAP)
  if(!is.null(x$final_rf$test)) {
    if(!is.null(x$final_rf$test$mse)) {
      cat(c("test set error: ", x$final_rf$test$mse[x$final_rf$ntree]))
    }
    if(!is.null(x$final_rf$test$err.rate)) {
      cat(c("test set error: ", x$final_rf$test$err.rate[x$final_rf$ntree]))
    }
  }
}

#' Predict method for fuzzy_forest object.
#' Obtains predictions from fuzzy forest algorithm.
#' @export
#' @param object   A fuzzy_forest object.
#' @param new_data A matrix or data.frame containing new_data.
#'                 Pay close attention to ensure feature names
#'                 match between training set and test set
#'                 data.frame.
#' @param ...      Additional arguments not in use.
#' @return A vector of predictions
#' @note This work was partially funded by NSF IIS 1251151.
predict.shap_fuzzy_forest <- function(object, new_data, ...) {
  out <- predict(object$final_rf, new_data)
  return(out)
}
