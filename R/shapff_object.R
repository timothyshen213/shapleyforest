#' Shapley Forest loading Message
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Entering the shapley forest ... branching to the Beta Version (06.15.25)--UNDER DEVELOPMENT")
}

#' SHAPley Forest Object
#'
#' Object return type for shapley forest algorithm.
#' 
#' @export
#' @param final_rf          The final random forest object output.
#' @param final_X           The final surviving `X`.
#' @param module_membership Module membership for each feature.
#' @param WGCNA_object      WGCNA object output, if \code{shapwff} was called.
#' @param survivor_list     Feature list of surviving features after screening step.
#' @param selection_list    Feature list of surviving features at each iteration
#'                          of the selection step.
#' @param final_shap        Feature list of final SHAP features.
#' @param shap_obj          Final shapley forest object from the final surviving features
#' @param runtime           Data frame of runtimes of screening step, selection step,
#'                          and final random forest.
#' @return An object of type `shapley_forest`.

shapley_forest <- function(final_rf, final_X, module_membership,
                           WGCNA_object=NULL, survivor_list, selection_list, 
                           final_shap, shap_obj, runtime) {
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
  out[[9]] <- runtime
  
  # define column name for each output
  names(out) <- c("final_rf", "final_X", "module_membership",
                  "WGCNA_object", "survivor_list", "selection_list",
                  "final_SHAP", "shap_obj", "runtimes")
  
  # defines class of object type
  class(out) <- "shapley_forest"
  
  return(out)
}



#' Prints a summarized output of shapley forest.
#' 
#' @export
#' @param x   `shapley_forest` object
#' @param ... Obsolete additional arguments.
#' 
#' @return data.frame with list of final selected features and its SHAP values.
#' 
print.shapley_forest <- function(x, ...) {
  print(x$final_SHAP)
  
  # prints test error, if applicable
  if(!is.null(x$final_rf$test)) {
    if(!is.null(x$final_rf$test$mse)) {
      cat(c("test set error: ", x$final_rf$test$mse[x$final_rf$ntree]))
    }
    if(!is.null(x$final_rf$test$err.rate)) {
      cat(c("test set error: ", x$final_rf$test$err.rate[x$final_rf$ntree]))
    }
  }
}

#' Prediction for shapley forest
#' 
#' Prediction function for shapley forest algorithm.
#' Runs prediction through the final random forest
#' after both screening and selection steps.
#' 
#' @export
#' @param object   A `shapley_forest` object.
#' @param new_data A matrix or data.frame containing new data.
#'                 Note: Feature names must match between 
#'                 training set and test set data.frame.
#' @param ...      Obsolete additional arguments.
#' 
#' @return A vector of predictions through random forest.
#' 
#' @examples
#' # see help(shapff) for example

predict.shapley_forest <- function(object, new_data, ...) {
  # obtains prediction through final_rf model
  out <- predict(object$final_rf, new_data)
  return(out)
}

