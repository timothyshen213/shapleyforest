#' Mossy Forest Object
#'
#' Constructor for the \code{mossy_forest} return type. Users do not call
#' this directly; it is called internally by \code{\link{mf}} and
#' \code{\link{wmf}}.
#'
#' @param final_rf        Final \code{ranger} object fitted on the stable
#'                        feature set.
#' @param final_X         Data frame of features in the stable set.
#' @param module_membership Data frame mapping each feature to its module.
#' @param WGCNA_object    WGCNA output object, or \code{NULL}.
#' @param survivor_list   Per-module screening survivors.
#' @param selection_list  Feature lists at each selection iteration.
#' @param final_shap      Data frame of stable features with mean absolute SHAP,
#'                        95 percent CIs, and shadow-stability frequencies.
#' @param shap_obj        List with elements \code{shapley_values} (signed SHAP
#'                        matrix, samples x stable features) and
#'                        \code{interaction_matrix} (mean absolute TreeSHAP
#'                        interaction matrix, stable features x stable features).
#'                        \code{NULL} when unavailable. For binary classification,
#'                        signed values (and \code{final_shap}'s \code{direction} /
#'                        \code{dir_corr} / \code{signed_importance}) describe the
#'                        probability of \code{levels(y)[1]} — the first level of
#'                        the outcome factor. Relevel \code{y} if you want a
#'                        different reference class.
#' @param runtime         Named list of runtimes (seconds): Screen, Selection,
#'                        Final_RF.
#' @param feature_list    Data frame of the stable feature set with importance
#'                        scores.
#' @param stability_freq  Named numeric vector: shadow-bootstrap selection
#'                        frequency per stable feature.
#' @param stability_data  Data frame with per-pool shadow-stability curves used
#'                        by \code{\link{plot_stability_elbow}}. Columns:
#'                        \code{pool}, \code{rank}, \code{feature}, \code{freq},
#'                        \code{selected}, \code{threshold}, \code{elbow_rank}.
#'                        \code{NULL} if unavailable.
#'
#' @return An object of class \code{mossy_forest}.
#' @keywords internal
mossy_forest <- function(final_rf, final_X, module_membership,
                            WGCNA_object    = NULL,
                            survivor_list,
                            selection_list,
                            final_shap,
                            shap_obj,
                            runtime,
                            feature_list    = NULL,
                            stability_freq  = NULL,
                            stability_data  = NULL) {
  out <- list(
    final_rf         = final_rf,
    final_X          = final_X,
    module_membership = module_membership,
    WGCNA_object     = WGCNA_object,
    survivor_list    = survivor_list,
    selection_list   = selection_list,
    final_SHAP       = final_shap,
    shap_obj         = shap_obj,
    runtimes         = runtime,
    feature_list     = feature_list,
    stability_freq   = stability_freq,
    stability_data   = stability_data
  )
  class(out) <- "mossy_forest"
  out
}


#' Print a Mossy Forest Object
#'
#' Prints the stable feature set with mean |SHAP| importance scores and
#' shadow-stability frequencies.
#'
#' @param x   A \code{mossy_forest} object.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.mossy_forest <- function(x, ...) {
  cat("--- Mossy Forest ---\n")
  cat(sprintf("Stable features selected: %d\n\n", nrow(x$final_SHAP)))
  print(x$final_SHAP)

  if (!is.null(x$final_rf$prediction.error)) {
    cat(sprintf("\nOOB prediction error: %.4f\n", x$final_rf$prediction.error))
  }
  invisible(x)
}


#' Predict from a Mossy Forest Object
#'
#' Generates predictions from the final random forest fitted on the stable
#' feature set.
#'
#' @param object   A \code{mossy_forest} object.
#' @param new_data A matrix or data frame of new observations. Column names
#'                 must match the training features.
#' @param ...      Additional arguments passed to \code{ranger:::predict}.
#'
#' @return A \code{ranger} prediction object.
#' @export
predict.mossy_forest <- function(object, new_data, ...) {
  predict(object$final_rf, data = new_data, ...)
}
