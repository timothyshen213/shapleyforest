#' Potential interaction variables for fuzzy_forest object
#' 
#' Calculates potential interaction variables from the fuzzy_forest
#' object for a given threshold
#' 
#' @param object  A fuzzy_forest object.
#' @param thresh  Threshold for interaction strength. Threshold should
#'                be between `0` and `1` where the higher values
#'                indicates strong interaction effects.
#' @param ...     Additional arguments not in use.
#' @return A dataframe of potential pairwise interactions
#' 
#' @export
interaction <- function(object, thresh) {
  names <- object$final_SHAP[[1]]
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)
  total_interactions <- list()
  
  for (p in 1:length(names)){
    feat <- names[p]
    interactions <- potential_interactions(shap_object, v = feat)
    sig_interactions <- sapply(interactions, function(x) if (x > thresh) x else NA)
    sig_interactions <- sig_interactions[!is.na(sig_interactions)]
    
    if (length(sig_interactions) > 0) {
      interactions_df <- data.frame(
        FeatureA = feat,
        FeatureB = names(sig_interactions),
        Interaction_Strength = sig_interactions
      )
      
      total_interactions[[feat]] <- interactions_df
    }
  }
  
  if (length(total_interactions) > 0) {
    potential_interactions <- do.call(rbind, total_interactions)
    potential_interactions <- potential_interactions %>%
      mutate(pair = pmap(list(FeatureA, FeatureB), ~ {
        pair <- sort(c(..1, ..2))
        paste(pair, collapse = "-")
      })) %>%
      arrange(pair) %>%
      distinct(pair, .keep_all = TRUE) %>%
      select(-pair)
    
    row.names(potential_interactions) <- NULL
    print(potential_interactions)
  } else {
    print("No significant interactions found.")
  }
  
}

