#' Potential interaction variables for fuzzy_forest object
#' 
#' Calculates potential interaction variables from the fuzzy_forest
#' object for a given threshold
#' 
#' @param object  A fuzzy_forest object.
#' @param thresh  Threshold for interaction strength. Threshold should
#'                be between `0` and `1` where the higher values
#'                indicates strong interaction effects.
#' @param all     Returns all interaction. strength values 
#'                Default is `FALSE`.
#' @param verbose Hide potential interaction printout. Defaul is `FALSE`
#' @param ...     Additional arguments not in use.
#' @return A dataframe of potential pairwise interactions
#' 
#' @export
detect_interaction <- function(object, thresh, all=FALSE, verbose=FALSE) {
  cat("Note: fastshap does not inherently calculate interaction. These are estimates. \n")
  names <- object$final_SHAP[[1]]
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)
  total_interactions <- list()
  total_sig_interactions <- list()
  
  
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
      
      total_sig_interactions[[feat]] <- interactions_df
    }
    interactions_df <- data.frame(
      FeatureA = feat,
      FeatureB = names(interactions),
      Interaction_Strength = interactions
    )
    total_interactions[[feat]] <- interactions_df
  }
  all_interactions <- do.call(rbind, total_interactions)

  row.names(all_interactions) <- NULL
  
  if (length(total_sig_interactions) > 0) {
    potential_interactions <- do.call(rbind, total_sig_interactions)
    potential_interactions <- potential_interactions %>%
      mutate(pair = pmap(list(FeatureA, FeatureB), ~ {
        pair <- sort(c(..1, ..2))
        paste(pair, collapse = "-")
      })) %>%
      arrange(pair) %>%
      distinct(pair, .keep_all = TRUE) %>%
      select(-pair)
    
    row.names(potential_interactions) <- NULL
    if (verbose){
      print(potential_interactions)
    }
  } else {
    if (verbose){
      print("No significant interactions found.")
    }
  }
  
  if (all == TRUE){
    return(all_interactions)
  }
}

