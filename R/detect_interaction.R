#' Potential interaction variables for `shapley_forest` object
#' 
#' Calculates potential interaction variables from the `shapley_forest`
#' object for a given threshold
#' 
#' @param object  A `shapley_forest` object.
#' @param thresh  Threshold for interaction strength. Threshold should
#'                be between `0` and `1` where the higher values
#'                indicates strong interaction effects.
#' @param all     Returns all interaction. strength values 
#'                Default is `FALSE`.
#' @param verbose Hide potential interaction printout. Defaul is `FALSE`
#' @param ...     Obsolete additional arguements.
#' @return A data.frame of potential pairwise interactions
#' 
#' @export
detect_interaction <- function(object, thresh, all=FALSE, verbose=TRUE) {
  cat("Note: fastshap does not inherently calculate interaction. These are estimates. \n")
  names <- object$final_SHAP[[1]] # feature names
  shap <- object$shap_obj # fastshap obejct
  shap_object <- shapviz(shap, X = object$final_X) # converts to shapviz object
  
  # initialize list
  total_interactions <- list()
  total_sig_interactions <- list()
  
  # for each feature detects significant interaction.
  for (p in 1:length(names)){
    feat <- names[p] # feature
    interactions <- potential_interactions(shap_object, v = feat) # calc potential interaction
    # if above threshold, denote as significant
    sig_interactions <- sapply(interactions, function(x) if (x > thresh) x else NA)
    sig_interactions <- sig_interactions[!is.na(sig_interactions)]
    
    # produces interactions dataframe of significant interactions
    if (length(sig_interactions) > 0) {
      interactions_df <- data.frame(
        FeatureA = feat,
        FeatureB = names(sig_interactions),
        Interaction_Strength = sig_interactions
      )
      
      total_sig_interactions[[feat]] <- interactions_df
    }
    # interaction of all interactions
    interactions_df <- data.frame(
      FeatureA = feat,
      FeatureB = names(interactions),
      Interaction_Strength = interactions
    )
    total_interactions[[feat]] <- interactions_df
  }
  
  # appends to all_interactions
  all_interactions <- do.call(rbind, total_interactions)

  row.names(all_interactions) <- NULL
  
  # produces dataframe of potential interactions
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

