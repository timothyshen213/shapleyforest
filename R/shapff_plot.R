#' Plot Feature Importance for shapley forest
#'
#' Generates a feature importance plot from `shapley_forest` object. Three types
#' of plots can be generated: bar, beeswarm, and both plots. Plots are created 
#' from \link[shapviz]{shapviz}.
#' 
#' @param object          A `shapley_forest` object.
#' @param kind            The type of plot to generate. Options are "bar", "beeswarm",
#'                        or "both". Default is "beeswarm".
#' @param color_bar_title Title for the color bar in the beeswarm plot. Default is "Feature Value".
#'                        Leave as is if \code{kind} = "bar".
#' @param max_display     The maximum number of features to display.
#'                        If set to `NULL` it is defaulted to maximum possible.
#'                        Default is `NULL`.
#' @param fill            The color to fill bars in the bar plot or color points in the beeswarm plot. 
#'                        Must be a valid \code{\link[ggplot2]{ggplot2}} color or hex code. 
#'                        Default is "#3e568a".
#' @param sort_features   Whether to sort features by importance. Default is `TRUE`.
#' @param show_numbers    Whether to show numeric SHAP values on the bar plot. Default is `TRUE`.
#' @param viridis_args    Additional arguments passed to `viridis` for customizing the color scale 
#'                        in the beeswarm plot. 
#'                        Default is `list(begin = 0.25, end = 0.75, option = "viridis")`.
#' @param ... Obsolete additional arguments.
#' 
#' @return A \code{\link[ggplot2]{ggplot2}} object representing the feature importance plot.
#' 
#' @import ggplot2
#' @import dplyr
#' @import shapviz
#' 
#' @export
plot_importance <- function(object, ...) {
    UseMethod("plot_importance")
}

#' @describeIn plot_importance
#'   Importance plot for an object of class `shapley_forest` through `shapviz`.
#' @export
plot_importance.shapley_forest <- function(object, kind = "beeswarm", 
                                              color_bar_title = "Feature Value",
                                              max_display = NULL, fill = "#3e568a",
                                              sort_features = TRUE,
                                              show_numbers = TRUE, 
                                              viridis_args = list(begin = 0.25, 
                                                                  end = 0.75, 
                                                                  option = "viridis"), ...){
  # stores fastshap object
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)
  
  ## validating prerequisites
  is_valid_color <- function(q) {
    hex_pattern <- "^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$"
    is_hex <- grepl(hex_pattern, q)
    is_ggplot_color <- q %in% colors()
    
    return(is_hex || is_ggplot_color)
  }
  
  if (is.null(max_display)){
    max_display = length(colnames(shap_object))
  }
  
  if (kind != "bar" && kind != "beeswarm" && kind != "both"){
    stop("kind must be `bar`, `beeswarm` or `both`")
  }
  
  if (kind != "bar" && !is.character(color_bar_title)) {
    stop("color_bar_title must be a string")
  }
  if (!is.numeric(max_display)){
    stop("max_display must be numeric")
  }
  if (max_display > length(colnames(shap_object))){
    stop("max_display must be less than total surviving features")
  }
  if (!is.logical(sort_features)){
    stop("sort_features must be boolean")
  }
  if (kind != "bar" && is_valid_color(fill)== FALSE) {
    stop("fill must be a valid ggplot color or hex code")
  }
  if (!is.logical(show_numbers)){
    stop("show_numbers must be boolean")
  }
    
  
  # calls shapviz for importance plot
  importance_plot <- sv_importance(shap_object, kind = kind, color_bar_title = color_bar_title,
                                   max_display = max_display, sort_features = sort_features,
                                   fill = fill, viridis_args = viridis_args)
  
    
  # titles and adjusts importance plot
  importance_plot <- importance_plot + 
      ggtitle("Feature Importance Plot") +
      theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  print(importance_plot)
  
}

#' Plot Waterfall for Specific Observations in shapley forest
#'
#' Generates a waterfall plot for specific observation(s) in a `shapley_forest` object.
#' Plots are  created from \link[shapviz]{shapviz}.
#' 
#' @param object          A `shapley_forest` object.
#' @param row_id          The row index (or indices) corresponding to the observation(s) to 
#'                        be generated on the waterfall plot. Can be character or list. 
#'                        If list, selected observations are averaged.
#' @param row_name        Optional. A name or vector of names for the observations. 
#'                        If not provided, the \code{row_id} will be used.
#' @param max_display     The maximum number of features to display in the plot. 
#'                        If set to `NULL` it is defaulted to maximum possible.
#'                        Default is `NULL`.
#' @param order_fun       Function used to order the features in the waterfall plot. Default orders 
#'                        by absolute SHAP values.
#' @param fill_colors     A vector of colors to fill the bars in the waterfall plot. 
#'                        Each color must be a valid ggplot color or hex code. Default is 
#'                        `c("#59c46b","#3b528b")`.
#' @param contrast        Whether to apply contrast between text and arrow. Default is `TRUE`.
#' @param show_connection Whether to show connections between consecutive bars in the waterfall 
#'                        plot. Default is `TRUE`.
#' @param show_annotation Whether to show annotations for each bar in the waterfall plot. 
#'                        Default is `TRUE`.
#' @param ... Obsolete additional arguments.
#' 
#' @return A \code{\link[ggplot2]{ggplot2}} object representing the waterfall plot.
#'
#' @import ggplot2
#' @import dplyr
#' @import shapviz
#' 
#' @export
plot_waterfall <- function(object, ...) {
  UseMethod("plot_waterfall")
}

#' @describeIn plot_waterfall
#'   Waterfall plot for an object of class `shapley_forest` through `shapviz`.
#' @export

plot_waterfall.shapley_forest <- function(object, row_id, row_name=NULL, max_display=NULL,
                                             order_fun = function(s) order(abs(s)),
                                             fill_colors = c("#59c46b","#3b528b"),
                                             contrast = TRUE, show_connection = TRUE,
                                             show_annotation = TRUE,...){
  
  # stores fastshap object
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)
  
  ## validating prerequisites
  is_valid_color <- function(q) {
    hex_pattern <- "^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$"
    is_hex <- grepl(hex_pattern, q)
    is_ggplot_color <- q %in% colors()
    
    return(is_hex || is_ggplot_color)
  }
  valid_colors <- sapply(fill_colors, is_valid_color)
  
  if (is.null(max_display)){
    max_display = length(colnames(shap_object))
  }
  if (!is.numeric(row_id)){
    stop("row_id must be a valid row index to its corresponding observation")
  }
  if (!is.character(row_name) && !is.null(row_name)){
    stop("row_name must be a character string")
  }
  if (!is.numeric(max_display)){
    stop("max_display must be numeric")
  }
  if (max_display > length(colnames(shap_object))){
    stop("max_display must be less than total surviving features")
  }  
  if (any(!valid_colors)){
    stop("fill_colors must be a valid ggplot color or hex code")
  }
  if (!is.logical(contrast)){
    stop("contrast must be boolean")
  }
  if (!is.logical(show_connection)){
    stop("show_connection must be boolean")
  }
  if (!is.logical(show_annotation)){
    stop("show_annotation must be boolean")
  }
  
  # calls shapviz for waterfall plot
  waterfall_plot <- sv_waterfall(shap_object,row_id,max_display = max_display, 
                                 order_fun = order_fun, fill_colors = fill_colors, 
                                 contrast = contrast, show_connection = show_connection, 
                                 show_annotation = show_annotation)
  
  # titles rows and samples for waterfall plot
  if (is.null(row_name)){
    if (length(row_id) > 1){
      row_name <- paste(row_id, collapse = ", ")
      sample_names <- paste0("Samples ", row_name)
    }
    if (length(row_id) == 1){
      sample_names <- paste0("Sample ", row_id)
    }
  }
  else{
    if (length(row_id) > 1){
      row_name <- paste(row_name, collapse = ", ")
      sample_names <- paste0("Samples ", row_name)
    }
    if (length(row_id) == 1){
      sample_names <- paste0("Sample ", row_name)
    }
  }
  
  # generates title based on sample name
  plot_name <- paste0("Waterfall Plot of ", sample_names)
  waterfall_plot <- waterfall_plot + 
    ggtitle(plot_name) +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  print(waterfall_plot)
}

#' Plot SHAP Force Plot for Specific Observations in shapley forest
#'
#' Generates a force plot for specific observations in a `shapley_forest` object.
#' Plots are created from \link[shapviz]{shapviz}.
#' @export
#' @param object          A `shapley_forest` object.
#' @param row_id          The row index (or indices) corresponding to the observation(s) to 
#'                        be generated on the waterfall plot. Can be character or list. 
#'                        If list, selected observations are averaged.
#' @param row_name        Optional. A name or vector of names for the observations. 
#'                        If not provided, the \code{row_id} will be used.
#' @param max_display     The maximum number of features to display in the plot. 
#'                        If set to `NULL` it is defaulted to maximum possible.
#'                        Default is `NULL`.
#' @param fill_colors     A vector of colors to fill the bars in the force plot. 
#'                        Each color must be a valid ggplot color or hex code. Default is 
#'                        `c("#59c46b","#3b528b")`.
#' @param contrast        Whether to apply contrast between text and arrow. Default is `TRUE`.
#' @param bar_label_size  The size of the labels on the bars. Default is `3.2`.
#' @param show_annotation Whether to show annotations for each bar in the waterfall plot. 
#'                        Default is `TRUE`.
#' @param ... Obsolete additional arguments.
#' 
#' @return A ggplot2 object representing the SHAP force plot.
#' 
#' @export
plot_force <- function(object, ...) {
  UseMethod("plot_force")
}

#' @describeIn plot_force
#'   Force plot for an object of class `shapley_forest` through `shapviz`.
#' @export

plot_force.shapley_forest <- function(object, row_id, row_name=NULL, max_display=NULL,
                                         fill_colors = c("#59c46b","#3b528b"),
                                         contrast = TRUE, bar_label_size = 3.2,
                                         show_annotation = TRUE,...){
  
  # stores fastshap object
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)

  ## validating prerequisites
  is_valid_color <- function(q) {
    hex_pattern <- "^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$"
    is_hex <- grepl(hex_pattern, q)
    is_ggplot_color <- q %in% colors()
    
    return(is_hex || is_ggplot_color)
  }
  valid_colors <- sapply(fill_colors, is_valid_color)
  
  if (is.null(max_display)){
    max_display = length(colnames(shap_object))
  }
  
  if (!is.numeric(row_id)){
    stop("row_id must be a valid row index to its corresponding observation")
  }
  if (!is.character(row_name) && !is.null(row_name)){
    stop("row_name must be a character string")
  }
  if (!is.numeric(max_display)){
    stop("max_display must be numeric")
  }
  if (max_display > length(colnames(shap_object))){
    stop("max_display must be less than total surviving features")
  }  
  if (any(!valid_colors)){
    stop("fill_colors must be a valid ggplot color or hex code")
  }
  if (!is.logical(contrast)){
    stop("contrast must be boolean")
  }
  if (!is.numeric(bar_label_size)){
    stop("bar_label_size must be numerical")
  }
  if (!is.logical(show_annotation)){
    stop("show_annotation must be boolean")
  }
  
  # class shapviz for force plot
  force_plot <- sv_force(shap_object, row_id, max_display = max_display, 
                         fill_colors = fill_colors, contrast = contrast, 
                         bar_label_size = bar_label_size, 
                         show_annotation = show_annotation)
  
  # titles row and sample for force plot
  if (is.null(row_name)){
    if (length(row_id) > 1){
      row_name <- paste(row_id, collapse = ", ")
      sample_names <- paste0("Samples ", row_name)
    }
    if (length(row_id) == 1){
      sample_names <- paste0("Sample ", row_id)
    }
  } else{
    if (length(row_id) > 1){
      row_name <- paste(row_name, collapse = ", ")
      sample_names <- paste0("Samples ", row_name)
    }
    if (length(row_id) == 1){
      sample_names <- paste0("Sample ", row_name)
    }
  }
  
  # generates title based on sample name
  plot_name <- paste0("Force Plot of ", sample_names)
  force_plot <- force_plot + 
    ggtitle(plot_name) +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  print(force_plot)
}

#' Obtains Potential Interaction Matrix in shapley forest.
#'
#' Generates a potential interaction strength matrix of all features.
#' Interactions is calculated from \code{shapviz}'s 
#' \link[shapviz]{potential_interactions}.
#' @export
#' @param object          A `shapley_forest` object.
#' @param ... Obsolete additional arguments.
#' 
#' @return A ggplot2 object representing the SHAP potential interaction matrix.
#' 
#' @export
plot_potential_interactions <- function(object,...){
  UseMethod("plot_potential_interactions")
}

#' @describeIn plot_potential_interactions
#'   Potential Interaction Matrix of class `shapley_forest` through `shapviz`.
#' @export
plot_potential_interactions.shapley_forest <- function(object,...){
  # calculates potential interaction via shapviz
  interaction_data <- detect_interaction(object, 0.1, all=TRUE, verbose=TRUE)
  
  # create interaction matrix
  ggplot(interaction_data, aes(x = FeatureA, y = FeatureB, fill = Interaction_Strength)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Interaction_Strength, 3)), color = "black", size = 3) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = median(interaction_data$Interaction_Strength), 
                         limit = c(min(interaction_data$Interaction_Strength), max(interaction_data$Interaction_Strength)),
                         space = "Lab", name="Interaction Strength") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "FeatureA", y = "FeatureB", title = "Potential Interactions")
  
}
  
#' Plot Decision Plot from shapley forest.
#'
#' Creates a decision plot from a `shapley_forest` object. Plot is adapted from
#' Python's `shap` package.
#' 
#' @export
#' @param object          A `shapley_forest` object.
#' @param highlight       A list of final surviving features to display for 
#'                        decision plot. If `highlight` is `NULL`, all final
#'                        surviving features will be used
#' @param plot_title      Optional. A name for decision plot. Default is `Decision Plot`.
#' @param geom_point      If `TRUE`, it will denote a point at each feature for each 
#'                        observation. Default is `FALSE`.
#' @param gradient        A vector denoting the starting color and ending color for
#'                        the gradient of the predicted outputs. Defaults is `blue`to
#'                        `red`.
#' @param ... Obsolete additional arguments.
#' 
#' @return A ggplot2 object representing the decision plot.
#' 
#' @export
plot_decisions <- function(object,...){
  UseMethod("plot_decisions")
}

#' @describeIn plot_decisions
#'   Decision plot for an object of class `shapley_forest` adapted from Python's `shap`.
#' @export
plot_decisions.shapley_forest <- function(object, highlight = NULL, plot_title = "Decision Plot", 
                                             geom_point = FALSE, 
                                             gradient = c("blue", "red"), ...){

  
  # store fastshap object
  shap_values <- object$shap_obj$shapley_values
  
  #if (object$shap_type == "tree"){
  #  shap_values <- object$shap_obj$shaps
  #}
  
  ## validating prerequisites
  # gathers feature names
  feature_names <- colnames(shap_values)
  
  # checks if highlight parameter contains valid features
  if (!is.null(highlight)){
    if (is.numeric(highlight)){
      if (all(highlight <= 1 & highlight >= nrow(shap_values))){
        stop("highlighted instance id(s) not within range")
      }
    }
    
    if (is.character(highlight)){
      if (!all(highlight %in% feature_name)){
        stop("highlighted instance(s) not valid variable name")
      }
    }
    
    if (!is.character(highlight) && !is.numeric(highlight)){
      stop("highlight must be a character or numeric vector")
    }
  }
  if (!is.character(plot_title)){
    stop("plot_title must be character string. see ggplot2")
  }
  
  if (!is.logical(geom_point)){
    stop("geom_point must be boolean. help(plot_decisions)")
  }
  
  # set up for decision plot
  base_value <- mean(predict(object$final_rf, object$final_X)) # average predicted value
  prediction_out <- predict(object$final_rf, object$final_X) # predictions
  feature_values <- object$final_X # features
  
  # dataframe of shap values
  shap_df <- as.data.frame(shap_values)
  colnames(shap_df) <- feature_names
  shap_df$Observation <- 1:nrow(shap_df)
  # converts to long formate for plotting
  shap_long <- melt(shap_df, id.vars = "Observation", variable.name = "Feature", value.name = "SHAP")
  
  # average SHAP value for each observation
  feature_importance <- colMeans(abs(shap_values))  
  
  # sorts features by average SHAP value
  feature_names <- names(sort(feature_importance, decreasing = FALSE))
  # factorizes feature by order of SHAP value
  shap_long$Feature <- factor(shap_long$Feature, levels = feature_names) 
  
  # groups by observation, calculate cumulative SHAP for each observation
  shap_long <- shap_long %>%
    group_by(Observation) %>%
    arrange(Feature) %>%  
    mutate(Cumulative_SHAP = cumsum(SHAP))  
  
  # starting point for decision plot
  start <- data.frame(
    Observation = 1:nrow(shap_df),
    Feature = "",  
    SHAP = 0,  
    Cumulative_SHAP = 0  
  )
  
  # adds starting point to shap_long
  shap_long <- bind_rows(start, shap_long)
  shap_long$Feature <- factor(shap_long$Feature, levels = c("", feature_names))
  
  # adds end point to shap_long
  last_shap_values <- shap_long %>%
    group_by(Observation) %>%
    summarize(Last_Cumulative_SHAP = last(Cumulative_SHAP))
  shap_long <- shap_long %>%
    left_join(last_shap_values, by = "Observation")
  
  # extracts highlighted features for decision plot
  if (!is.null(highlight)){
    shap_long <- shap_long %>% filter(Observation %in% highlight)  
  }
  
  # generates decision plot
  if (geom_point == FALSE){ # without geom_point
    decision_plot <- ggplot(shap_long, aes(x = Cumulative_SHAP, y = Feature, group = Observation)) +
      geom_path(aes(color = Last_Cumulative_SHAP), size = 1) +  
      geom_vline(xintercept = base_value, color = "#999999", linetype = "dashed") +  
      scale_color_gradient(low = gradient[1], high = gradient[2]) + 
      theme_minimal() +
      labs(title = "Decision Plot",
           x = "Cumulative SHAP Value",
           y = "Feature",
           color = "Predicted Output") +
      theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 10))
  }
  if (geom_point == TRUE){ # with geom_point
    decision_plot <- ggplot(shap_long, aes(x = Cumulative_SHAP, y = Feature, group = Observation)) +
      geom_path(aes(color = Last_Cumulative_SHAP), size = 1) +  
      geom_point(size = 2) + 
      geom_vline(xintercept = base_value, color = "#999999", linetype = "dashed") +  
      scale_color_gradient(low = gradient[1], high = gradient[2]) +  
      theme_minimal() +
      labs(title = "Decision Plot",
           x = "Cumulative SHAP Value",
           y = "Feature",
           color = "Predicted Output") +
      theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 10))
  }
  
  
  print(decision_plot)
}
#' Plot Module Plot from shapley forest
#'
#' Generates a module distribution plot and highlights important features that survived
#' after shapley forest algorithm.
#' The plot function is derived from \link[fuzzyforest]{fuzzyforest}'s modplot function.
#' 
#' @export
#' @param object          A `shapley_forest` object.
#' @param main            Main title for module plot. Default is `NULL` and will title
#'                        "Module Membership Distribution."
#' @param xlab            Title for x-axis. Default is `NULL` and will be titled "Module".
#' @param ylab            Title for y-axis. Default is `NULL` and will be titled
#'                        "Percentage of features in module."
#' @param module_labels   Optional. A list of labels for the modules. Default is `NULL`.
#' @param ... Obsolete additional arguments.
#' 
#' @return A ggplot object representing the module membership distribution colored by importance.
#' 
#' @export
plot_modules <- function(object,...){
  UseMethod("plot_modules")
}
#' @describeIn plot_modules
#'   Module plot for an object of class `shapley_forest` adapted from \code{fuzzyforest}.
#' @export
plot_modules.shapley_forest <- function(object, main=NULL, xlab=NULL, ylab=NULL,
                    module_labels=NULL, ...) {
  # Generates default labels, if applicable
  if(is.null(main)) {
    main <- "Module Membership Distribution"
  }
  if(is.null(xlab)) {
    xlab <- "Module"
  }
  if(is.null(ylab)) {
    ylab <- "Percentage of features in module"
  }
  
  # replaces module names with user defined module_labels
  if(!is.null(module_labels)) {
    old_labels <- object$module_membership$module #old labels
    
    # sorts and factorizes new labels in alphabetical label to match old labels
    module_labels <- module_labels[order(module_labels[, 1]), ] 
    new_labels <- as.character(factor(old_labels, labels=module_labels[, 2])) #
    object$module_membership$module <- new_labels
    
    # adds new labels to module table
    select_mods <- as.factor(object$feature_list$module_membership)
    select_module_table <- module_labels[which(module_labels[, 1] %in%
                                                 levels(select_mods)), ,drop=FALSE]
    
    # handles "." in levels
    if( "." %in% levels(select_mods)) {
      dot_index <- which(levels(select_mods) == ".")
      levels(select_mods)[-dot_index] <- select_module_table[, 2]
    }
    else {
      levels(select_mods) <- select_module_table[, 2]
    }
    object$final_SHAP$module_membership <- as.character(select_mods)
  }
  
  # handles modules for plotting
  mods <- object$module_membership$module
  mod_length <- length(mods)
  mod_tab <- table(mods)
  mod_name <- names(mod_tab)
  final_shap <- object$final_SHAP$module_membership
  
  # handles covariates not in the module
  mod_feature_list <- final_shap[final_shap != "."]
  
  # amount of important feautres in each module
  imp_feature_tab <- table(mod_feature_list)
  imp_names <- names(imp_feature_tab)
  
  # stores count of important features in each module
  feature_tab <- rep(0, length(mod_tab))
  names(feature_tab) <- mod_name
  for(i in 1:length(feature_tab)) {
    if(mod_name[i] %in% names(imp_feature_tab)) {
      feature_tab[i] <- imp_feature_tab[which(imp_names == mod_name[i])]
    }
  }
  
  # calcualtes percent important/unimportant for each module
  unimportant_pct <- (mod_tab - feature_tab)/mod_length
  important_pct <- feature_tab/mod_length
  
  # generates table of the percetnages
  mod_name <- rep(mod_name, 2)
  pct <- c(unimportant_pct, important_pct)
  pct_type <- rep(c("% Unimportant", "% Important"), each=length(mod_tab))
  importance_pct <- data.frame(Module=mod_name, Status=pct_type,
                               Percentage=pct)
  
  # if labels are numeric, it reorders them
  num_mods <- suppressWarnings(as.numeric(object$module_membership[, 2]))
  num_test <- sum(is.na(num_mods))
  if(num_test == 0) {
    importance_pct[, 1] <- as.factor(importance_pct[, 1])
    levels(importance_pct[, 1]) <- sort(unique(num_mods))
  }
  
  # work-around to get rid of notes in R CMD Check
  Module <- NULL
  Percentage <- NULL
  Status <- NULL

  # plots importance plot
  imp_plot <- ggplot(importance_pct, aes(x=Module, y=Percentage, fill=Status)) +
    geom_bar(stat="identity") +
    ggtitle(main) + labs(x = xlab, y = ylab) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),
          legend.title = element_blank())
  
  plot(imp_plot)
}


