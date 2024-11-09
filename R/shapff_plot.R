#' Plot Feature Importance for SHAP Fuzzy Forests
#'
#' Generates a feature importance plot from SHAP Fuzzy Forest object. Three types
#' of plots can be generated: bar, beeswarm, and both plots. Plots are created 
#' from \link[shapviz]{shapviz}.
#' 
#' @param object          A SHAP Fuzzy Forest object.
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
#'   Importance plot for an object of class "shap_fuzzy_forest" through `shapviz`.
#' @export
plot_importance.shap_fuzzy_forest <- function(object, kind = "beeswarm", 
                                              color_bar_title = "Feature Value",
                                              max_display = NULL, fill = "#3e568a",
                                              sort_features = TRUE,
                                              show_numbers = TRUE, 
                                              viridis_args = list(begin = 0.25, 
                                                                  end = 0.75, 
                                                                  option = "viridis"), ...){
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)
  
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
    
  
  importance_plot <- sv_importance(shap_object, kind = kind, color_bar_title = color_bar_title,
                                   max_display = max_display, sort_features = sort_features,
                                   fill = fill, viridis_args = viridis_args)
  
    
  importance_plot <- importance_plot + 
      ggtitle("Feature Importance Plot") +
      theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  print(importance_plot)
  
}

#' Plot Dependence for Selected Features in SHAP Fuzzy Forests
#'
#' Generates a SHAP dependence plot for selected feature(s) in a SHAP Fuzzy Forest object.
#' Scatter plot of SHAP values against the selected feature(s). Plots are created 
#' from \link[shapviz]{shapviz}.
#' 
#' @param object        A SHAP Fuzzy Forest object.
#' @param features      The features for which the SHAP dependence plot should be generated.
#'                      Can be character or list.
#' @param color_var     The variable used to color the points in the plot. Default is "auto", which 
#'                      automatically selects a variable. For no color axis, set to \code{NULL}.
#'                      Can be character or list. See \code{shapviz} \link[shapviz]{sv_dependence()}
#'                      for more details.
#' @param viridis_args    Additional arguments passed to `viridis` for customizing the color scale 
#'                        in the beeswarm plot. 
#'                        Default is `list(begin = 0.25, end = 0.75, option = "viridis")`.
#' @param interaction  To plot interaction values between `features` and `color_var`. Now note, including
#'                     `features` in `color_var`, those corresponding will tree `color_var` = NULL and plot
#'                     the main effect.
#' @return A \code{\link[ggplot2]{ggplot2}} object representing the dependence plot.
#'
#' @import ggplot2
#' @import dplyr
#' @import shapviz
#' 
#' @export
plot_dependence <- function(object, ...) {
  UseMethod("plot_dependence")
}

#' @describeIn plot_dependence
#'   Dependence plot for an object of class "shap_fuzzy_forest" through `shapviz`.
#' @export
plot_dependence.shap_fuzzy_forest <- function(object, features, color_var = "auto",
                                              interaction = FALSE,
                                              viridis_args = list(begin = 0.25, 
                                                                  end = 0.75, 
                                                                  option = "viridis"), ...){
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X, interactions = TRUE)
  
  if (class(shap) != "treeshap"){
    stop("shapley method must be treeshap. set shap_type `tree` in shapff()/shapwff()")
  }
  if (!all(features %in% colnames(shap_object))){
    stop("feature(s) are not in final surviving features")
  }
  
  if (!is.logical(interaction)){
    stop("interaction must be boolean")
  }

  dependence_plot <- sv_dependence(shap_object, v = features, color_var = color_var, 
                                   interactions = interaction, viridis_args = viridis_args)
  print(dependence_plot)
}

#' Plot Waterfall for Specific Observations in SHAP Fuzzy Forests
#'
#' Generates a waterfall plot for specific observation(s) in a SHAP Fuzzy Forest object.
#' Plots are  created from \link[shapviz]{shapviz}.
#' 
#' @param object          A SHAP Fuzzy Forest object.
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
#'   Waterfall plot for an object of class "shap_fuzzy_forest" through `shapviz`.
#' @export

plot_waterfall.shap_fuzzy_forest <- function(object, row_id, row_name=NULL, max_display=NULL,
                                             order_fun = function(s) order(abs(s)),
                                             fill_colors = c("#59c46b","#3b528b"),
                                             contrast = TRUE, show_connection = TRUE,
                                             show_annotation = TRUE,...){
  
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)
  
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
  
  waterfall_plot <- sv_waterfall(shap_object,row_id,max_display = max_display, 
                                 order_fun = order_fun, fill_colors = fill_colors, 
                                 contrast = contrast, show_connection = show_connection, 
                                 show_annotation = show_annotation)
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
  plot_name <- paste0("Waterfall Plot of ", sample_names)
  waterfall_plot <- waterfall_plot + 
    ggtitle(plot_name) +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  print(waterfall_plot)
}

#' Plot SHAP Force Plot for Specific Observations in SHAP Fuzzy Forests
#'
#' Generates a force plot for specific observations in a SHAP Fuzzy Forest object.4
#' Plots are created from \link[shapviz]{shapviz}.
#' @export
#' @param object          A SHAP Fuzzy Forest object.
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
#'                        
#' @return A ggplot2 object representing the SHAP force plot.
#' 
#' @export
plot_force <- function(object, ...) {
  UseMethod("plot_force")
}

#' @describeIn plot_force
#'   Force plot for an object of class "shap_fuzzy_forest" through `shapviz`.
#' @export

plot_force.shap_fuzzy_forest <- function(object, row_id, row_name=NULL, max_display=NULL,
                                         fill_colors = c("#59c46b","#3b528b"),
                                         contrast = TRUE, bar_label_size = 3.2,
                                         show_annotation = TRUE,...){
  
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X)

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
  force_plot <- sv_force(shap_object, row_id, max_display = max_display, 
                         fill_colors = fill_colors, contrast = contrast, 
                         bar_label_size = bar_label_size, 
                         show_annotation = show_annotation)
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
  plot_name <- paste0("Force Plot of ", sample_names)
  force_plot <- force_plot + 
    ggtitle(plot_name) +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  print(force_plot)
}

#' Plot Potential Interaction Matrix in SHAP Fuzzy Forests
#'
#' Generates a potential interaction strength matrix of all features.
#' Interactions is calculated from \code{shapviz}'s 
#' \link[shapviz]{potential_interactions}.
#' @export
#' @param object          A SHAP Fuzzy Forest object.
#'                        
#' @return A ggplot2 object representing the SHAP potential interaction matrix.
#' 
#' @export
plot_potential_interactions <- function(object,...){
  UseMethod("plot_potential_interactions")
}

#' @describeIn plot_potential_interactions
#'   Potential Interaction Matrix of class "shap_fuzzy_forest" through `shapviz`.
#' @export
plot_potential_interactions.shap_fuzzy_forest <- function(object,...){
  interaction_data <- detect_interaction(object, 0.1, all=TRUE, verbose=TRUE)
  
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

#' Plot Interaction Plot from SHAP Fuzzy Forest
#'
#' Generates an interaction plot from SHAP interaction values.
#' Interactions is generated from \code{shapviz}. Note this requires
#' setting \code{shap_type} = `tree` when running \code{shapff} or \code{shapwff}.
#' @export
#' @param object          A SHAP Fuzzy Forest object.
#'                        
#' @return A ggplot2 object representing the SHAP potential interaction matrix.
#' 
#' @export
plot_interactions <- function(object,...){
  UseMethod("plot_interactions")
}

#' @describeIn plot_potential_interactions
#'   Potential Interaction Matrix of class "shap_fuzzy_forest" through `shapviz`.
#' @export
plot_interactions.shap_fuzzy_forest <- function(object, kind = "beeswarm", max_display = "Inf", 
                                                alpha=0.3, bee_width = 0.3, bee_adjust = 0.5,
                                                viridis_args = list(begin = 0.25, 
                                                                    end = 0.75, 
                                                                    option = "viridis"),
                                                color_bar_title = "Row feature value",
                                                sort_features = TRUE, ...){
  shap <- object$shap_obj
  shap_object <- shapviz(shap, X = object$final_X, interactions = TRUE)
  
  if (class(shap) != "treeshap"){
    stop("shapley method must be treeshap. set shap_type `tree` in shapff()/shapwff()")
  }
  if (kind != "beeswarm" && kind != "matrix"){
    stop("kind must be `beeswarm` or `matrix`")
  }
  
  if (!is.character(color_bar_title)) {
    stop("color_bar_title must be a string")
  }
  
  if (!is.numeric(max_display) && max_display != "Inf"){
    stop("max_display must be numeric or `Inf`")
  }
  if (max_display != "Inf" && max_display > length(colnames(shap_object))){
    stop("max_display must be less than total surviving features")
  }
  
  if (!is.logical(sort_features)){
    stop("sort_features must be boolean")
  }
  
  if (kind == "beeswarm"){
    interaction_plot <- sv_interaction(shap_object,kind = "beeswarm", max_display, 
                                       alpha, bee_width, bee_adjust,
                                       viridis_args, color_bar_title,
                                       sort_features) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    plot_name <- paste0("Interaction Plot")
    interaction_plot <- interaction_plot + 
      ggtitle(plot_name) +
      theme(plot.title = element_text(size = 14, hjust = 0.5))
    
    print(interaction_plot)
  }
  
  if (kind == "matrix"){
    interaction_matrix <- sv_interaction(shap_object,kind = "no", max_display, 
                                       sort_features)
    diag(interaction_matrix) <- NA
    interaction_matrix <- melt(interaction_matrix)
    
    max <- max(na.omit(interaction_matrix$value))
    matrix <- ggplot(interaction_matrix, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red",
                           midpoint = 0, limit = c(0, max), space = "Lab", 
                           name="SHAP Interaction Value", na.value = "white") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  
      ggtitle("Interaction Matrix") + 
      xlab(NULL) +  
      ylab(NULL) +  
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.title = element_text(hjust = 0.5))
    print(matrix)
    return(interaction_matrix)
  }
}
  
#' Plot Decision Plot from SHAP Fuzzy Forest
#'
#' TO DO
#' 
#' @export
#' @param object          A SHAP Fuzzy Forest object.
#'                        
#' @return A ggplot2 object representing the SHAP potential interaction matrix.
#' 
#' @export
plot_decisions <- function(object,...){
  UseMethod("plot_decisions")
}

#' @describeIn plot_decisions
#'   TO DO.
#' @export
plot_decisions.shap_fuzzy_forest <- function(object, highlight = NULL, plot_title = "Decision Plot", 
                                             geom_point = FALSE, 
                                             gradient = c("blue", "red"), ...){

  if (object$shap_type == "shapley"){
    shap_values <- object$shap_obj$shapley_values
  }
  
  if (object$shap_type == "tree"){
    shap_values <- object$shap_obj$shaps
  }
  feature_names <- colnames(shap_values)
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
  
  base_value <- mean(predict(object$final_rf, object$final_X))
  prediction_out <- predict(object$final_rf, object$final_X)
  feature_values <- object$final_X
  
  shap_df <- as.data.frame(shap_values)
  colnames(shap_df) <- feature_names
  
  shap_df$Observation <- 1:nrow(shap_df)
  shap_long <- melt(shap_df, id.vars = "Observation", variable.name = "Feature", value.name = "SHAP")
  
  feature_importance <- colMeans(abs(shap_values))  
  feature_names <- names(sort(feature_importance, decreasing = FALSE))
  shap_long$Feature <- factor(shap_long$Feature, levels = feature_names)
  
  shap_long <- shap_long %>%
    group_by(Observation) %>%
    arrange(Feature) %>%  
    mutate(Cumulative_SHAP = cumsum(SHAP))  
  
  start <- data.frame(
    Observation = 1:nrow(shap_df),
    Feature = "",  
    SHAP = 0,  
    Cumulative_SHAP = 0  
  )
  
  shap_long <- bind_rows(start, shap_long)
  shap_long$Feature <- factor(shap_long$Feature, levels = c("", feature_names))
  last_shap_values <- shap_long %>%
    group_by(Observation) %>%
    summarize(Last_Cumulative_SHAP = last(Cumulative_SHAP))
  shap_long <- shap_long %>%
    left_join(last_shap_values, by = "Observation")
  
  if (!is.null(highlight)){
    shap_long <- shap_long %>% filter(Observation %in% highlight)  
  }
  if (geom_point == FALSE){
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
  
  if (geom_point == TRUE){
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

#' Plots relative importance of modules.
#'
#' The plot is designed
#' to depict the size of each module and what percentage of selected
#' features fall into each module.  In particular, it is easy to
#' determine which module is over-represented in the group of selected
#' features.
#' @export
#' @param object   A shap_fuzzy_forest object.
#' @param main Title of plot.
#' @param xlab Title for the x axis.
#' @param ylab Title for the y axis.
#' @param module_labels Labels for the modules.  A data.frame
#'                      or character matrix with first column giving
#'                      the current name of module and second column giving
#'                      the assigned name of each module.
modplot <- function(object, main=NULL, xlab=NULL, ylab=NULL,
                    module_labels=NULL) {
  if(is.null(main)) {
    main <- "Module Membership Distribution"
  }
  if(is.null(xlab)) {
    xlab <- "Module"
  }
  if(is.null(ylab)) {
    ylab <- "Percentage of features in module"
  }
  if(!is.null(module_labels)) {
    old_labels <- object$module_membership[, 2]
    new_labels <- as.factor(old_labels)
    module_labels <- module_labels[order(module_labels[, 1]), ]
    levels(new_labels) <- module_labels[, 2]
    new_labels <- as.character(new_labels)
    object$module_membership[, 2] <- new_labels
    
    select_mods <- as.factor(object$feature_list[, 3])
    select_key <- module_labels[which(module_labels[, 1] %in% levels(select_mods)), ,drop=FALSE]
    if( "." %in% levels(select_mods)) {
      levels(select_mods)[-1] <- select_key[, 2]
    }
    else {
      levels(select_mods) <- select_key[, 2]
    }
    object$feature_list[, 3] <- as.character(select_mods)
  }
  shap_fuzzy_forest <- object
  us_modules <- shap_fuzzy_forest$final_shap$module_membership
  us_modules <- us_modules[us_modules != "."]
  us_modules = as.data.frame(prop.table(table(us_modules))*100)
  us_modules = cbind(us_modules, rep("us", nrow(us_modules)))
  names(us_modules) = c("module", "percent", "type")
  df = as.data.frame(prop.table(table(shap_fuzzy_forest$module_membership[, 2]))*100)
  df = cbind(df, rep("overall", nrow(df)))
  names(df) = c("module", "percent", "type")
  df = rbind(df
             , us_modules
  )
  #check to see if module names are numeric, if so put them in correct order
  num_test <- suppressWarnings(as.numeric(object$module_membership$module))
  if(sum(is.na(num_test))==0) {
    levels(df[,1]) <- as.character(sort(as.numeric(levels(df[,1]))))
  }
  module=5
  percent=5
  type=5
  p_module_dist = ggplot(df
                         , aes(x = module
                               , y = percent
                               , fill = type)
  ) +
    geom_bar(stat = "identity"
             , position="dodge"
             , colour = "#999999"
    ) +
    labs(list(title = main
              , x = xlab
              , y = ylab
    )) +
    theme(axis.line = element_line(colour = "black")
          , panel.grid.major = element_blank()
          , panel.grid.minor = element_blank()
          , panel.border = element_blank()
          , panel.background = element_blank()
          , axis.text.y = element_text(size=10)
          , axis.text.x = element_text(size=10)
          , axis.title = element_text(size=12, face="bold")
          , plot.title = element_text(size=14, face="bold")
    ) +
    scale_fill_manual(values = c("#CDC9C9", "#95C9FF")
                      , name = "Category"
                      , breaks=c("overall", "us")
                      , labels = c("Overall", "Selected Features")
    ) +
    scale_y_continuous(expand=c(0,0))
  plot(p_module_dist)
}


