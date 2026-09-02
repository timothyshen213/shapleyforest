#' Plot Feature Importance for a Mossy Forest Object
#'
#' Generates a feature importance plot from a \code{mossy_forest} object.
#' Three plot types are available: bar, beeswarm, and both. Plots are built
#' via \link[shapviz]{shapviz}.
#'
#' @param object          A \code{mossy_forest} object.
#' @param kind            Plot type: \code{"bar"}, \code{"beeswarm"}, or
#'                        \code{"both"}. Default \code{"beeswarm"}.
#' @param color_bar_title Title for the colour bar in the beeswarm plot.
#'                        Default \code{"Feature Value"}.
#' @param max_display     Maximum number of features to display. \code{NULL}
#'                        shows all. Default \code{NULL}.
#' @param fill            Bar/point colour. Any valid \code{ggplot2} colour or
#'                        hex code. Default \code{"#3e568a"}.
#' @param sort_features   Sort features by importance? Default \code{TRUE}.
#' @param show_numbers    Show numeric SHAP values on the bar plot?
#'                        Default \code{TRUE}.
#' @param viridis_args    Named list of arguments forwarded to \code{viridis}
#'                        for the beeswarm colour scale. Default
#'                        \code{list(begin = 0.25, end = 0.75, option = "viridis")}.
#' @param ...             Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @import shapviz
#' @export
plot_importance <- function(object, ...) {
  UseMethod("plot_importance")
}

#' @describeIn plot_importance Importance plot for a \code{mossy_forest} object.
#' @export
plot_importance.mossy_forest <- function(object,
                                           kind            = "beeswarm",
                                           color_bar_title = "Feature Value",
                                           max_display     = NULL,
                                           fill            = "#3e568a",
                                           sort_features   = TRUE,
                                           show_numbers    = TRUE,
                                           viridis_args    = list(begin  = 0.25,
                                                                   end    = 0.75,
                                                                   option = "viridis"),
                                           ...) {

  shap_object <- shapviz::shapviz(object$shap_obj$shapley_values,
                                   X = object$final_X)

  is_valid_color <- function(x) {
    grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", x) || x %in% colors()
  }

  if (is.null(max_display)) max_display <- length(colnames(shap_object))

  if (!kind %in% c("bar", "beeswarm", "both"))
    stop("kind must be 'bar', 'beeswarm', or 'both'")
  if (kind != "bar" && !is.character(color_bar_title))
    stop("color_bar_title must be a string")
  if (!is.numeric(max_display))
    stop("max_display must be numeric")
  if (max_display > length(colnames(shap_object)))
    stop("max_display must not exceed the number of stable features")
  if (!is.logical(sort_features))
    stop("sort_features must be logical")
  if (kind != "bar" && !is_valid_color(fill))
    stop("fill must be a valid ggplot2 colour or hex code")
  if (!is.logical(show_numbers))
    stop("show_numbers must be logical")

  p <- shapviz::sv_importance(shap_object,
                               kind            = kind,
                               color_bar_title = color_bar_title,
                               max_display     = max_display,
                               sort_features   = sort_features,
                               fill            = fill,
                               viridis_args    = viridis_args)

  p <- p +
    ggtitle("Feature Importance") +
    theme(plot.title = element_text(size = 14, hjust = 0.5))

  print(p)
  invisible(p)
}


#' Plot Waterfall for Specific Observations in a Mossy Forest Object
#'
#' Generates a SHAP waterfall plot for one or more observations. Plots are
#' built via \link[shapviz]{shapviz}.
#'
#' @param object        A \code{mossy_forest} object.
#' @param row_id        Row index (or indices) for the observation(s) to plot.
#'                      When multiple indices are given, SHAP values are averaged.
#' @param row_name      Optional name(s) for the observation(s).
#' @param max_display   Maximum number of features to display. \code{NULL}
#'                      shows all. Default \code{NULL}.
#' @param order_fun     Function used to order features. Default orders by
#'                      absolute SHAP value.
#' @param fill_colors   Length-2 vector of colours for positive/negative bars.
#'                      Default \code{c("#59c46b", "#3b528b")}.
#' @param contrast      Apply text/arrow contrast? Default \code{TRUE}.
#' @param show_connection Show connections between bars? Default \code{TRUE}.
#' @param show_annotation Annotate bars? Default \code{TRUE}.
#' @param ...           Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @import shapviz
#' @export
plot_waterfall <- function(object, ...) {
  UseMethod("plot_waterfall")
}

#' @describeIn plot_waterfall Waterfall plot for a \code{mossy_forest} object.
#' @export
plot_waterfall.mossy_forest <- function(object,
                                          row_id,
                                          row_name        = NULL,
                                          max_display     = NULL,
                                          order_fun       = function(s) order(abs(s)),
                                          fill_colors     = c("#59c46b", "#3b528b"),
                                          contrast        = TRUE,
                                          show_connection = TRUE,
                                          show_annotation = TRUE,
                                          ...) {

  shap_object <- shapviz::shapviz(object$shap_obj$shapley_values,
                                   X = object$final_X)

  is_valid_color <- function(x) {
    grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", x) || x %in% colors()
  }

  if (is.null(max_display)) max_display <- length(colnames(shap_object))
  if (!is.numeric(row_id))
    stop("row_id must be a valid row index")
  if (!is.null(row_name) && !is.character(row_name))
    stop("row_name must be a character string")
  if (!is.numeric(max_display))
    stop("max_display must be numeric")
  if (max_display > length(colnames(shap_object)))
    stop("max_display must not exceed the number of stable features")
  if (any(!sapply(fill_colors, is_valid_color)))
    stop("fill_colors must be valid ggplot2 colours or hex codes")
  if (!is.logical(contrast))
    stop("contrast must be logical")
  if (!is.logical(show_connection))
    stop("show_connection must be logical")
  if (!is.logical(show_annotation))
    stop("show_annotation must be logical")

  p <- shapviz::sv_waterfall(shap_object, row_id,
                              max_display     = max_display,
                              order_fun       = order_fun,
                              fill_colors     = fill_colors,
                              contrast        = contrast,
                              show_connection = show_connection,
                              show_annotation = show_annotation)

  if (is.null(row_name)) {
    label <- if (length(row_id) > 1) paste0("Samples ", paste(row_id, collapse = ", "))
             else                    paste0("Sample ", row_id)
  } else {
    label <- if (length(row_id) > 1) paste0("Samples ", paste(row_name, collapse = ", "))
             else                    paste0("Sample ", row_name)
  }

  p <- p +
    ggtitle(paste0("Waterfall Plot — ", label)) +
    theme(plot.title = element_text(size = 14, hjust = 0.5))

  print(p)
  invisible(p)
}


#' Plot SHAP Force Plot for Specific Observations in a Mossy Forest Object
#'
#' Generates a SHAP force plot for one or more observations. Plots are built
#' via \link[shapviz]{shapviz}.
#'
#' @param object          A \code{mossy_forest} object.
#' @param row_id          Row index (or indices) for the observation(s) to plot.
#' @param row_name        Optional name(s) for the observation(s).
#' @param max_display     Maximum number of features to display. \code{NULL}
#'                        shows all. Default \code{NULL}.
#' @param fill_colors     Length-2 vector of colours. Default
#'                        \code{c("#59c46b", "#3b528b")}.
#' @param contrast        Apply text/arrow contrast? Default \code{TRUE}.
#' @param bar_label_size  Label size on bars. Default \code{3.2}.
#' @param show_annotation Annotate bars? Default \code{TRUE}.
#' @param ...             Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @import shapviz
#' @export
plot_force <- function(object, ...) {
  UseMethod("plot_force")
}

#' @describeIn plot_force Force plot for a \code{mossy_forest} object.
#' @export
plot_force.mossy_forest <- function(object,
                                      row_id,
                                      row_name        = NULL,
                                      max_display     = NULL,
                                      fill_colors     = c("#59c46b", "#3b528b"),
                                      contrast        = TRUE,
                                      bar_label_size  = 3.2,
                                      show_annotation = TRUE,
                                      ...) {

  shap_object <- shapviz::shapviz(object$shap_obj$shapley_values,
                                   X = object$final_X)

  is_valid_color <- function(x) {
    grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", x) || x %in% colors()
  }

  if (is.null(max_display)) max_display <- length(colnames(shap_object))
  if (!is.numeric(row_id))
    stop("row_id must be a valid row index")
  if (!is.null(row_name) && !is.character(row_name))
    stop("row_name must be a character string")
  if (!is.numeric(max_display))
    stop("max_display must be numeric")
  if (max_display > length(colnames(shap_object)))
    stop("max_display must not exceed the number of stable features")
  if (any(!sapply(fill_colors, is_valid_color)))
    stop("fill_colors must be valid ggplot2 colours or hex codes")
  if (!is.logical(contrast))
    stop("contrast must be logical")
  if (!is.numeric(bar_label_size))
    stop("bar_label_size must be numeric")
  if (!is.logical(show_annotation))
    stop("show_annotation must be logical")

  p <- shapviz::sv_force(shap_object, row_id,
                          max_display     = max_display,
                          fill_colors     = fill_colors,
                          contrast        = contrast,
                          bar_label_size  = bar_label_size,
                          show_annotation = show_annotation)

  if (is.null(row_name)) {
    label <- if (length(row_id) > 1) paste0("Samples ", paste(row_id, collapse = ", "))
             else                    paste0("Sample ", row_id)
  } else {
    label <- if (length(row_id) > 1) paste0("Samples ", paste(row_name, collapse = ", "))
             else                    paste0("Sample ", row_name)
  }

  p <- p +
    ggtitle(paste0("Force Plot — ", label)) +
    theme(plot.title = element_text(size = 14, hjust = 0.5))

  print(p)
  invisible(p)
}


#' Plot Shadow-Stability Frequencies for a Mossy Forest Object
#'
#' Displays the bootstrap shadow-stability selection frequency for each stable
#' feature as a horizontal bar chart. Features are sorted by frequency.
#'
#' @param object  A \code{mossy_forest} object.
#' @param fill    Bar colour. Default \code{"#3e568a"}.
#' @param main    Plot title. Default \code{"Shadow-Stability Frequencies"}.
#' @param ...     Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @export
plot_stability <- function(object, ...) {
  UseMethod("plot_stability")
}

#' @describeIn plot_stability Stability frequency plot for a \code{mossy_forest} object.
#' @export
plot_stability.mossy_forest <- function(object,
                                          fill = "#3e568a",
                                          main = "Shadow-Stability Frequencies",
                                          ...) {

  freq <- object$stability_freq
  if (is.null(freq) || length(freq) == 0)
    stop("No stability frequencies found in this mossy_forest object.")

  df <- data.frame(
    feature   = names(freq),
    frequency = as.numeric(freq),
    stringsAsFactors = FALSE
  )
  df$feature <- factor(df$feature,
                        levels = df$feature[order(df$frequency)])

  p <- ggplot(df, aes(x = frequency, y = feature)) +
    geom_col(fill = fill) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    labs(x = "Selection frequency", y = NULL, title = main) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 14, hjust = 0.5))

  print(p)
  invisible(p)
}


#' Plot Shadow-Stability Selection Curves for a Mossy Forest Object
#'
#' Displays the shadow-bootstrap stability frequency curve for each selection
#' pool with a horizontal threshold line and coloured points indicating which
#' features were selected.
#'
#' For \code{shadow_mode = "split"} two panels are shown: one for the
#' correlated pool (fixed \eqn{\pi} threshold) and one for the unassigned
#' pool (elbow or fixed threshold).  For \code{shadow_mode = "within_module"}
#' one panel is shown per module.
#'
#' The dashed red line marks the stability threshold applied. For unassigned
#' pools with elbow detection the elbow rank is marked with a vertical dotted
#' line.
#'
#' @param object         A \code{mossy_forest} object.
#' @param col_selected   Point colour for selected features.
#'                       Default \code{"#3e568a"}.
#' @param col_unselected Point colour for unselected features.
#'                       Default \code{"gray65"}.
#' @param col_threshold  Colour of the horizontal threshold line.
#'                       Default \code{"tomato"}.
#' @param col_elbow      Colour of the vertical elbow-rank line.
#'                       Default \code{"#e07b00"}.
#' @param show_labels    Show feature names on points? Default \code{FALSE}.
#'                       Set \code{TRUE} only for small pools (\eqn{p \le 20}).
#' @param ncol           Number of columns in the facet grid. \code{NULL}
#'                       lets ggplot choose. Default \code{NULL}.
#' @param main           Plot title. Default \code{"Shadow-stability selection"}.
#' @param ...            Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @export
plot_stability_elbow <- function(object, ...) {
  UseMethod("plot_stability_elbow")
}

#' @describeIn plot_stability_elbow Stability elbow plot for a
#'   \code{mossy_forest} object.
#' @export
plot_stability_elbow.mossy_forest <- function(object,
                                                 col_selected   = "#3e568a",
                                                 col_unselected = "gray65",
                                                 col_threshold  = "tomato",
                                                 col_elbow      = "#e07b00",
                                                 show_labels    = FALSE,
                                                 ncol           = NULL,
                                                 main           = "Shadow-stability selection",
                                                 ...) {

  sd <- object$stability_data
  if (is.null(sd) || nrow(sd) == 0L)
    stop(paste0(
      "No stability selection data found in this mossy_forest object.\n",
      "  Refit with mf() to populate stability_data."
    ), call. = FALSE)

  # ── ensure factor order: pools in the order they appear ─────────────────────
  pool_order    <- unique(sd$pool)
  sd$pool_fct   <- factor(sd$pool, levels = pool_order)
  sd$selected_f <- factor(sd$selected, levels = c(TRUE, FALSE),
                           labels = c("Selected", "Not selected"))

  # ── one threshold-line value per pool (same for all rows in that pool) ───────
  thr_df <- unique(sd[, c("pool_fct", "threshold")])

  # ── elbow-rank lines (NA pools are skipped) ───────────────────────────────────
  elbow_df <- unique(sd[!is.na(sd$elbow_rank), c("pool_fct", "elbow_rank")])

  # ── base plot ─────────────────────────────────────────────────────────────────
  p <- ggplot(sd, aes(x = rank, y = freq)) +
    geom_line(colour = "gray70", linewidth = 0.5) +
    geom_point(aes(colour = selected_f), size = 2.5, alpha = 0.9) +
    # threshold line per facet
    geom_hline(data = thr_df,
               aes(yintercept = threshold),
               colour   = col_threshold,
               linetype = "dashed",
               linewidth = 0.7) +
    scale_colour_manual(
      values = c("Selected" = col_selected, "Not selected" = col_unselected),
      name   = NULL
    ) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1L)) +
    facet_wrap(~ pool_fct, scales = "free_x", ncol = ncol) +
    labs(
      x     = "Feature rank (stability freq descending)",
      y     = "Stability frequency",
      title = main
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text      = element_text(face = "bold", size = 11),
      plot.title      = element_text(size = 14, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  # ── elbow rank: vertical dotted line + annotation ────────────────────────────
  if (nrow(elbow_df) > 0L) {
    p <- p +
      geom_vline(data = elbow_df,
                 aes(xintercept = elbow_rank),
                 colour   = col_elbow,
                 linetype = "dotted",
                 linewidth = 0.8) +
      geom_label(data = elbow_df,
                 aes(x = elbow_rank, y = 0.97,
                     label = paste0("elbow\nrank ", elbow_rank)),
                 colour    = col_elbow,
                 size      = 2.8,
                 hjust     = 0,
                 label.size = 0,
                 fill      = "white",
                 alpha     = 0.85,
                 inherit.aes = FALSE)
  }

  # ── optional feature labels ───────────────────────────────────────────────────
  if (show_labels) {
    p <- p +
      geom_text(aes(label = feature),
                size  = 2.5,
                hjust = -0.1,
                vjust = 0.5,
                check_overlap = TRUE)
  }

  print(p)
  invisible(p)
}


#' Plot Potential Interaction Matrix for a Mossy Forest Object
#'
#' Generates a heatmap of pairwise SHAP-based potential interaction strengths
#' for all stable features.
#'
#' @param object A \code{mossy_forest} object.
#' @param ...    Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @export
plot_potential_interactions <- function(object, ...) {
  UseMethod("plot_potential_interactions")
}

#' @describeIn plot_potential_interactions Interaction heatmap for a
#'   \code{mossy_forest} object (uses exact TreeSHAP interaction values).
#' @export
plot_potential_interactions.mossy_forest <- function(object, ...) {

  int_mat <- object$shap_obj$interaction_matrix
  if (is.null(int_mat) || !is.matrix(int_mat))
    stop(paste0(
      "No interaction matrix found. ",
      "Refit with mf() to compute TreeSHAP interaction values."
    ), call. = FALSE)

  feats <- rownames(int_mat)
  p_n   <- length(feats)

  # upper-triangle long format (off-diagonal pairs only)
  rows <- vector("list", p_n * (p_n - 1L) / 2L)
  k    <- 1L
  for (j in seq_len(p_n - 1L)) {
    for (l in (j + 1L):p_n) {
      rows[[k]] <- data.frame(
        FeatureA             = feats[j],
        FeatureB             = feats[l],
        Interaction_Strength = int_mat[j, l],
        stringsAsFactors     = FALSE
      )
      k <- k + 1L
    }
  }
  plot_df <- do.call(rbind, rows)

  # symmetric: add the mirror so heatmap is square
  mirror  <- data.frame(
    FeatureA             = plot_df$FeatureB,
    FeatureB             = plot_df$FeatureA,
    Interaction_Strength = plot_df$Interaction_Strength,
    stringsAsFactors     = FALSE
  )
  # diagonal (main effects)
  diag_df <- data.frame(
    FeatureA             = feats,
    FeatureB             = feats,
    Interaction_Strength = diag(int_mat),
    stringsAsFactors     = FALSE
  )
  full_df <- rbind(plot_df, mirror, diag_df)
  full_df$FeatureA <- factor(full_df$FeatureA, levels = feats)
  full_df$FeatureB <- factor(full_df$FeatureB, levels = rev(feats))

  mid_val  <- median(plot_df$Interaction_Strength)
  lim_vals <- range(full_df$Interaction_Strength)

  p <- ggplot(full_df,
              aes(x = FeatureA, y = FeatureB,
                  fill = Interaction_Strength)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.3f", Interaction_Strength)),
              colour = "black", size = 2.5) +
    scale_fill_gradient2(
      low      = "#4575b4",
      mid      = "white",
      high     = "#d73027",
      midpoint = mid_val,
      limits   = lim_vals,
      name     = "Mean |SHAP\ninteraction|"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1),
      plot.title   = element_text(size = 14, hjust = 0.5)
    ) +
    labs(x = NULL, y = NULL,
         title = "TreeSHAP Interaction Values",
         subtitle = "Off-diagonal: pairwise interactions  |  Diagonal: main effects")

  print(p)
  invisible(p)
}


#' Plot Decision Plot for a Mossy Forest Object
#'
#' Creates a cumulative SHAP decision plot adapted from the Python \code{shap}
#' package.
#'
#' @param object     A \code{mossy_forest} object. If \code{just_shap = TRUE},
#'                   pass a data frame of SHAP values instead.
#' @param highlight  Row indices or feature names to highlight. \code{NULL}
#'                   shows all observations. Default \code{NULL}.
#' @param just_shap  If \code{TRUE}, \code{object} is treated as a plain SHAP
#'                   data frame. Default \code{FALSE}.
#' @param plot_title Plot title. Default \code{"Decision Plot"}.
#' @param geom_point Add a point at each feature step? Default \code{FALSE}.
#' @param gradient   Length-2 colour vector for the prediction gradient.
#'                   Default \code{c("blue", "red")}.
#' @param ...        Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot_decisions <- function(object, ...) {
  UseMethod("plot_decisions")
}

#' @describeIn plot_decisions Decision plot for a \code{mossy_forest} object.
#' @export
plot_decisions.mossy_forest <- function(object,
                                          highlight  = NULL,
                                          just_shap  = FALSE,
                                          plot_title = "Decision Plot",
                                          geom_point = FALSE,
                                          gradient   = c("blue", "red"),
                                          ...) {

  if (!just_shap) {
    shap_values <- object$shap_obj$shapley_values
  } else {
    if (!is.data.frame(object))
      stop("object is not a data frame. Should just_shap be FALSE?")
    shap_values <- object
  }

  if (!is.character(plot_title))
    stop("plot_title must be a character string")
  if (!is.logical(geom_point))
    stop("geom_point must be logical")

  feature_names <- colnames(shap_values)

  if (!is.null(highlight)) {
    if (!is.character(highlight) && !is.numeric(highlight))
      stop("highlight must be a character or numeric vector")
  }

  # base value: mean prediction on training data (not available in just_shap mode)
  base_value <- if (!just_shap)
    mean(predict(object$final_rf, data = object$final_X)$predictions)
  else
    NULL

  shap_df           <- as.data.frame(shap_values)
  colnames(shap_df) <- feature_names
  shap_df$Observation <- seq_len(nrow(shap_df))

  # long format
  shap_long <- reshape(shap_df,
                        varying   = feature_names,
                        v.names   = "SHAP",
                        timevar   = "Feature",
                        times     = feature_names,
                        direction = "long")
  shap_long$Feature <- factor(shap_long$Feature)
  shap_long$id      <- NULL

  # sort features by mean |SHAP|
  feat_order <- names(sort(colMeans(abs(shap_values)), decreasing = FALSE))
  shap_long$Feature <- factor(shap_long$Feature, levels = feat_order)

  # cumulative SHAP per observation
  shap_long <- shap_long %>%
    dplyr::group_by(Observation) %>%
    dplyr::arrange(Feature) %>%
    dplyr::mutate(Cumulative_SHAP = cumsum(SHAP)) %>%
    dplyr::ungroup()

  # starting baseline row
  start <- data.frame(
    Observation      = seq_len(nrow(shap_df)),
    Feature          = factor("", levels = c("", feat_order)),
    SHAP             = 0,
    Cumulative_SHAP  = 0
  )
  shap_long$Feature <- factor(as.character(shap_long$Feature),
                               levels = c("", feat_order))
  shap_long <- dplyr::bind_rows(start, shap_long)

  # terminal cumulative value for colour scale
  last_vals <- shap_long %>%
    dplyr::group_by(Observation) %>%
    dplyr::summarise(Last_Cumulative_SHAP = dplyr::last(Cumulative_SHAP),
                     .groups = "drop")
  shap_long <- dplyr::left_join(shap_long, last_vals, by = "Observation")

  if (!is.null(highlight))
    shap_long <- dplyr::filter(shap_long, Observation %in% highlight)

  p <- ggplot(shap_long,
              aes(x = Cumulative_SHAP, y = Feature,
                  group = Observation)) +
    geom_path(aes(colour = Last_Cumulative_SHAP), linewidth = 0.8) +
    scale_colour_gradient(low = gradient[1], high = gradient[2]) +
    theme_minimal() +
    labs(title  = plot_title,
         x      = "Cumulative SHAP Value",
         y      = "Feature",
         colour = "Predicted\nOutput") +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          plot.title  = element_text(size = 14, hjust = 0.5))

  # add baseline vline only when the model is available
  if (!is.null(base_value))
    p <- p + geom_vline(xintercept = base_value,
                         colour = "#999999", linetype = "dashed")

  if (geom_point)
    p <- p + geom_point(size = 1.5)

  print(p)
  invisible(p)
}


#' Plot Module Membership Distribution for a Mossy Forest Object
#'
#' Shows the proportion of features per module that were selected as stable,
#' coloured by selection status.
#'
#' @param object        A \code{mossy_forest} object.
#' @param main          Plot title. Default \code{"Module Membership Distribution"}.
#' @param xlab          X-axis label. Default \code{"Module"}.
#' @param ylab          Y-axis label. Default \code{"Proportion of features"}.
#' @param module_labels Optional two-column data frame mapping old module labels
#'                      to new display labels. Default \code{NULL}.
#' @param ...           Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @export
plot_modules <- function(object, ...) {
  UseMethod("plot_modules")
}

#' @describeIn plot_modules Module membership plot for a \code{mossy_forest} object.
#' @export
plot_modules.mossy_forest <- function(object,
                                        main          = NULL,
                                        xlab          = NULL,
                                        ylab          = NULL,
                                        module_labels = NULL,
                                        ...) {

  if (is.null(main)) main <- "Module Membership Distribution"
  if (is.null(xlab)) xlab <- "Module"
  if (is.null(ylab)) ylab <- "Proportion of features"

  if (!is.null(module_labels)) {
    old_labels <- object$module_membership$module
    module_labels <- module_labels[order(module_labels[, 1]), ]
    new_labels <- as.character(factor(old_labels, labels = module_labels[, 2]))
    object$module_membership$module <- new_labels

    select_mods <- as.factor(object$feature_list$module_membership)
    select_module_table <- module_labels[
      module_labels[, 1] %in% levels(select_mods), , drop = FALSE]

    if ("." %in% levels(select_mods)) {
      dot_idx <- which(levels(select_mods) == ".")
      levels(select_mods)[-dot_idx] <- select_module_table[, 2]
    } else {
      levels(select_mods) <- select_module_table[, 2]
    }
    object$final_SHAP$module_membership <- as.character(select_mods)
  }

  mods       <- object$module_membership$module
  mod_length <- length(mods)
  mod_tab    <- table(mods)
  mod_name   <- names(mod_tab)

  final_shap      <- object$final_SHAP$module_membership
  mod_feature_list <- final_shap[final_shap != "."]
  imp_feature_tab  <- table(mod_feature_list)
  imp_names        <- names(imp_feature_tab)

  feature_tab <- rep(0, length(mod_tab))
  names(feature_tab) <- mod_name
  for (i in seq_along(feature_tab)) {
    if (mod_name[i] %in% imp_names)
      feature_tab[i] <- imp_feature_tab[mod_name[i]]
  }

  unimportant_pct <- (mod_tab - feature_tab) / mod_length
  important_pct   <- feature_tab / mod_length

  importance_pct <- data.frame(
    Module     = rep(mod_name, 2),
    Status     = rep(c("Not selected", "Selected"), each = length(mod_tab)),
    Percentage = c(unimportant_pct, important_pct)
  )

  num_mods <- suppressWarnings(as.numeric(object$module_membership[, 2]))
  if (sum(is.na(num_mods)) == 0) {
    importance_pct[, 1] <- factor(importance_pct[, 1],
                                   levels = as.character(sort(unique(num_mods))))
  }

  Module <- Percentage <- Status <- NULL  # R CMD check

  p <- ggplot(importance_pct,
              aes(x = Module, y = Percentage, fill = Status)) +
    geom_col() +
    scale_fill_manual(values = c("Not selected" = "#cccccc", "Selected" = "#3e568a")) +
    ggtitle(main) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(plot.title   = element_text(lineheight = 0.8, face = "bold"),
          legend.title = element_blank())

  plot(p)
  invisible(p)
}


#' Plot Elbow Plot of Absolute SHAP Values for a Mossy Forest Object
#'
#' Ranks stable features by mean |SHAP| importance and displays an elbow
#' (scree) plot, optionally coloured by module membership.
#'
#' @param object        A \code{mossy_forest} object.
#' @param main          Plot title. Default \code{"Elbow Plot of Absolute SHAP Values"}.
#' @param xlab          X-axis label. Default \code{"Features"}.
#' @param ylab          Y-axis label. Default \code{"Absolute SHAP Values"}.
#' @param colorName     Legend title for module colour. Default \code{"Module"}.
#' @param color         Colour points by module membership? Default \code{TRUE}.
#' @param module_labels Optional two-column data frame of label mappings. Default \code{NULL}.
#' @param ...           Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot_elbow <- function(object, ...) {
  UseMethod("plot_elbow")
}

#' @describeIn plot_elbow Elbow plot for a \code{mossy_forest} object.
#' @export
plot_elbow.mossy_forest <- function(object,
                                      main          = NULL,
                                      xlab          = NULL,
                                      ylab          = NULL,
                                      colorName     = NULL,
                                      color         = TRUE,
                                      module_labels = NULL,
                                      ...) {

  if (is.null(main))      main      <- "Elbow Plot of Absolute SHAP Values"
  if (is.null(xlab))      xlab      <- "Features"
  if (is.null(ylab))      ylab      <- "Absolute SHAP Values"
  if (is.null(colorName)) colorName <- "Module"

  if (!is.null(module_labels)) {
    old_labels <- object$module_membership$module
    module_labels <- module_labels[order(module_labels[, 1]), ]
    new_labels <- as.character(factor(old_labels, labels = module_labels[, 2]))
    object$module_membership$module <- new_labels

    select_mods <- as.factor(object$feature_list$module_membership)
    select_module_table <- module_labels[
      module_labels[, 1] %in% levels(select_mods), , drop = FALSE]

    if ("." %in% levels(select_mods)) {
      dot_idx <- which(levels(select_mods) == ".")
      levels(select_mods)[-dot_idx] <- select_module_table[, 2]
    } else {
      levels(select_mods) <- select_module_table[, 2]
    }
    object$final_SHAP$module_membership <- as.character(select_mods)
  }

  plot_df <- object$final_SHAP %>%
    dplyr::mutate(rank = dplyr::row_number())

  if (color) {
    is_color <- function(x) x %in% colors()
    if (all(is_color(plot_df$module_membership))) {
      p <- ggplot(plot_df,
                  aes(x = rank, y = variable_importance,
                      colour = module_membership)) +
        geom_point(size = 3) +
        geom_line(aes(group = 1), colour = "gray60") +
        scale_colour_identity(guide = "legend")
    } else {
      p <- ggplot(plot_df,
                  aes(x = rank, y = variable_importance,
                      colour = module_membership)) +
        geom_point(size = 3) +
        geom_line(aes(group = 1), colour = "gray60") +
        scale_colour_discrete()
    }
  } else {
    p <- ggplot(plot_df,
                aes(x = rank, y = variable_importance)) +
      geom_point(size = 3) +
      geom_line(aes(group = 1), colour = "gray60")
  }

  p <- p +
    scale_x_continuous(breaks = plot_df$rank, labels = plot_df$feature_name) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(size = 14, hjust = 0.5)) +
    labs(x = xlab, y = ylab, colour = colorName, title = main)

  plot(p)
  invisible(p)
}


#' Plot Directional (Signed) Feature Importance for a Mossy Forest Object
#'
#' Diverging bar chart of signed importance = direction \eqn{\times} mean|SHAP|.
#' Magnitude is mean|SHAP|; the sign is the global effect direction, taken from
#' the Spearman correlation between each feature's value and its SHAP value
#' (\code{mean(signed SHAP)} \eqn{\approx 0} because contributions cancel, so it
#' is unusable for direction). Bars to the right increase the outcome, to the
#' left decrease it. Features whose direction is ambiguous (non-monotone:
#' \eqn{|\rho| < 0.1}) are greyed — a single direction is misleading for those.
#'
#' @param object    A \code{mossy_forest} object.
#' @param top_n     Show only the \code{top_n} features with the largest
#'   \eqn{|}signed importance\eqn{|}. Default \code{20L}; use \code{NULL} or
#'   \code{Inf} to plot every stable feature.
#' @param labels    Optional named character vector mapping feature names (as in
#'   \code{final_SHAP$feature_name}) to human-readable labels for the axis, e.g.
#'   a gene-symbol or indicator-name lookup. Names not found are left unchanged.
#'   Default \code{NULL} (show the raw feature names).
#' @param pos_color Bar colour for positive direction. Default \code{"#4393c3"}.
#' @param neg_color Bar colour for negative direction. Default \code{"#d6604d"}.
#' @param amb_color Bar colour for ambiguous features. Default \code{"#bbbbbb"}.
#' @param pos_label Legend text for the positive-direction group. Default
#'   \code{"increase"}.
#' @param neg_label Legend text for the negative-direction group. Default
#'   \code{"decrease"}.
#' @param amb_label Legend text for the ambiguous group. Default
#'   \code{"ambiguous"}.
#' @param main      Plot title. Default \code{"Directional feature importance"}.
#' @param ...       Ignored.
#'
#' @return A \code{ggplot2} object (invisibly).
#' @import ggplot2
#' @export
plot_signed_importance <- function(object, ...) {
  UseMethod("plot_signed_importance")
}

#' @describeIn plot_signed_importance Signed-importance plot for a
#'   \code{mossy_forest} object.
#' @export
plot_signed_importance.mossy_forest <- function(object,
                                                  top_n     = 20L,
                                                  labels    = NULL,
                                                  pos_color = "#4393c3",
                                                  neg_color = "#d6604d",
                                                  amb_color = "#bbbbbb",
                                                  pos_label = "increase",
                                                  neg_label = "decrease",
                                                  amb_label = "ambiguous",
                                                  main = "Directional feature importance",
                                                  ...) {
  df <- object$final_SHAP
  if (is.null(df) || !("signed_importance" %in% names(df)) || nrow(df) == 0L)
    stop("No directional SHAP data found. Refit with backend = \"python\", ",
         "or backend = \"R\" with r_shap = \"treeshap\" or \"fastshap\" ",
         "(r_shap = \"permutation\" has no per-observation SHAP values).")

  # Keep the top_n features by |signed importance| (largest effect either way)
  if (!is.null(top_n) && is.finite(top_n) && nrow(df) > top_n)
    df <- df[order(-abs(df$signed_importance)), , drop = FALSE][seq_len(top_n), , drop = FALSE]

  df <- df[order(abs(df$signed_importance)), , drop = FALSE]

  # Annotate: swap raw feature names for human-readable labels where supplied
  if (!is.null(labels)) {
    mapped <- unname(labels[as.character(df$feature_name)])
    df$feature_name <- ifelse(is.na(mapped), as.character(df$feature_name), mapped)
  }

  df$feature_name <- factor(df$feature_name, levels = df$feature_name)
  df$fill_grp <- ifelse(df$direction_ambiguous, "ambiguous",
                 ifelse(df$signed_importance >= 0, "positive", "negative"))

  signed_importance <- feature_name <- fill_grp <- NULL  # R CMD check

  p <- ggplot(df, aes(x = signed_importance, y = feature_name, fill = fill_grp)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(
      values = c(positive = pos_color, negative = neg_color, ambiguous = amb_color),
      breaks = c("positive", "negative", "ambiguous"),
      labels = c(paste0(pos_label, " (+)"), paste0(neg_label, " (−)"), amb_label),
      name   = NULL
    ) +
    labs(x = "Signed importance  (direction × mean |SHAP|)",
         y = NULL, title = main) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom")

  print(p)
  invisible(p)
}
