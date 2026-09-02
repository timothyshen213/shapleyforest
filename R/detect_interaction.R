#' Detect Pairwise Interactions from a Mossy Forest Object
#'
#' Returns exact TreeSHAP pairwise interaction strengths for all stable
#' features.  Interaction values are computed via
#' \code{TreeExplainer.shap_interaction_values()} during model fitting and
#' stored in the object — no re-computation is needed here.
#'
#' The \eqn{p \times p} \strong{interaction matrix} stored in
#' \code{object$shap_obj$interaction_matrix} is defined as:
#' \deqn{M_{jk} = \frac{1}{n}\sum_{i=1}^{n}|\phi_{jk}(x_i)|}
#' where \eqn{\phi_{jk}(x_i)} is the TreeSHAP interaction value for the pair
#' \eqn{(j, k)} at observation \eqn{i} (Lundberg et al., 2018).  Diagonal
#' entries are main effects; off-diagonal entries are pairwise interactions.
#' The matrix is symmetric: \eqn{M_{jk} = M_{kj}}.
#'
#' @param object  A \code{mossy_forest} object.
#' @param thresh  Interaction-strength threshold.  Only pairs whose mean
#'                \eqn{|\phi_{jk}|} exceeds \code{thresh} are returned as
#'                significant.  Set to \code{0} to return all pairs.
#' @param all     If \code{TRUE}, return every off-diagonal pair (upper
#'                triangle) regardless of \code{thresh}.  Default \code{FALSE}.
#' @param verbose If \code{TRUE}, print the significant-pair table (or a
#'                message if none are found). Default \code{TRUE}.
#' @param ...     Ignored.
#'
#' @return A data frame with columns \code{FeatureA}, \code{FeatureB}, and
#'   \code{Interaction_Strength} (mean absolute TreeSHAP interaction), sorted
#'   by interaction strength descending.  Invisibly when \code{all = TRUE}.
#'
#' @references
#' Lundberg, S. M., Erion, G., Chen, H., DeGrave, A., Prutkin, J. M.,
#'   Nair, B., Katz, R., Himmelfarb, J., Bansal, N., & Lee, S.-I. (2020).
#'   From local explanations to global understanding with explainable AI for
#'   trees. \emph{Nature Machine Intelligence}, 2, 56–67.
#'
#' @export
detect_interaction <- function(object, thresh = 0, all = FALSE,
                                verbose = TRUE, ...) {

  int_mat <- object$shap_obj$interaction_matrix
  if (is.null(int_mat) || !is.matrix(int_mat))
    stop(paste0(
      "No interaction matrix found in this mossy_forest object.\n",
      "  Refit with mf() — interaction values are computed automatically ",
      "during the final SHAP pass."
    ), call. = FALSE)

  feats <- rownames(int_mat)
  p     <- length(feats)
  if (p < 2L)
    stop("At least two stable features are required to report interactions.",
         call. = FALSE)

  # ── build upper-triangle long data frame ────────────────────────────────────
  rows <- vector("list", p * (p - 1L) / 2L)
  k    <- 1L
  for (j in seq_len(p - 1L)) {
    for (l in (j + 1L):p) {
      rows[[k]] <- data.frame(
        FeatureA             = feats[j],
        FeatureB             = feats[l],
        Interaction_Strength = int_mat[j, l],
        stringsAsFactors     = FALSE
      )
      k <- k + 1L
    }
  }
  all_pairs <- do.call(rbind, rows)
  all_pairs <- all_pairs[order(-all_pairs$Interaction_Strength), ]
  rownames(all_pairs) <- NULL

  if (all) return(invisible(all_pairs))

  # ── filter by threshold ──────────────────────────────────────────────────────
  sig_pairs <- all_pairs[all_pairs$Interaction_Strength > thresh, , drop = FALSE]
  rownames(sig_pairs) <- NULL

  if (verbose) {
    if (nrow(sig_pairs) > 0L)
      print(sig_pairs)
    else
      message(sprintf(
        "No interactions found above threshold (%.4g). ",
        thresh
      ), "Use a lower thresh or all=TRUE to see all pairs.")
  }

  invisible(sig_pairs)
}
