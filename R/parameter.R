#' WGCNA Parameter Organization
#' 
#' TO DO
#' 
WGCNA_control <- function(power=6, min_features=20, ...) {
  extra_args <- list(...)
  obj <- list()
  obj$power <- power
  obj$min_features <- min_features
  obj$extra_args <- extra_args
  class(obj) <- "WGCNA_control"
  return(obj)
}
