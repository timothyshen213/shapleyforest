# ── Internal Python environment ───────────────────────────────────────────────
.sf_py <- new.env(parent = emptyenv())
.sf_py$loaded <- FALSE

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Shapley Forest loaded. Run sf_setup() to configure the Python environment."
  )
}

#' Configure the Python environment for Shapley Forest
#'
#' Tells Shapley Forest which Python installation to use for its backend.
#' Call this once at the start of each R session before using \code{\link{sf}}.
#'
#' @param python     Path to a Python binary. Passed to
#'                   \code{reticulate::use_python()}.
#' @param condaenv   Name of a conda environment. Passed to
#'                   \code{reticulate::use_condaenv()}.
#' @param virtualenv Path to a virtual environment. Passed to
#'                   \code{reticulate::use_virtualenv()}.
#'
#' @details
#' The Python environment must have \code{scikit-learn}, \code{shap}, and
#' \code{numpy} installed. Example setup:
#' \preformatted{
#'   conda create -n sfenv python=3.11
#'   conda activate sfenv
#'   pip install scikit-learn shap numpy
#' }
#'
#' @return Invisible \code{NULL}.
#' @export
#' @examples
#' \dontrun{
#'   sf_setup(condaenv = "sfenv")
#'   sf_setup(python = "/usr/local/bin/python3")
#' }
sf_setup <- function(python = NULL, condaenv = NULL, virtualenv = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required. Install with: install.packages('reticulate')",
         call. = FALSE)

  if (!is.null(condaenv))
    reticulate::use_condaenv(condaenv, required = TRUE)
  else if (!is.null(virtualenv))
    reticulate::use_virtualenv(virtualenv, required = TRUE)
  else if (!is.null(python))
    reticulate::use_python(python, required = TRUE)

  .load_python_backend()
  invisible(NULL)
}

# Internal: load the bundled Python backend into .sf_py
.load_python_backend <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required. Install with: install.packages('reticulate')",
         call. = FALSE)

  # backend verbose output contains non-ASCII characters; without this,
  # ASCII-locale sessions crash with UnicodeEncodeError at verbose >= 2
  if (!nzchar(Sys.getenv("PYTHONIOENCODING")))
    Sys.setenv(PYTHONIOENCODING = "utf-8")

  py_file <- system.file("python", "sf_python_backend.py",
                          package = "shapleyforest")
  if (!nzchar(py_file) || !file.exists(py_file))
    stop("Shapley Forest Python backend not found in package installation.",
         call. = FALSE)

  reticulate::source_python(py_file, envir = .sf_py)
  .sf_py$loaded <- TRUE
  invisible(NULL)
}

# Internal: ensure Python backend is ready, loading lazily on first call
.ensure_python <- function() {
  if (isTRUE(.sf_py$loaded)) return(invisible(NULL))
  .load_python_backend()
}
