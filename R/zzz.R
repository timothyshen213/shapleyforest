# ── Internal Python environment ───────────────────────────────────────────────
.mf_py <- new.env(parent = emptyenv())
.mf_py$loaded <- FALSE

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Mossy Forest loaded. Run mf_setup() to configure the Python environment."
  )
}

#' Configure the Python environment for Mossy Forest
#'
#' Tells Mossy Forest which Python installation to use for its backend.
#' Call this once at the start of each R session before using \code{\link{mf}}.
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
#'   mf_setup(condaenv = "sfenv")
#'   mf_setup(python = "/usr/local/bin/python3")
#' }
mf_setup <- function(python = NULL, condaenv = NULL, virtualenv = NULL) {
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

# Internal: load the bundled Python backend into .mf_py
.load_python_backend <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required. Install with: install.packages('reticulate')",
         call. = FALSE)

  # backend verbose output contains non-ASCII characters; without this,
  # ASCII-locale sessions crash with UnicodeEncodeError at verbose >= 2
  if (!nzchar(Sys.getenv("PYTHONIOENCODING")))
    Sys.setenv(PYTHONIOENCODING = "utf-8")

  py_file <- system.file("python", "mf_python_backend.py",
                          package = "mossyforest")
  if (!nzchar(py_file) || !file.exists(py_file))
    stop("Mossy Forest Python backend not found in package installation.",
         call. = FALSE)

  reticulate::source_python(py_file, envir = .mf_py)
  .mf_py$loaded <- TRUE
  invisible(NULL)
}

# Internal: ensure Python backend is ready, loading lazily on first call
.ensure_python <- function() {
  if (isTRUE(.mf_py$loaded)) return(invisible(NULL))
  .load_python_backend()
}
