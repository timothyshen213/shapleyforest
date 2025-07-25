% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shap_rfe.R
\name{screen_RFE}
\alias{screen_RFE}
\title{Screening step of shapley forest through SHAP values.}
\usage{
screen_RFE(
  X,
  y,
  module_list,
  module_membership,
  screen_control,
  select_control,
  shap_model = "full",
  CLASSIFICATION,
  min_features,
  nsim,
  nodesize,
  num_processors,
  verbose,
  debug,
  parallel,
  seed
)
}
\arguments{
\item{X}{A data.frame. With columns denoting a feature vector.
Could include additional covariates not a part of
the original modules.}

\item{y}{A Response vector.}

\item{module_list}{A list containing the modules.}

\item{module_membership}{A vector that specifies the module membership for each
each feature. See \code{wsf} for possible method.}

\item{screen_control}{Defines the parameter settings for the screening step.}

\item{select_control}{Defines the parameter setting for the selection step.}

\item{shap_model}{Binary indicator for \code{sf} model. If \code{full}, \code{sf}
runs SHAPley values at both screening and selection step.
If \code{after}, \code{sf} model runs SHAPley values at the end
of final model and keeps permutation VIMs usage at other steps.
\code{full} is default.}

\item{CLASSIFICATION}{Binary indicator if the model is a classification model.}

\item{min_features}{Defines minimum feature allowed for each module.}

\item{nsim}{Number of Monte Carlo repetitions for estimating SHAP
values in the screening step. Default is \code{1}. Increasing
\code{nsim} leads to more accurate results, but at the cost
of computational cost.}

\item{nodesize}{Minimum \code{nodesize}}

\item{num_processors}{Number of processors used to fit random forests.}

\item{verbose}{Defines the warning message protocol. See \code{sf} for information.}

\item{debug}{Sets the debugging procedures. See \code{sf} for more information}

\item{parallel}{Binary indicator if parallel is called.}

\item{seed}{Seed to be set for reproducibility. Note: for a given seed,
parallel and non-parallel versions will differ.}
}
\value{
A data.frame with final surviving features.
}
\description{
Runs the screening step of shapley forest algorithm.
Returns data.frame with the final surviving SHAP values and features.
}
