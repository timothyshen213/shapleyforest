% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shapff_object.R
\name{predict.shapley_forest}
\alias{predict.shapley_forest}
\title{Prediction for shapley forest}
\usage{
\method{predict}{shapley_forest}(object, new_data, ...)
}
\arguments{
\item{object}{A \code{shapley_forest} object.}

\item{new_data}{A matrix or data.frame containing new data.
Note: Feature names must match between
training set and test set data.frame.}

\item{...}{Obsolete additional arguments.}
}
\value{
A vector of predictions through random forest.
}
\description{
Prediction function for shapley forest algorithm.
Runs prediction through the final random forest
after both screening and selection steps.
}
