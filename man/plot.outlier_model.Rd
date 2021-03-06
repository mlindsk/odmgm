% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_outlier.R
\name{plot.outlier_model}
\alias{plot.outlier_model}
\title{Plot the probability mass function of the deviances}
\usage{
\method{plot}{outlier_model}(x, alpha = 0.05, sig_col = "#FF0000A0", ...)
}
\arguments{
\item{x}{A \code{outlier_model} object}

\item{alpha}{Significance level}

\item{sig_col}{Color of the significance level area (default is red)}

\item{...}{Not used (for S3 compatibility)}
}
\description{
A plot method to show the pmf of the approximated pmf of \code{D(Y)}
}
