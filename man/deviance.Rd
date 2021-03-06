% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_deviance.R
\name{deviance}
\alias{deviance}
\alias{deviance.pure_disc}
\alias{deviance.pure_cont}
\alias{deviance.mixed}
\title{Calculate deviance}
\usage{
deviance(x, z, ...)

\method{deviance}{pure_disc}(x, z, ...)

\method{deviance}{pure_cont}(x, z, ...)

\method{deviance}{mixed}(x, z, ...)
}
\arguments{
\item{x}{A \code{outlier_model} object}

\item{z}{An observation (data.frame or named character vector)}

\item{...}{Not used (for S3 compatability)}
}
\description{
This function calculates the affine value \code{-2 log} likelihood-ratio statistic; also called the deviance statistic
}
