% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{MR_Robin_resample}
\alias{MR_Robin_resample}
\title{Obtain p-value from MR-Robin results}
\usage{
MR_Robin_resample(MR_Robin_res, nsamp = 10000)
}
\arguments{
\item{MR_Robin_res}{a list object returned by \code{MR_Robin}.}

\item{nsamp}{integer of the number of samples to use in estimating \eqn{P}-value
using a null distribution.}
}
\value{
A list of two elements:
\tabular{ll}{
\code{pvalue} \tab numeric of the estimated \eqn{P}-value.\cr
\code{nsamp_used} \tab integer of the number of samples used in estimating the \eqn{P}-value.\cr
}

\code{nsamp_used} is returned because not all \code{nsamp} may be used (samples will be dropped
if the model does not converge or results in a singular fit of the random slope).
}
\description{
Uses a resampling procedure to estimate a \eqn{P}-value for a MR-Robin object.
}
