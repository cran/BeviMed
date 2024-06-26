% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{tune_temperatures}
\alias{tune_temperatures}
\title{Tune temperatures}
\usage{
tune_temperatures(number_of_temperatures, return_temperatures = FALSE, ...)
}
\arguments{
\item{number_of_temperatures}{Integer value giving number of tuned temperatures (including 0 and 1) to obtain.}

\item{return_temperatures}{Logical value determining whether to return just the numeric vector of tuned temperatures or to return the \code{BeviMed_m}-classed object containing the output of the MCMC sampling.}

\item{...}{Other arguments to pass to \code{call_cpp}.}
}
\value{
If \code{return_temperatures == TRUE}, a numeric vector of tuned temperatures, otherwise an object of class \code{BeviMed_m}.
}
\description{
Tune temperatures using interval bisection to minimimise Kullback-Liebler divergence between adjacent power posteriors
}
