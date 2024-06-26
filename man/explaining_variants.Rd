% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience.R
\name{explaining_variants}
\alias{explaining_variants}
\title{Calculate expected number of pathogenic variants in cases}
\usage{
explaining_variants(...)
}
\arguments{
\item{...}{Arguments to pass to \code{\link{bevimed_m}}.}
}
\value{
Numeric value.
}
\description{
Use \code{\link{bevimed_m}} to perform inference under model gamma = 1 and return only the expected number of pathogenic variants in cases.
}
\seealso{
\code{\link{extract_explaining_variants}}, \code{\link{bevimed_m}}
}
