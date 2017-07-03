\name{BeviMed-package}
\alias{BeviMed-package}
\alias{BeviMed}
\docType{package}
\title{
Bayesian Evaluation of Variant Involvement in Mendelian Disease 
}
\description{
A fast integrative genetic association test for rare diseases.
}
\details{
BeviMed estimates a probability of association between a case/control label and allele counts at rare variant sites in a genomic locus and also, given that there is an association, the probabilities that each variant is involved in the disease. It does so by estimating the evidence for a model where the case/control label is independent of the allele configurations, and a model in which the probability of the case/control label depends on the corresponding allele configuration and a latent partition of variants into pathogenic and non-pathogenic groups.
}
\author{
Daniel Greene.

Maintainer: Daniel Greene <dg333@cam.ac.uk>
}
\references{
Greene et al., A Fast Association Test for Identifying Pathogenic Variants Involved in Rare Diseases, The American Journal of Human Genetics (2017), http://dx.doi.org/10.1016/j.ajhg.2017.05.015.
}
\keyword{ package, statistical association test, Bayesian inference, rare disease  }
\seealso{
  \code{\link{bevimed}}
}
\examples{
# generate random allele count matrix for 10 individuals at 5 rare variant sites
G <- matrix(rbinom(n=50, size=2, prob=0.1), nrow=10, ncol=5)
# generate case/control label whereby individuals with any alleles in the first 
# three variants are labelled cases
y <- apply(G[,1:3], 1, sum) > 0
# apply BeviMed
bevimed(y=y, G=G)
}
