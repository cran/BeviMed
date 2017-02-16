#' Calculate marginal probability of observed case-control status y under model gamma = 0
#' 
#' Marginal probability calculated exactly by integration.
#' @template y
#' @param tau0_shape Beta shape hyper-priors for prior on rate of case labels
#' @export
gamma0_evidence <- function(
	y,
	tau0_shape=c(1, 1)
) {
	lbeta(sum(y)+tau0_shape[1], length(y)-sum(y)+tau0_shape[2]) -
	lbeta(tau0_shape[1], tau0_shape[2])
}

#' Extract evidence for gamma = 1 model from \code{BeviMed} object.
#' 
#' @template x_BeviMed
#' @export
extract_gamma1_evidence <- function(x) {
	stopifnot(class(x) == "BeviMed")
	sum_ML_over_PP(x[["y_log_lik_t_equals_1"]], x[["temperatures"]])
}

#' Extract expected number of cases explained by pathogenic conigurations of alleles from \code{BeviMed} object.
#' 
#' @template x_BeviMed
#' @export
extract_expected_explained <- function(x) {
	stopifnot(class(x) == "BeviMed")
	if (("x" %in% names(x))*prod(dim(x[["x"]])) == 0)
		stop("Must make sure to set 'return_x_trace=TRUE' in call to 'bevimed' to use this function")
	mean(apply(x$x[,x$y,drop=FALSE], 1, sum))
}

#' Get expected number of cases explained by pathogenic configurations of alleles
#'
#' @param ... Arguments to be passed to \code{\link{bevimed}}.
#' @export
expected_explained <- function(...) {
	extract_expected_explained(bevimed(
		return_z_trace=FALSE,
		return_x_trace=TRUE,
		...
	))
}

#' Extract expected number of variants involved in cases explained by pathogenic conigurations of alleles from \code{BeviMed} object.
#' 
#' @template x_BeviMed
#' @export
extract_explaining_variants <- function(x) {
	stopifnot(class(x) == "BeviMed")
	if (("z" %in% names(x))*prod(dim(x[["z"]])) == 0)
		stop("Must make sure to set 'return_z_trace=TRUE' in call to 'bevimed' to use this function")
	mean(apply(x$z[,x$k * (length(x$temperatures)-1) + unique(x$variant_table$variant[x$variant_table$case %in% which(x$y)]),drop=FALSE], 1, sum))
}

#' Get expected number of variants involved in cases explained by pathogenic conigurations of alleles
#'
#' @param ... Arguments to be passed to \code{\link{bevimed}}.
#' @export
explaining_variants <- function(...) {
	extract_explaining_variants(bevimed(return_z_trace=TRUE, return_x_trace=FALSE, ...))
}

#' Calculate marginal probability of observed case-control status under model gamma = 1 
#' 
#' @template y
#' @template G_matrix
#' @template min_ac
#' @param ... Other arguments to pass to \code{\link{bevimed}}.
#' @export
gamma1_evidence <- function(
	y,
	G,
	min_ac=1L,
	...
) {
	bv <- bevimed(
		return_z_trace=FALSE,
		return_x_trace=FALSE,
		y=y,
		G=G,
		min_ac=min_ac,
		...
	)
	extract_gamma1_evidence(bv)
}

#' Calculate log Bayes factor between models gamma=1 and gamma=0 given the data
#' 
#' @template y
#' @param ... Other arguments to pass to \code{\link{bevimed}}.
#' @template tau0_shape
#' @return Log Bayes factor.
#' @export
log_BF <- function(
	y,
	tau0_shape=c(1, 1),
	...
) {
	gamma1_evidence(
		y=y,
		...
	) - 
	gamma0_evidence(y, tau0_shape=tau0_shape)
}

#' Calculate probability of an association between presence/absence of local genotype configuration and case-control label
#'
#' @template y
#' @template G_matrix
#' @template ploidy 
#' @template prior_prob_association
#' @template prior_prob_dominant
#' @template tau0_shape
#' @template by_MOI
#' @param ... Other arguments to pass to \code{\link{log_BF}}. 
#' @return Probability of association.
#' @export
prob_association <- function(
	y,
	G,
	ploidy=rep(2L, length(y)),
	prior_prob_association=0.01,
	prior_prob_dominant=0.5,
	tau0_shape=c(1, 1),
	by_MOI=FALSE,
	...
) {
	stopifnot(is.matrix(G) & is.numeric(G))
	stopifnot(nrow(G) == length(y))
	stopifnot(prior_prob_association <= 1 & prior_prob_association >= 0)
	stopifnot(prior_prob_dominant <= 1 & prior_prob_dominant >= 0)

	dominant_prior <- prior_prob_association * prior_prob_dominant
	recessive_prior <- prior_prob_association * (1-prior_prob_dominant)

	dom_ev <- if (dominant_prior > 0) gamma1_evidence(y=y, G=G, min_ac=1L, ...) else -Inf
	rec_ev <- if (recessive_prior > 0) gamma1_evidence(y=y, G=G, min_ac=ploidy, ...) else -Inf
	g0_ev <- gamma0_evidence(y=y, tau0_shape=tau0_shape)

	modal_model <- max(c(dom_ev, rec_ev, g0_ev))
	num_dom <- dominant_prior * exp(dom_ev-modal_model)
	num_rec <- recessive_prior * exp(rec_ev-modal_model)
	num_g0 <- (1-dominant_prior-recessive_prior) * exp(g0_ev-modal_model)
	denom <- num_dom + num_rec + num_g0

	#posterior
	moi_probs <- c(Dominant=num_dom/denom, Recessive=num_rec/denom)
	if (by_MOI) moi_probs
	else sum(moi_probs)
}

#' Calculate probability of an association
#'
#'Synonym of \code{\link{prob_association}}
#'
#' @param ... Arguments passed to \code{\link{prob_association}}
#' @return Probability of association.
#' @export
gamma1_prob <- function(...) prob_association(...)

#' Calculate probability of pathogencity for variants in region given a prior probability of association between case label and the region
#'
#' @template y
#' @template G_matrix
#' @template ploidy 
#' @template prior_prob_association
#' @template prior_prob_dominant
#' @template tau0_shape
#' @template by_MOI
#' @param ... Other arguments to pass to \code{\link{bevimed}}. 
#' @return Probabilities of pathogenicity.
#' @export
prob_pathogenic <- function(
	y,
	G,
	ploidy=rep(2L, length(y)),
	prior_prob_association=0.01,
	prior_prob_dominant=0.5,
	tau0_shape=c(1, 1),
	by_MOI=FALSE,
	...
) {
	stopifnot(is.matrix(G) & is.numeric(G))
	stopifnot(nrow(G) == length(y))
	stopifnot(prior_prob_association <= 1 & prior_prob_association >= 0)
	stopifnot(prior_prob_dominant <= 1 & prior_prob_dominant >= 0)

	dominant_prior <- prior_prob_association * prior_prob_dominant
	recessive_prior <- prior_prob_association * (1-prior_prob_dominant)

	dom_bv <- bevimed(
		return_z_trace=TRUE,
		return_x_trace=FALSE,
		y=y,
		G=G,
		min_ac=1L,
		...
	)

	rec_bv <- bevimed(
		return_z_trace=TRUE,
		return_x_trace=TRUE,
		y=y,
		G=G,
		min_ac=ploidy,
		...
	)

	dom_ev <- sum_ML_over_PP(dom_bv[["y_log_lik_t_equals_1"]], dom_bv[["temperatures"]])
	rec_ev <- sum_ML_over_PP(rec_bv[["y_log_lik_t_equals_1"]], rec_bv[["temperatures"]])
	g0_ev <- gamma0_evidence(y, tau0_shape=tau0_shape)

	modal_model <- max(c(dom_ev, rec_ev, g0_ev))
	num_dom <- dominant_prior * exp(dom_ev-modal_model)
	num_rec <- recessive_prior * exp(rec_ev-modal_model)
	num_g0 <- (1-dominant_prior-recessive_prior) * exp(g0_ev-modal_model)
	denom <- num_dom + num_rec + num_g0

	dom_conditional_prob_pathogenic <- apply(dom_bv[["z"]][,dom_bv[["k"]]*(length(dom_bv[["temperatures"]])-1)+1:dom_bv[["k"]],drop=FALSE], 2, mean)
	rec_conditional_prob_pathogenic <- apply(rec_bv[["z"]][,rec_bv[["k"]]*(length(rec_bv[["temperatures"]])-1)+1:rec_bv[["k"]],drop=FALSE], 2, mean)

	prob_by_moi <- cbind(
		Dominant=dom_conditional_prob_pathogenic*num_dom/denom,
		Recessive=rec_conditional_prob_pathogenic*num_rec/denom
	)
	if (by_MOI) prob_by_moi
	else apply(prob_by_moi, 1, sum)
}

#' Calculate probability of pathogencity for variants in region given an association between case label and the region
#'
#' @template y
#' @template G_matrix
#' @template min_ac
#' @param ... Other arguments to pass to \code{\link{bevimed}}. 
#' @return Probabilities of pathogenicity.
#' @export
conditional_prob_pathogenic <- function(
	y,
	G,
	min_ac=1L,
	...
) {
	bv <- bevimed(
		return_z_trace=TRUE,
		return_x_trace=FALSE,
		y=y,
		G=G,
		min_ac=min_ac,
		...
	)

	apply(bv[["z"]][,bv[["k"]]*(length(bv[["temperatures"]])-1)+1:bv[["k"]],drop=FALSE], 2, mean)
}

