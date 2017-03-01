#' @title Bayesian Evaluation of Variant Involvement in Mendelian Disease 
#'
#' @description Infer probabilities of association between disease label and locus and posterior parameter values under BeviMed model.
#' 
#' @template y
#' @template G_matrix
#' @template ploidy
#' @template prior_prob_association
#' @template prior_prob_dominant
#' @template tau0_shape
#' @param dominant_args Arguments to pass to \code{\link{bevimed_m}} conditioning on dominant inheritance.
#' @param recessive_args Arguments to pass to \code{\link{bevimed_m}} conditioning on recessive inheritance.
#' @param ... Arguments to be passed to \code{\link{bevimed_m}} for both modes of inheritance.
#' @return \code{BeviMed} object containing results of inference.
#' @export
#' @seealso \code{\link{prob_association}}, \code{\link{bevimed_m}}, \code{\link{summary.BeviMed}}
bevimed <- function(
	y,
	G,
	ploidy=rep(2L, length(y)),
	prior_prob_association=0.01,
	prior_prob_dominant=0.5,
	tau0_shape=c(1, 1),
	dominant_args=NULL,	
	recessive_args=NULL,
	...
) {
	if (prior_prob_association > 1 | prior_prob_association < 0)
		stop("'prior_prob_association' must be between 0 and 1")
	if (prior_prob_dominant > 1 | prior_prob_dominant < 0)
		stop("'prior_prob_dominant' must be between 0 and 1")

	dominant_prior <- prior_prob_association * prior_prob_dominant
	recessive_prior <- prior_prob_association * (1-prior_prob_dominant)
	shared <- list(...)

	G_args <- get_G_args(G)

	structure(
		class="BeviMed",
		list(
			parameters=list(
				tau0_shape=tau0_shape,
				prior_prob_association=prior_prob_association,
				prior_prob_dominant=prior_prob_dominant,
				ploidy=ploidy,
				y=y,
				variant_table=to_var_tab(G_args)
			),
			moi=Map(
				f=function(min_ac, args) do.call(what=bevimed_m, c(list(y=y, G=G, min_ac=min_ac), args, shared)),
				list(dominant=1L, recessive=ploidy),
				list(dominant_args, recessive_args)
			)
		)
	)
}

#' Calculate marginal probability of observed case-control status y under model gamma = 0
#' 
#' Marginal probability calculated exactly by integration.
#'
#' @template y
#' @param tau0_shape Beta shape hyper-priors for prior on rate of case labels
#' @return Log marginal likelihood.
#' @export
#' @seealso \code{\link{bevimed}}, \code{\link{gamma1_evidence}}
gamma0_evidence <- function(
	y,
	tau0_shape=c(1, 1)
) {
	lbeta(sum(y)+tau0_shape[1], length(y)-sum(y)+tau0_shape[2]) -
	lbeta(tau0_shape[1], tau0_shape[2])
}

#' @title Extract evidence for model gamma = 1
#' 
#' @description Extract evidence from \code{BeviMed_m} object.
#' @template x_BeviMed_m
#' @return Log marginal likelihood.
#' @export
#' @seealso \code{\link{gamma1_evidence}}, \code{\link{bevimed_m}}
extract_gamma1_evidence <- function(x) {
	stopifnot(class(x) == "BeviMed_m")
	sum_ML_over_PP(x[["traces"]][["y_log_lik_t_equals_1"]], x[["parameters"]][["temperatures"]])
}

#' @title Calculate evidence under model gamma = 1 
#'
#' @description Use \code{\link{bevimed_m}} to perform inference under model gamma = 1 and return only the log evidence/integrated likelihood.
#' 
#' @template dots_to_bevimed_m
#' @return Log marginal likelihood.
#' @export
#' @seealso \code{\link{bevimed_m}}, \code{\link{extract_gamma1_evidence}}
gamma1_evidence <- function(
	...
) {
	bv <- bevimed_m(
		return_z_trace=FALSE,
		return_x_trace=FALSE,
		...
	)
	extract_gamma1_evidence(bv)
}

#' @title Extract expected number of explained cases
#'
#' @description Extract expected number of cases explained by pathogenic configurations of alleles from \code{BeviMed_m} object.
#' 
#' @template x_BeviMed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{expected_explained}}, \code{\link{bevimed_m}}
extract_expected_explained <- function(x) {
	stopifnot(class(x) == "BeviMed_m")
	if (("x" %in% names(x[["traces"]]))*prod(dim(x[["traces"]][["x"]])) == 0)
		stop("Must make sure to set 'return_x_trace=TRUE' in call to 'bevimed_m' to use this function")
	mean(apply(x$traces$x[,x$parameters$y,drop=FALSE], 1, sum))
}

#' @title Calculate expected number of explained cases
#'
#' @description Use \code{\link{bevimed_m}} to perform inference under model gamma = 1 and return only the expected number of cases explained by pathogenic allele configurations.
#'
#' @template dots_to_bevimed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{bevimed_m}}, \code{\link{extract_expected_explained}}
expected_explained <- function(...) {
	extract_expected_explained(bevimed_m(
		return_z_trace=FALSE,
		return_x_trace=TRUE,
		...
	))
}

#' @title Extract expected number of pathogenic variants in cases
#' @description Extract expected number of variants involved in cases explained by pathogenic conigurations of alleles from \code{BeviMed_m} object.
#' 
#' @template x_BeviMed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{explaining_variants}}, \code{\link{bevimed_m}}
extract_explaining_variants <- function(x) {
	stopifnot(class(x) == "BeviMed_m")
	if (("z" %in% names(x[["traces"]]))*prod(dim(x[["traces"]][["z"]])) == 0)
		stop("Must make sure to set 'return_z_trace=TRUE' in call to 'bevimed_m' to use this function")
	mean(apply(x$traces$z[,x$parameters$k * (length(x$parameters$temperatures)-1) + unique(x$parameters$variant_table$variant[x$parameters$variant_table$case %in% which(x$parameters$y)]),drop=FALSE], 1, sum))
}

#' Calculate expected number of pathogenic variants in cases
#'
#' @description Use \code{\link{bevimed_m}} to perform inference under model gamma = 1 and return only the expected number of pathogenic variants in cases.
#'
#' @template dots_to_bevimed_m
#' @return Numeric value.
#' @export
#' @seealso \code{\link{extract_explaining_variants}}, \code{\link{bevimed_m}}
explaining_variants <- function(...) {
	extract_explaining_variants(bevimed_m(return_z_trace=TRUE, return_x_trace=FALSE, ...))
}

#' Calculate log Bayes factor between model gamma = 1 conditional on a given mode of inheritance and model gamma = 0
#' 
#' @description Compute log Bayes factor of model gamma = 1 against model gamma = 0 for just one mode of inheritance. 
#'
#' @template y
#' @template dots_to_bevimed_m
#' @template tau0_shape
#' @return Log Bayes factor.
#' @export
#' @seealso \code{\link{bevimed_m}}, \code{\link{prob_association_m}}
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

#' @title Calculate probability of association for one mode of inheritance
#'
#' @description Equivalent to \code{\link{prob_association}} where the prior probability of one mode of inheritance is 1. This function is faster, as it only calls \code{\link{bevimed_m}} once.
#'
#' @template y
#' @template min_ac
#' @template prior_prob_association
#' @param ... Other arguments to pass to \code{\link{log_BF}}.
#' @return Probability value.
#' @export
#' @seealso \code{\link{log_BF}}, \code{\link{prob_association}}, \code{\link{bevimed_m}}
prob_association_m <- function(
	y,
	min_ac=1L,
	prior_prob_association=0.01,
	...
) {
	bf <- log_BF(y, min_ac=min_ac, ...)
	num <- prior_prob_association*exp(bf)
	num/(1-prior_prob_association+num)
}

#' @title Extract the posterior probability of association
#'
#' @description Get posterior probability of association as numeric value, or optionally as numeric vector of length two with probabilities broken down by mode of inheritance (by passing \code{by_MOI=TRUE}), from a \code{BeviMed} object.
#' @template x_BeviMed
#' @template by_MOI
#' @return Probability values.
#' @export
#' @seealso \code{\link{prob_association}}, \code{\link{bevimed}}
extract_prob_association <- function(x, by_MOI=FALSE) {
	stopifnot(class(x) == "BeviMed")
	dominant_prior <- x[["parameters"]][["prior_prob_association"]] * x[["parameters"]][["prior_prob_dominant"]]
	recessive_prior <- x[["parameters"]][["prior_prob_association"]] - dominant_prior 

	dom_ev <- extract_gamma1_evidence(x[["moi"]][["dominant"]])
	rec_ev <- extract_gamma1_evidence(x[["moi"]][["recessive"]])
	g0_ev <- gamma0_evidence(y=x[["parameters"]][["y"]], tau0_shape=x[["parameters"]][["tau0_shape"]])

	modal_model <- max(c(dom_ev, rec_ev, g0_ev))
	num_dom <- dominant_prior * exp(dom_ev-modal_model)
	num_rec <- recessive_prior * exp(rec_ev-modal_model)
	num_g0 <- (1-dominant_prior-recessive_prior) * exp(g0_ev-modal_model)
	denom <- num_dom + num_rec + num_g0

	moi_probs <- c(dominant=num_dom/denom, recessive=num_rec/denom)
	if (by_MOI) moi_probs
	else sum(moi_probs)
}

#' @title Calculate probability of association
#'
#' @description Calculate probability of an association between case/control label and allele configuration, optionally broken down by mode of inheritance.
#' @template by_MOI
#' @param ... Arguments to pass to \code{\link{bevimed}}. 
#' @return Probability of association.
#' @export
#' @seealso \code{\link{bevimed}}, \code{\link{extract_prob_association}}
prob_association <- function(
	by_MOI=FALSE,
	...
) {
	bv <- bevimed(..., return_z_trace=FALSE, return_x_trace=FALSE)
	extract_prob_association(bv, by_MOI=by_MOI)
}

#' @title Extract variant marginal probabilities of pathogenicity
#'
#' @description Extract the marginal probability of pathogenicity for individual variants from \code{BeviMed} object, optionally broken down by mode of inheritance.
#'
#' @template x_BeviMed
#' @template by_MOI
#' @return A vector of probabilities of pathogenicity for individual variants, or if \code{by_MOI} is \code{TRUE}, then a matrix of probabilities, with rows corresponding to modes of inheritance and columns to variants.
#' @export
#' @seealso \code{\link{prob_pathogenic}}, \code{\link{bevimed}}
extract_prob_pathogenic <- function(x, by_MOI=TRUE) {
	probs <- extract_prob_association(x, by_MOI=by_MOI)
	apply(mapply(FUN="*", probs, lapply(x$moi, extract_conditional_prob_pathogenic)), 1, if (by_MOI) identity else sum)
}

#' Calculate variant marginal probabilities of pathogencity 
#'
#' Calls \code{\link{bevimed}} and \code{\link{extract_prob_pathogenic}} to obtain marginal probabilities of pathogenicity.
#'
#' @template by_MOI
#' @param ... Arguments to pass to \code{\link{bevimed}}. 
#' @return If \code{by_MOI} is \code{FALSE}, a vector of probabilities of pathogenicity for each variant, otherwise a matrix of probabilities of pathogenicity conditional on mode of inheritance, with a column for each variant and two rows corresponding to each mode of inheritance.
#' @export
#' @seealso \code{\link{extract_prob_pathogenic}}, \code{\link{bevimed}}
prob_pathogenic <- function(
	by_MOI=FALSE,
	...
) {
	bv <- bevimed(..., return_z_trace=TRUE, return_x_trace=FALSE)
	extract_prob_pathogenic(bv, by_MOI=by_MOI)
}

#' @title Extract probability of pathogenicity for variant conditional on a given mode of inheritance
#'
#' @description Extract the probability of pathogenicity for individual variants from a \code{BeviMed_m} object.
#'
#' @template x_BeviMed_m
#' @return Vector of probabilities of pathogenicity for individual variants.
#' @export
#' @seealso \code{\link{conditional_prob_pathogenic}}, \code{\link{bevimed_m}}
extract_conditional_prob_pathogenic <- function(x) {
	k <- x[["parameters"]][["k"]]
	apply(x[["traces"]][["z"]][,k*(length(x[["parameters"]][["temperatures"]])-1L)+seq(length.out=k),drop=FALSE], 2, mean)
}

#' Calculate probability of pathogencity for variants conditional on mode of inheritance.
#'
#' Calls \code{\link{bevimed_m}} and \code{\link{extract_conditional_prob_pathogenic}} to obtain probabilities of pathogenicity.
#' @template dots_to_bevimed_m
#' @return Probabilities of pathogenicity.
#' @export
#' @seealso \code{\link{extract_conditional_prob_pathogenic}}, \code{\link{bevimed_m}}
conditional_prob_pathogenic <- function(
	...
) {
	bv <- bevimed_m(
		return_z_trace=TRUE,
		return_x_trace=FALSE,
		...
	)
	extract_conditional_prob_pathogenic(bv)
}

