#' Create summary of \code{BeviMed} classed-object
#'
#' @param object Object of class \code{BeviMed}.
#' @param gamma1_prior Prior probability of that gamma = 1.
#' @template tau0_shape
#' @template confidence
#' @template simulations
#' @param ... Non-used arguments.
#' @return Object of class \code{BeviMed_summary}.
#' @method summary BeviMed
#' @export
summary.BeviMed <- function(object, gamma1_prior=0.01, tau0_shape=object[["tau_shape"]], confidence=0.95, simulations=1000, ...) {
	vt <- object[["variant_table"]]
	y <- object[["y"]]
	variant_counts <- lapply(setNames(nm=c(F,T)), function(y_is) as.integer(table(factor(vt$variant[y[vt$case]==y_is], levels=seq(object[["k"]])))))

	gamma1_evidence <- sum_ML_over_PP(object[["y_log_lik_t_equals_1"]], object[["temperatures"]])
	gamma0_evidence <- gamma0_evidence(object[["y"]], tau0_shape=tau0_shape)
	lbf <- gamma1_evidence-gamma0_evidence

	conditional_prob_pathogenic <- apply(object[["z"]][,object[["k"]]*(length(object[["temperatures"]])-1)+seq(object[["k"]]),drop=FALSE], 2, mean)
	prob_gamma1 <- gamma1_prior * exp(lbf) / (gamma1_prior * exp(lbf) + 1 - gamma1_prior)

	structure(list(
		prob_gamma1=prob_gamma1,
		conditional_prob_pathogenic=conditional_prob_pathogenic,
		marginal_prob_pathogenic=prob_gamma1 * conditional_prob_pathogenic,
		phi=mean(exp(object[["log_phi"]])),
		omega=mean(1-1/(1+exp(object[["logit_omega"]]))),
		confidence_interval=CI_gamma1_evidence(
			temperatures=object[["temperatures"]],
			y_log_lik_t_equals_1_traces=object[["y_log_lik_t_equals_1"]],
			confidence=confidence,
			simulations=simulations
		),
		expected_explained=if (prod(dim(object$x)) > 0) extract_expected_explained(object) else NULL,
		explaining_variants=if (prod(dim(object$z)) > 0) extract_explaining_variants(object) else NULL,
		log_gamma1_evidence=gamma1_evidence,
		log_gamma0_evidence=gamma0_evidence,
		log_BF=lbf,
		gamma1_prior=gamma1_prior,
		number_of_posterior_samples=nrow(object[["y_log_lik_t_equals_1"]]),
		omega_estimated=object[["estimate_omega"]],
		phi_estimated=object[["estimate_phi"]],
		phi_acceptance_rate=apply(object[["log_phi"]], 2, function(log_phis) mean(log_phis[-length(log_phis)] != log_phis[-1])), 
		omega_acceptance_rate=apply(object[["logit_omega"]], 2, function(logit_omegas) mean(logit_omegas[-length(logit_omegas)] != logit_omegas[-1])), 
		n=length(object[["y"]]),
		k=object[["k"]],
		variant_counts=variant_counts,
		temperatures=object[["temperatures"]]
	), class="BeviMed_summary")
}

#' Print readable summary of \code{BeviMed_summary} object.
#'
#' @param x \code{BeviMed_summary} object.
#' @param print_prob_pathogenic Logical value indicating whether to print list of marginal probabilities of \code{z_j = 1} for all variants \code{j}.
#' @param ... Not-used arguments 
#' @return Prints a summary
#' @method print BeviMed_summary
#' @export
print.BeviMed_summary <- function(x, print_prob_pathogenic=TRUE, ...) {
	stopifnot(class(x) == "BeviMed_summary")
	dashed <- paste0(rep("-", getOption("width")), collapse="")
	cat(dashed, "\n")
	cat("The probability of association is ", round(x[["prob_gamma1"]], digits=2), " [prior: ", round(x[["gamma1_prior"]], digits=2), "]\n", sep="")
	if (!is.null(x$expected_explained)) cat("\nThe expected number of cases explained is: ", round(x$expected_explained, digits=2), sep="")
	if (!is.null(x$explaining_variants)) cat("\nThe expected number of variants involved in explained cases is: ", round(x$explaining_variants, digits=2), sep="")
	cat("\n")
	cat("\nLog Bayes factor between gamma 1 model and gamma 0 model is ", round(x[["log_BF"]], digits=2), sep="")
	cat("\nA confidence interval for the log Bayes factor is:\n")
	print(round(digits=2, x[["confidence_interval"]] - x[["log_gamma0_evidence"]]))

	cat(dashed, "\n")
	if (x[["omega_estimated"]]) {
		cat("Estimate of omega: ", round(digits=2, x[["omega"]]), "\n", sep="")
		cat("\tAcceptance rate in sequential chains: \n\t", paste0(collapse=" : ", round(digits=2, x[["omega_acceptance_rate"]])), "\n", sep="")
		if (x[["phi_estimated"]]) {
			cat("Estimate of phi: ", round(digits=2, x[["phi"]]), "\n", sep="")
			cat("\tAcceptance rate in sequential chains: \n\t", paste0(collapse=" : ", round(digits=2, x[["phi_acceptance_rate"]])), "\n", sep="")
		}
		cat(dashed, "\n")
	}

	#only show this for the highest temperature...
	if (print_prob_pathogenic) {
		cat("Estimated probabilities of pathogenicity of individual variants\n")
		cat("(conditional on gamma = 1)\n\n")

		print(row.names=FALSE, data.frame(
			check.names=FALSE,
			stringsAsFactors=FALSE,
			Variant=if (is.null(names(x[["conditional_prob_pathogenic"]]))) seq(length(x[["conditional_prob_pathogenic"]])) else names(x[["conditional_prob_pathogenic"]]),
			Controls=x[["variant_counts"]][["FALSE"]],
			Cases=x[["variant_counts"]][["TRUE"]],
			`P(z_j=1|y,gamma=1)`=round(digits=2, x[["conditional_prob_pathogenic"]]),
			`Bar Chart`=sapply(seq(length(x[["conditional_prob_pathogenic"]])), function(j) paste0("[", paste0(collapse="", rep("=", as.integer(x[["conditional_prob_pathogenic"]][j]*20))), paste0(collapse="", rep(" ", 20-as.integer(x[["conditional_prob_pathogenic"]][j]*20))), "]"))))

	}
}

#' Print readable summary of \code{BeviMed} object.
#'
#' @param x \code{BeviMed} object.
#' @param ... Arguments passed to \code{\link{summary.BeviMed}} 
#' @return Prints a summary
#' @method print BeviMed
#' @export
print.BeviMed <- function(x, ...) {
	stopifnot(class(x) == "BeviMed")
	print.BeviMed_summary(summary(x, ...))
}


