#include <Rcpp.h>
#include <cmath>

#ifndef INCLUDED_SAMPLING_H
#define INCLUDED_SAMPLING_H

using namespace Rcpp;
using namespace std;

const double MINUS_LOG_SQRT_2_PI = -0.9189385;

inline double sq(double x) { return x * x; }

inline double log_likelihood_normal(
	double mean,
	double sd,
	double value
){ return MINUS_LOG_SQRT_2_PI - log(sd) - 0.5 * sq((value - mean)/sd); }

inline int random_integer(int exc_max)
{
	return (int)(unif_rand() * (double)exc_max) % exc_max;
}

inline double logit(double value) {
	return log(value) - log(1.0 - value);
}

inline double expit(double value) {
	return 1.0 - 1.0/(1.0+exp(value));
}

inline double logit_beta(
	double shape1,
	double shape2,
	double value
){ return value - 2.0 * log(1.0 + exp(value)) + log(R::dbeta(expit(value), shape1, shape2, true)); }

List bevimed_mc(
	int its,
	LogicalVector y,
	IntegerVector var_block_start_index,
	IntegerVector var_block_stop_index,
	IntegerVector cases,
	IntegerVector counts,
	IntegerVector min_ac,
	double q_shape1,
	double q_shape2,
	double p_shape1,
	double p_shape2,
	double Z_shape1,
	double Z_shape2,
	LogicalMatrix Z0,
	bool estimate_logit_Z_rate,
	NumericVector logit_Z_rates,
	NumericVector logit_Z_rate_proposal_sds,
	NumericVector Z_weights,
	bool estimate_phi,
	NumericVector log_phis,
	double log_phi_mean,
	double log_phi_sd,
	NumericVector log_phi_proposal_sds,
	NumericVector t,
	int swaps,
	bool annealing,
	int tandem_variant_updates,
	IntegerVector y1_case_block_start_index,
	IntegerVector y1_case_block_stop_index,
	IntegerVector y1_variants,
	bool store_Z_trace
);

#endif
