#include "RcppExports.h"

using namespace Rcpp;

RcppExport SEXP R_bevimed_mc(
	SEXP its, 
	SEXP y, 
	SEXP var_block_start_index, 
	SEXP var_block_stop_index, 
	SEXP cases, 
	SEXP counts, 
	SEXP min_ac, 
	SEXP q_shape1, 
	SEXP q_shape2, 
	SEXP p_shape1, 
	SEXP p_shape2, 
	SEXP Z_shape1, 
	SEXP Z_shape2, 
	SEXP Z0,
	SEXP estimate_logit_Z_rate,
	SEXP logit_Z_rates,
	SEXP logit_Z_rate_proposal_sds,
	SEXP Z_weights,
	SEXP estimate_phi,
	SEXP log_phis,
	SEXP log_phi_mean,
	SEXP log_phi_sd,
	SEXP log_phi_proposal_sds,
	SEXP t,
	SEXP swaps,
	SEXP annealing,
	SEXP tandem_variant_updates,
	SEXP y1_case_block_start_index,
	SEXP y1_case_block_stop_index,
	SEXP y1_variants,
	SEXP store_Z_trace
) {
BEGIN_RCPP
	return bevimed_mc(
		as<int>(its),
		as<LogicalVector>(y),
		as<IntegerVector>(var_block_start_index),
		as<IntegerVector>(var_block_stop_index),
		as<IntegerVector>(cases),
		as<IntegerVector>(counts),
		as<IntegerVector>(min_ac),
		as<double>(q_shape1),
		as<double>(q_shape2),
		as<double>(p_shape1),
		as<double>(p_shape2),
		as<double>(Z_shape1),
		as<double>(Z_shape2),
		as<LogicalMatrix>(Z0),
		as<bool>(estimate_logit_Z_rate),
		as<NumericVector>(logit_Z_rates),
		as<NumericVector>(logit_Z_rate_proposal_sds),
		as<NumericVector>(Z_weights),
		as<bool>(estimate_phi),
		as<NumericVector>(log_phis),
		as<double>(log_phi_mean),
		as<double>(log_phi_sd),
		as<NumericVector>(log_phi_proposal_sds),
		as<NumericVector>(t),
		as<int>(swaps),
		as<bool>(annealing),
		as<int>(tandem_variant_updates),
		as<IntegerVector>(y1_case_block_start_index),
		as<IntegerVector>(y1_case_block_stop_index),
		as<IntegerVector>(y1_variants),
		as<bool>(store_Z_trace)
	);
END_RCPP
}


