#include <Rcpp.h>
#include <cmath>
#include "sampling.h"

#ifndef INCLUDED_RCPP_EXPORTS_H
#define INCLUDED_RCPP_EXPORTS_H

using namespace Rcpp;
using namespace std;

RcppExport SEXP R_parallel_tempered_markov_chain(
	SEXP its, 
	SEXP y, 
	SEXP var_block_start_index, 
	SEXP var_block_stop_index, 
	SEXP cases, 
	SEXP counts, 
	SEXP min_ac, 
	SEXP null_shape1, 
	SEXP null_shape2, 
	SEXP patho_shape1, 
	SEXP patho_shape2, 
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
);

#endif

