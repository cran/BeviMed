#include "sampling.h"

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
) {
	RNGScope scope;
	
	double logit_Z_rate_mean = Z_shape1;
	double logit_Z_rate_sd = Z_shape2;

	int n = y.length();

	int k = var_block_start_index.length();
	int num_temps = t.length();


	IntegerVector chain_temperature_reference(num_temps);
	IntegerVector temperature_chain_reference(num_temps);
	for (int temp = 0; temp < num_temps; temp++) {
		chain_temperature_reference[temp] = temp;
		temperature_chain_reference[temp] = temp;
	}

	LogicalMatrix Z_trace(store_Z_trace ? its : 0, store_Z_trace ? (k * num_temps) : 0);
	NumericMatrix y_log_lik_trace(its, num_temps);
	NumericMatrix y_log_lik_t_equals_1_trace(its, num_temps);
	LogicalMatrix Z(num_temps, k);

	NumericMatrix logit_Z_rates_trace(its, num_temps);
	NumericMatrix log_phis_trace(its, num_temps);

	LogicalVector swap_accept_trace(swaps * its);
	IntegerVector swap_temp1_trace(swaps * its);

	IntegerVector count_Z1(num_temps, 0);
	for (int temp = 0; temp < num_temps; temp++) {
		for (int v = 0; v < k; v++) {
			Z(temp, v) = Z0(temp, v); 
			count_Z1[temp] += (int)Z(temp, v);	
		}
	}

	IntegerMatrix pathogenic_var_count(num_temps, n);
	IntegerVector temporary_pathogenic_var_count(n, 0);
	LogicalVector temporary_counted_indicator(n, false);

	LogicalMatrix x(num_temps, n);
	for (int temp = 0; temp < num_temps; temp++) {
		for (int i = 0; i < n; i++) {
			pathogenic_var_count(temp, i) = 0;
			x(temp, i) = false;
		}
	}

	IntegerVector count_y1(num_temps, 0);
	IntegerVector count_x1(num_temps, 0);
	IntegerVector count_y1x1(num_temps, 0);

	for (int temp = 0; temp < num_temps; temp++) {
		for (int v = 0; v < k; v++) {
			for (int case_with_v = var_block_start_index[v]; case_with_v < var_block_stop_index[v]; case_with_v++) {
				pathogenic_var_count(temp, cases[case_with_v]) += (int)Z(temp, v) * counts[case_with_v];
			}
		}

		for (int i = 0; i < n; i++) {
			x(temp, i) = pathogenic_var_count(temp, i) >= min_ac[i];
			count_x1[temp] += (int)x(temp, i);
			if (y[i]) { 
				count_y1x1[temp] += (int)x(temp, i);
				count_y1[temp] += 1;
			}
		}
	}

	NumericVector q1_tab(n+1);
	NumericVector q2_tab(n+1);
	NumericVector p1_tab(n+1);
	NumericVector p2_tab(n+1);
	NumericVector q_tot_tab(n+1);
	NumericVector p_tot_tab(n+1);
	q1_tab[0] = 0.0;
	q2_tab[0] = 0.0;
	p1_tab[0] = 0.0;
	p2_tab[0] = 0.0;
	q_tot_tab[0] = 0.0;
	p_tot_tab[0] = 0.0;
	for (int i = 1; i <= n; i++) {
		q1_tab[i] = q1_tab[i-1] + log(i-1 + q_shape1);
		q2_tab[i] = q2_tab[i-1] + log(i-1 + q_shape2);
		p1_tab[i] = p1_tab[i-1] + log(i-1 + p_shape1);
		p2_tab[i] = p2_tab[i-1] + log(i-1 + p_shape2);
		q_tot_tab[i] = q_tot_tab[i-1] + log(i-1 + q_shape1 + q_shape2);
		p_tot_tab[i] = p_tot_tab[i-1] + log(i-1 + p_shape1 + p_shape2);
	}

	NumericVector Z_shape1_gamma_table(k+1);
	NumericVector Z_shape2_gamma_table(k+1);
	NumericVector Z_total_gamma_table(k+1);
	Z_shape1_gamma_table[0] = 0.0;
	Z_shape2_gamma_table[0] = 0.0;
	Z_total_gamma_table[0] = 0.0;
	for (int v = 0; v <= k; v++) {
		Z_shape1_gamma_table[v] = Z_shape1_gamma_table[v-1] + log(v-1 + Z_shape1);
		Z_shape2_gamma_table[v] = Z_shape2_gamma_table[v-1] + log(v-1 + Z_shape2);
		Z_total_gamma_table[v] = Z_total_gamma_table[v-1] + log(v-1 + Z_shape1 + Z_shape2);
	}

	for (int it = 0; it < its; it++) {
		double annealing_factor = annealing ? (double)(its-it) : 1.0;
		for (int chain_number = 0; chain_number < num_temps; chain_number++) {
			if (estimate_logit_Z_rate) {
				double proposal = logit_Z_rates[chain_number] + norm_rand() * logit_Z_rate_proposal_sds[chain_number];

				double ll_cur = logit_beta(logit_Z_rate_mean, logit_Z_rate_sd, logit_Z_rates[chain_number]);
				double ll_upd = logit_beta(logit_Z_rate_mean, logit_Z_rate_sd, proposal);
				for (int v = 0; v < k; v++) {
					double pv_cur = expit(exp(log_phis[chain_number]) * Z_weights[v] + logit_Z_rates[chain_number]);
					double pv_upd = expit(exp(log_phis[chain_number]) * Z_weights[v] + proposal);
					ll_cur += Z(chain_number, v) ? log(pv_cur) : log(1.0-pv_cur);
					ll_upd += Z(chain_number, v) ? log(pv_upd) : log(1.0-pv_upd);
				}
				if (log(unif_rand()) < (ll_upd-ll_cur)) {
					logit_Z_rates[chain_number] = proposal;
				}

				if (estimate_phi) {
					double phi_proposal = log_phis[chain_number] + norm_rand() * log_phi_proposal_sds[chain_number];

					double ll_phi_cur = log_likelihood_normal(log_phi_mean, log_phi_sd, log_phis[chain_number]);
					double ll_phi_upd = log_likelihood_normal(log_phi_mean, log_phi_sd, phi_proposal);
					for (int v = 0; v < k; v++) {
						double pv_cur = expit(exp(log_phis[chain_number]) * Z_weights[v] + logit_Z_rates[chain_number]);
						double pv_upd = expit(exp(phi_proposal) * Z_weights[v] + logit_Z_rates[chain_number]);
						ll_phi_cur += Z(chain_number, v) ? log(pv_cur) : log(1.0-pv_cur);
						ll_phi_upd += Z(chain_number, v) ? log(pv_upd) : log(1.0-pv_upd);
					}
					if (log(unif_rand()) < (ll_phi_upd-ll_phi_cur)) {
						log_phis[chain_number] = phi_proposal;
					}

				}
			}
			for (int v = 0; v < k; v++) {
				int dx1 = 0;
				int dy1x1 = 0;

				int c_s01 = count_y1[chain_number] - count_y1x1[chain_number];
				int c_s00 = (n - count_x1[chain_number]) - (count_y1[chain_number] - count_y1x1[chain_number]);
				int c_s11 = count_y1x1[chain_number];
				int c_s10 = count_x1[chain_number] - count_y1x1[chain_number];

				for (int i = var_block_start_index[v]; i < var_block_stop_index[v]; i++) {
					int alt_pathogenic_var_count = pathogenic_var_count(chain_number, cases[i]) + (Z(chain_number, v) ? (-counts[i]) : counts[i]);
					bool alt_x = alt_pathogenic_var_count >= min_ac[cases[i]];
					if (alt_x != x(chain_number, cases[i])) {
						dx1 += (int)alt_x - (int)x(chain_number, cases[i]);
						if (y[cases[i]]) dy1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
					}
				}

				int alt_count_x1 = count_x1[chain_number] + dx1;
				int alt_count_y1x1 = count_y1x1[chain_number] + dy1x1;

				int u_s01 = count_y1[chain_number] - alt_count_y1x1;
				int u_s00 = (n - alt_count_x1) - (count_y1[chain_number] - alt_count_y1x1);
				int u_s11 = alt_count_y1x1;
				int u_s10 = alt_count_x1 - alt_count_y1x1;

				double Z_contribution;
				if (estimate_logit_Z_rate) {
					int dZ1 = Z(chain_number, v) ? -1 : 1;
					double vp = expit(logit_Z_rates[chain_number] + exp(log_phis[chain_number]) * Z_weights[v]);
					Z_contribution = + dZ1 * log(1.0-vp) - dZ1 * log(vp); 
				} 
				else {
					Z_contribution = + ((Z(chain_number, v)) ? (log(count_Z1[chain_number] + Z_shape1 - 1.0) - log(k - count_Z1[chain_number] + Z_shape2)) : (-log(count_Z1[chain_number] + Z_shape1) + log(k - count_Z1[chain_number] + Z_shape2 - 1.0)));
				}

				double current_Zv_log_odds = 
					+ Z_contribution
					+ (
						+ (q1_tab[c_s01]-q1_tab[u_s01])
						+ (q2_tab[c_s00]-q2_tab[u_s00])
						+ (p1_tab[c_s11]-p1_tab[u_s11])
						+ (p2_tab[c_s10]-p2_tab[u_s10])
						+ (q_tot_tab[u_s01 + u_s00]-q_tot_tab[c_s01 + c_s00])
						+ (p_tot_tab[u_s11 + u_s10]-p_tot_tab[c_s11 + c_s10])
					) * t[chain_temperature_reference[chain_number]] * annealing_factor
				;

				if (unif_rand() < (1.0-1.0/(1.0+exp(-current_Zv_log_odds)))) {
					count_x1[chain_number] += dx1;
					count_y1x1[chain_number] += dy1x1;
					for (int i = var_block_start_index[v]; i < var_block_stop_index[v]; i++) {
						pathogenic_var_count(chain_number, cases[i]) += (Z(chain_number, v) ? (-counts[i]) : counts[i]);
						x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
					}
					count_Z1[chain_number] += Z(chain_number, v) ? (-1) : 1;
					Z(chain_number, v) = !Z(chain_number, v);
				}
			}
		}

		if (tandem_variant_updates > 0) {
			for (int chain_number = 0; chain_number < num_temps; chain_number++) {
				for (int double_variant_update_number = 0; double_variant_update_number < tandem_variant_updates; double_variant_update_number ++) {
					


					int v1 = y1_variants[random_integer(y1_variants.length())];
					int v2;
					do {
						v2 = y1_variants[random_integer(y1_variants.length())];
					} while (v1 == v2);


					int no_change_s01 = count_y1[chain_number] - count_y1x1[chain_number];
					int no_change_s00 = (n - count_x1[chain_number]) - (count_y1[chain_number] - count_y1x1[chain_number]);
					int no_change_s11 = count_y1x1[chain_number];
					int no_change_s10 = count_x1[chain_number] - count_y1x1[chain_number];


					int v1_change_count_x1 = count_x1[chain_number];
					int v1_change_count_y1x1 = count_y1x1[chain_number];

					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
						int alt_pathogenic_var_count = pathogenic_var_count(chain_number, cases[i]) + (Z(chain_number, v1) ? (-counts[i]) : counts[i]);
						bool alt_x = alt_pathogenic_var_count >= min_ac[cases[i]];
						if (alt_x != x(chain_number, cases[i])) {
							v1_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v1_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
					}

					int v1_change_s01 = count_y1[chain_number] - v1_change_count_y1x1;
					int v1_change_s00 = (n - v1_change_count_x1) - (count_y1[chain_number] - v1_change_count_y1x1);
					int v1_change_s11 = v1_change_count_y1x1;
					int v1_change_s10 = v1_change_count_x1 - v1_change_count_y1x1;

					int v2_change_count_x1 = count_x1[chain_number];
					int v2_change_count_y1x1 = count_y1x1[chain_number];

					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
						int alt_pathogenic_var_count = pathogenic_var_count(chain_number, cases[i]) + (Z(chain_number, v2) ? (-counts[i]) : counts[i]);
						bool alt_x = alt_pathogenic_var_count >= min_ac[cases[i]];
						if (alt_x != x(chain_number, cases[i])) {
							v2_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v2_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
					}

					int v2_change_s01 = count_y1[chain_number] - v2_change_count_y1x1;
					int v2_change_s00 = (n - v2_change_count_x1) - (count_y1[chain_number] - v2_change_count_y1x1);
					int v2_change_s11 = v2_change_count_y1x1;
					int v2_change_s10 = v2_change_count_x1 - v2_change_count_y1x1;

					int v1_and_v2_change_count_x1 = count_x1[chain_number];
					int v1_and_v2_change_count_y1x1 = count_y1x1[chain_number];
					
					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++)
						temporary_pathogenic_var_count[cases[i]] = pathogenic_var_count(chain_number, cases[i]);
					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++)
						temporary_pathogenic_var_count[cases[i]] = pathogenic_var_count(chain_number, cases[i]);
					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++)
						temporary_pathogenic_var_count[cases[i]] += (Z(chain_number, v1) ? (-counts[i]) : counts[i]);
					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++)
						temporary_pathogenic_var_count[cases[i]] += (Z(chain_number, v2) ? (-counts[i]) : counts[i]);

					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
						bool alt_x = temporary_pathogenic_var_count[cases[i]] >= min_ac[cases[i]];
						if ((!temporary_counted_indicator[cases[i]]) && (alt_x != x(chain_number, cases[i]))) {
							v1_and_v2_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v1_and_v2_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
						temporary_counted_indicator[cases[i]] = true;
					}

					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
						bool alt_x = temporary_pathogenic_var_count[cases[i]] >= min_ac[cases[i]];
						if ((!temporary_counted_indicator[cases[i]]) && (alt_x != x(chain_number, cases[i]))) {
							v1_and_v2_change_count_x1 += (int)alt_x - (int)x(chain_number, cases[i]);
							if (y[cases[i]]) v1_and_v2_change_count_y1x1 += (int)alt_x - (int)x(chain_number, cases[i]);
						}
						temporary_counted_indicator[cases[i]] = true;
					}

					int v1_and_v2_change_s01 = count_y1[chain_number] - v1_and_v2_change_count_y1x1;
					int v1_and_v2_change_s00 = (n - v1_and_v2_change_count_x1) - (count_y1[chain_number] - v1_and_v2_change_count_y1x1);
					int v1_and_v2_change_s11 = v1_and_v2_change_count_y1x1;
					int v1_and_v2_change_s10 = v1_and_v2_change_count_x1 - v1_and_v2_change_count_y1x1;

					double v1_Z_log_odds_favour_current;
					if (estimate_logit_Z_rate) {
						int dZ1 = Z(chain_number, v1) ? -1 : 1;
						double vp = expit(logit_Z_rates[chain_number] + exp(log_phis[chain_number]) * Z_weights[v1]);
						v1_Z_log_odds_favour_current = + dZ1 * log(1.0-vp) - dZ1 * log(vp); 
					}
					else {
						v1_Z_log_odds_favour_current = ((Z(chain_number, v1)) ? (log(count_Z1[chain_number] + Z_shape1 - 1.0) - log(k - count_Z1[chain_number] + Z_shape2)) : (-log(count_Z1[chain_number] + Z_shape1) + log(k - count_Z1[chain_number] + Z_shape2 - 1.0)));
					}
					double v1_change_odds_ratio = exp((
						+ v1_Z_log_odds_favour_current 
						+ (
							+ (q1_tab[no_change_s01]-q1_tab[v1_change_s01])
							+ (q2_tab[no_change_s00]-q2_tab[v1_change_s00])
							+ (p1_tab[no_change_s11]-p1_tab[v1_change_s11])
							+ (p2_tab[no_change_s10]-p2_tab[v1_change_s10])
							+ (q_tot_tab[v1_change_s01 + v1_change_s00]-q_tot_tab[no_change_s01 + no_change_s00])
							+ (p_tot_tab[v1_change_s11 + v1_change_s10]-p_tot_tab[no_change_s11 + no_change_s10])
						) * t[chain_temperature_reference[chain_number]] * annealing_factor
					) * (-1.0));

					double v2_Z_log_odds_favour_current;
					if (estimate_logit_Z_rate) {
						int dZ1 = Z(chain_number, v2) ? -1 : 1;
						double vp = expit(logit_Z_rates[chain_number] + exp(log_phis[chain_number]) * Z_weights[v2]);
						v2_Z_log_odds_favour_current = + dZ1 * log(1.0-vp) - dZ1 * log(vp); 
					}
					else {
						v2_Z_log_odds_favour_current = ((Z(chain_number, v2)) ? (log(count_Z1[chain_number] + Z_shape1 - 1.0) - log(k - count_Z1[chain_number] + Z_shape2)) : (-log(count_Z1[chain_number] + Z_shape1) + log(k - count_Z1[chain_number] + Z_shape2 - 1.0)));
					}
					double v2_change_odds_ratio = exp((
						+ v2_Z_log_odds_favour_current
						+ (
							+ (q1_tab[no_change_s01]-q1_tab[v2_change_s01])
							+ (q2_tab[no_change_s00]-q2_tab[v2_change_s00])
							+ (p1_tab[no_change_s11]-p1_tab[v2_change_s11])
							+ (p2_tab[no_change_s10]-p2_tab[v2_change_s10])
							+ (q_tot_tab[v2_change_s01 + v2_change_s00]-q_tot_tab[no_change_s01 + no_change_s00])
							+ (p_tot_tab[v2_change_s11 + v2_change_s10]-p_tot_tab[no_change_s11 + no_change_s10])
						) * t[chain_temperature_reference[chain_number]] * annealing_factor
					) * (-1.0));



					double v1_and_v2_Z_log_odds_favour_current;
					if (estimate_logit_Z_rate) {
						int dZ1_v1 = Z(chain_number, v1) ? -1 : 1;
						int dZ1_v2 = Z(chain_number, v2) ? -1 : 1;
						double v1p = expit(logit_Z_rates[chain_number] + exp(log_phis[chain_number]) * Z_weights[v1]);
						double v2p = expit(logit_Z_rates[chain_number] + exp(log_phis[chain_number]) * Z_weights[v2]);
						v1_and_v2_Z_log_odds_favour_current = 
							+ dZ1_v1 * log(1.0-v1p) - dZ1_v1 * log(v1p)
							+ dZ1_v2 * log(1.0-v2p) - dZ1_v2 * log(v2p)
						;
					}
					else {
						if (Z(chain_number, v1) != Z(chain_number, v2)) {
							v1_and_v2_Z_log_odds_favour_current = 0;
						}
						else {
							if (Z(chain_number, v1)) {
								v1_and_v2_Z_log_odds_favour_current = (
									+ log(Z_shape1 + count_Z1[chain_number] - 2.0)
									+ log(Z_shape1 + count_Z1[chain_number] - 1.0)
									- log(Z_shape2 + k - count_Z1[chain_number] + 1.0)
									- log(Z_shape2 + k - count_Z1[chain_number])
								);
							}
							else {
								v1_and_v2_Z_log_odds_favour_current = (
									- log(Z_shape1 + count_Z1[chain_number] + 1.0) 
									- log(Z_shape1 + count_Z1[chain_number])
									+ log(Z_shape2 + k - count_Z1[chain_number] - 1.0)
									+ log(Z_shape2 + k - count_Z1[chain_number] - 2.0)
								);
							}
						}
					}
					double v1_and_v2_change_odds_ratio = exp((
						+ v1_and_v2_Z_log_odds_favour_current
						+ (
							+ (q1_tab[no_change_s01]-q1_tab[v1_and_v2_change_s01])
							+ (q2_tab[no_change_s00]-q2_tab[v1_and_v2_change_s00])
							+ (p1_tab[no_change_s11]-p1_tab[v1_and_v2_change_s11])
							+ (p2_tab[no_change_s10]-p2_tab[v1_and_v2_change_s10])
							+ (q_tot_tab[v1_and_v2_change_s01 + v1_and_v2_change_s00]-q_tot_tab[no_change_s01 + no_change_s00])
							+ (p_tot_tab[v1_and_v2_change_s11 + v1_and_v2_change_s10]-p_tot_tab[no_change_s11 + no_change_s10])
						) * t[chain_temperature_reference[chain_number]] * annealing_factor
					) * (-1.0));

					double sum_ratios = 
						+ v1_change_odds_ratio
						+ v2_change_odds_ratio
						+ v1_and_v2_change_odds_ratio
					;

					double p_no_change = 1.0/(1.0+sum_ratios);
					double p_v1_change = p_no_change * v1_change_odds_ratio;
					double p_v2_change = p_no_change * v2_change_odds_ratio;


					double p_v2_threshold = p_no_change + p_v1_change + p_v2_change;
					double random_draw = unif_rand();



					if (random_draw >= p_no_change) {
						if (random_draw < (p_no_change + p_v1_change)) {
							count_x1[chain_number] = v1_change_count_x1;
							count_y1x1[chain_number] = v1_change_count_y1x1;
							for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
								pathogenic_var_count(chain_number, cases[i]) += (Z(chain_number, v1) ? (-counts[i]) : counts[i]);
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_Z1[chain_number] += Z(chain_number, v1) ? (-1) : 1;
							Z(chain_number, v1) = !Z(chain_number, v1);
						}
						else if (random_draw < p_v2_threshold) {
							count_x1[chain_number] = v2_change_count_x1;
							count_y1x1[chain_number] = v2_change_count_y1x1;
							for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
								pathogenic_var_count(chain_number, cases[i]) += (Z(chain_number, v2) ? (-counts[i]) : counts[i]);
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_Z1[chain_number] += Z(chain_number, v2) ? (-1) : 1;
							Z(chain_number, v2) = !Z(chain_number, v2);
						}
						else {
							count_x1[chain_number] = v1_and_v2_change_count_x1;
							count_y1x1[chain_number] = v1_and_v2_change_count_y1x1;
							for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
								pathogenic_var_count(chain_number, cases[i]) = temporary_pathogenic_var_count[cases[i]];
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_Z1[chain_number] += Z(chain_number, v1) ? (-1) : 1;
							Z(chain_number, v1) = !Z(chain_number, v1);
							for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
								pathogenic_var_count(chain_number, cases[i]) = temporary_pathogenic_var_count[cases[i]];
								x(chain_number, cases[i]) = pathogenic_var_count(chain_number, cases[i]) >= min_ac[cases[i]];
							}
							count_Z1[chain_number] += Z(chain_number, v2) ? (-1) : 1;
							Z(chain_number, v2) = !Z(chain_number, v2);
						}
					}

					for (int i = var_block_start_index[v1]; i < var_block_stop_index[v1]; i++) {
						temporary_counted_indicator[cases[i]] = false;
						temporary_pathogenic_var_count[cases[i]] = 0;
					}
					for (int i = var_block_start_index[v2]; i < var_block_stop_index[v2]; i++) {
						temporary_counted_indicator[cases[i]] = false;
						temporary_pathogenic_var_count[cases[i]] = 0;
					}
				}
			}
		}

		if ((swaps > 0) && (num_temps > 1)) {
			for (int swap_no = 0; swap_no < swaps; swap_no++) {
				int temperature_number1 = random_integer(num_temps-1);
				int chain1 = temperature_chain_reference[temperature_number1];
				int chain2 = temperature_chain_reference[temperature_number1+1];

				int s01_1 = count_y1[chain1] - count_y1x1[chain1];
				int s00_1 = (n - count_x1[chain1]) - (count_y1[chain1] - count_y1x1[chain1]);
				int s11_1 = count_y1x1[chain1];
				int s10_1 = count_x1[chain1] - count_y1x1[chain1];

				int s01_2 = count_y1[chain2] - count_y1x1[chain2];
				int s00_2 = (n - count_x1[chain2]) - (count_y1[chain2] - count_y1x1[chain2]);
				int s11_2 = count_y1x1[chain2];
				int s10_2 = count_x1[chain2] - count_y1x1[chain2];


				double chain1_y_log_lik_t_equals_1 =
					+ (q1_tab[s01_1]-q1_tab[0])
					+ (q2_tab[s00_1]-q2_tab[0])
					+ (p1_tab[s11_1]-p1_tab[0])
					+ (p2_tab[s10_1]-p2_tab[0])
					+ (q_tot_tab[0 + 0]-q_tot_tab[s01_1 + s00_1])
					+ (p_tot_tab[0 + 0]-p_tot_tab[s11_1 + s10_1])
				;
				
				double chain2_y_log_lik_t_equals_1 =
					+ (q1_tab[s01_2]-q1_tab[0])
					+ (q2_tab[s00_2]-q2_tab[0])
					+ (p1_tab[s11_2]-p1_tab[0])
					+ (p2_tab[s10_2]-p2_tab[0])
					+ (q_tot_tab[0 + 0]-q_tot_tab[s01_2 + s00_2])
					+ (p_tot_tab[0 + 0]-p_tot_tab[s11_2 + s10_2])
				;

				
				double logA = (
					- (t[chain_temperature_reference[chain1]] * chain1_y_log_lik_t_equals_1 + t[chain_temperature_reference[chain2]] * chain2_y_log_lik_t_equals_1)
					+ (t[chain_temperature_reference[chain2]] * chain1_y_log_lik_t_equals_1 + t[chain_temperature_reference[chain1]] * chain2_y_log_lik_t_equals_1)
				) * annealing_factor;

				swap_temp1_trace(it * swaps + swap_no) = chain_temperature_reference[chain1]+1;

				if (log(unif_rand()) < logA) {
					int temp_ref1 = chain_temperature_reference[chain1];
					chain_temperature_reference[chain1] = chain_temperature_reference[chain2];
					chain_temperature_reference[chain2] = temp_ref1;

					temperature_chain_reference[chain_temperature_reference[chain2]] = chain2;
					temperature_chain_reference[chain_temperature_reference[chain1]] = chain1;

					swap_accept_trace(it * swaps + swap_no) = true;

				} else {
					swap_accept_trace(it * swaps + swap_no) = false;
				}

			}
		}

		for (int chain_number = 0; chain_number < num_temps; chain_number++) {
			double y_log_lik_t_equals_1;
			
			int s01 = count_y1[chain_number] - count_y1x1[chain_number];
			int s00 = (n - count_x1[chain_number]) - (count_y1[chain_number] - count_y1x1[chain_number]);
			int s11 = count_y1x1[chain_number];
			int s10 = count_x1[chain_number] - count_y1x1[chain_number];

			y_log_lik_t_equals_1 =
				+ (q1_tab[s01]-q1_tab[0])
				+ (q2_tab[s00]-q2_tab[0])
				+ (p1_tab[s11]-p1_tab[0])
				+ (p2_tab[s10]-p2_tab[0])
				+ (q_tot_tab[0 + 0]-q_tot_tab[s01 + s00])
				+ (p_tot_tab[0 + 0]-p_tot_tab[s11 + s10])
			;

			double y_log_lik = t[chain_temperature_reference[chain_number]] * y_log_lik_t_equals_1 * annealing_factor;

			y_log_lik_trace(it, chain_temperature_reference[chain_number]) = y_log_lik;
			y_log_lik_t_equals_1_trace(it, chain_temperature_reference[chain_number]) = y_log_lik_t_equals_1;
			logit_Z_rates_trace(it, chain_temperature_reference[chain_number]) = logit_Z_rates[chain_number];
			log_phis_trace(it, chain_temperature_reference[chain_number]) = log_phis[chain_number];

			if (store_Z_trace) {
				for (int v = 0; v < k; v++) {
					Z_trace(it, chain_temperature_reference[chain_number] * k + v) = Z(chain_number, v);
				}
			}
		}
	}

	LogicalMatrix terminal_Z(num_temps, k);
	for (int chain_number = 0; chain_number < num_temps; chain_number++)
		for (int j = 0; j < k; j++)
			terminal_Z(chain_temperature_reference[chain_number], j) = Z(chain_number, j);

	NumericVector terminal_log_phi(num_temps);
	NumericVector terminal_logit_omega(num_temps);
	for (int chain_number = 0; chain_number < num_temps; chain_number++) {
		terminal_log_phi[chain_number] = log_phis_trace(its-1, chain_number);
		terminal_logit_omega[chain_number] = logit_Z_rates_trace(its-1, chain_number);
	}

	return List::create(
		Named("y_log_lik")=y_log_lik_trace,
		Named("y_log_lik_t_equals_1")=y_log_lik_t_equals_1_trace,
		Named("terminal_Z")=terminal_Z,
		Named("terminal_log_phi")=terminal_log_phi,
		Named("terminal_logit_omega")=terminal_logit_omega,
		Named("Z")=Z_trace,
		Named("swap_accept")=swap_accept_trace,
		Named("swap_at_temperature")=swap_temp1_trace,
		Named("logit_omega")=logit_Z_rates_trace,
		Named("log_phi")=log_phis_trace
	);
}
