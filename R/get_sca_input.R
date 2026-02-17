get_sca_input = function(stan_data,
                         true_catch_new, 
                         true_nya_new,   
                         waa,
                         selectivity,
                         sigma_index = 0.05,
                         sigma_catch = 0.15,
                         ess = 100) {
  
  # --- 1. Simulate New Observations ---
  
  if (true_catch_new > 0) {
    obs_catch_new = true_catch_new
  } else {
    obs_catch_new = 0
  }
  
  b_vuln_new = sum(true_nya_new * waa * selectivity)
  q_sim = 1e-3
  obs_index_new = (q_sim * b_vuln_new) * rlnorm(1, -0.5 * sigma_index^2, sigma_index)
  
  n_vuln_new = true_nya_new * selectivity
  
  if (sum(n_vuln_new) > 0) {
    prob_a = n_vuln_new / sum(n_vuln_new)
    prob_a[is.na(prob_a)] = 0
    
    obs_nya_new = as.integer(rmultinom(1, size = ess, prob = prob_a))
  } else {
    obs_nya_new = rep(0L, length(selectivity)) # 0L ensures integer type
  }
  
  # --- 2. Update Stan Input Data ---
  
  # Append new observations
  stan_data$catches = as.numeric(c(stan_data$catches, obs_catch_new))
  stan_data$index   = as.numeric(c(stan_data$index, obs_index_new))
  
  # CRITICAL CHANGE 2: rbind upgrades types to numeric (double). 
  # We must bind, then FORCE IT BACK to integer.
  stan_data$nya_obs = rbind(stan_data$nya_obs, obs_nya_new)
  storage.mode(stan_data$nya_obs) <- "integer"
  
  # Increment dimensions (Force Integer)
  stan_data$nyears = as.integer(stan_data$nyears + 1)
  
  # Expand weight-at-age matrix
  stan_data$waa = rbind(stan_data$waa, waa)
  
  return(stan_data)
}