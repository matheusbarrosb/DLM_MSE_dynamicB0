get_spm_input = function(stan_data,
                         true_catch_new, 
                         true_nya_new,   
                         waa,
                         selectivity,
                         sigma_index = 0.05) {
  
  obs_catch_new = true_catch_new
  
  b_vuln_new = sum(true_nya_new * waa * selectivity)
  q_sim = 1e-3
  obs_index_new = (q_sim * b_vuln_new) * rlnorm(1, -0.5 * sigma_index^2, sigma_index)
  
  stan_data$catches = c(stan_data$catches, as.numeric(obs_catch_new))
  stan_data$index   = c(stan_data$index, as.numeric(obs_index_new))
  
  stan_data$nyears = as.integer(stan_data$nyears + 1)
  
  return(stan_data)
}