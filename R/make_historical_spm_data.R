make_historical_spm_data = function(sim, burn_in_length, waa, selectivity, 
                                    sigma_index = 0.2) {
  
  years = 1:burn_in_length
  
  obs_catch = as.numeric(sim$catch[years])
  nya_true  = sim$nya_mat[years, ]
  
  b_vuln_true = numeric(burn_in_length)
  for(y in 1:burn_in_length) {
    b_vuln_true[y] = sum(nya_true[y, ] * waa * selectivity)
  }
  
  q_sim = 1e-3 
  index_obs = as.numeric((q_sim * b_vuln_true) * rlnorm(burn_in_length, -0.5 * sigma_index^2, sigma_index))
  
  stan_data = list(
    nyears  = as.integer(burn_in_length),
    nages   = 1, # unused but required by data block
    catches = as.numeric(obs_catch),
    index   = as.numeric(index_obs)
  )
  
  return(stan_data)
}