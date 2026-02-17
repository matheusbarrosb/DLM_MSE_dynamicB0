make_historical_sca_data = function(sim, burn_in_length, waa, maturity, selectivity, 
                                    nages, M = 0.2, sigma_index = 0.2, 
                                    sigma_catch = 0.15, ess = 100) {
  
  years = 1:burn_in_length
  
  obs_catch = as.numeric(sim$catch[years])
  nya_true  = sim$nya_mat[years, ]
  
  b_vuln_true = numeric(burn_in_length)
  for(y in 1:burn_in_length) {
    b_vuln_true[y] = sum(nya_true[y, ] * waa * selectivity)
  }
  
  q_sim = 1e-3 
  index_obs = as.numeric((q_sim * b_vuln_true) * rlnorm(burn_in_length, -0.5 * sigma_index^2, sigma_index))
  
  # Create matrix
  obs_nya = matrix(0, nrow = burn_in_length, ncol = nages)
  
  for(y in 1:burn_in_length) {
    n_vuln = nya_true[y, ] * selectivity
    prob_a = n_vuln / sum(n_vuln)
    prob_a[is.na(prob_a)] = 0
    obs_nya[y, ] = t(rmultinom(1, size = ess, prob = prob_a))
  }
  
  # CRITICAL: Force R matrix to be strictly integer type for Stan 'int' array
  storage.mode(obs_nya) <- "integer"
  
  waa_matrix = matrix(rep(waa, each = burn_in_length), nrow = burn_in_length, ncol = nages)
  
  stan_data = list(
    nyears      = as.integer(burn_in_length),
    nages       = as.integer(nages),
    maturity    = as.numeric(maturity),
    waa         = waa_matrix,
    M           = as.numeric(M),
    catches     = as.numeric(obs_catch),
    selectivity = as.numeric(selectivity),
    index       = as.numeric(index_obs),
    nya_obs     = obs_nya,   
    catch_sigma = as.numeric(sigma_catch),
    index_sigma = as.numeric(sigma_index)
  )
  
  return(stan_data)
}