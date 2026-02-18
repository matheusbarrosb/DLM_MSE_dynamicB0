get_lbspr_input = function(stan_data, current_obs_len) {
  
  # scale to effective sample size
  ess = 100
  obs_prop = current_obs_len / sum(current_obs_len)
  obs_counts_scaled = round(obs_prop * ess)
  
  # update observed counts
  stan_data$obs_counts = as.integer(obs_counts_scaled)
  
  return(stan_data)
}