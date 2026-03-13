get_lbspr_input = function(stan_data, current_obs_len) {
  
  ess = 100
  if (sum(current_obs_len) > 0) {
    obs_prop = current_obs_len / sum(current_obs_len)
    obs_counts_scaled = as.integer(round(obs_prop * ess))
  } else {
    obs_counts_scaled = integer(length(current_obs_len))
  }
  
  stan_data$obs_counts = obs_counts_scaled
  return(stan_data)
}