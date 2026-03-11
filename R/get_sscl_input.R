get_sscl_input = function(stan_data, current_obs_catch, current_obs_len) {
  
  # scale new length observation to ESS
  ess = 100
  prop = current_obs_len / sum(current_obs_len)
  scaled_len = as.integer(round(prop * ess))
  
  # append data
  stan_data$catches = c(stan_data$catches, as.numeric(current_obs_catch))
  stan_data$obs_len_counts = rbind(stan_data$obs_len_counts, scaled_len)
  stan_data$nyears = as.integer(stan_data$nyears + 1)
  
  return(stan_data)
}