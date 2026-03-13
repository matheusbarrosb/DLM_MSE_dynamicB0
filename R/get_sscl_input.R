get_sscl_input = function(stan_data, current_obs_catch, current_obs_len) {
  
  if (current_obs_catch < 1e-4) current_obs_catch = 1e-4
  
  ess = 100
  if (sum(current_obs_len) > 0) {
    prop = current_obs_len / sum(current_obs_len)
    scaled_len = as.integer(round(prop * ess))
  } else {
    scaled_len = integer(length(current_obs_len))
  }
  
  # append data
  stan_data$catches = c(stan_data$catches, as.numeric(current_obs_catch))
  stan_data$obs_len_counts = rbind(stan_data$obs_len_counts, scaled_len)
  stan_data$nyears = as.integer(stan_data$nyears + 1)
  
  return(stan_data)
}