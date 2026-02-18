make_historical_lbspr_data = function(sim, burn_in_length, vb_params, 
                                      waa, maturity, selectivity, nages, M, bin_size = 1) {
  
  # get length composition from last year of burn-in
  obs_len_comp = sim$length_mat[burn_in_length, ]
  
  # determine bin midpoints
  len_mids = as.numeric(colnames(sim$length_mat))
  
  if(any(is.na(len_mids))) {
    max_l = vb_params$linf * 1.5
    len_mids = seq(bin_size/2, max_l, by = bin_size)[1:ncol(sim$length_mat)]
  }
  
  # scale to effective sample size
  ess = 100
  obs_prop = obs_len_comp / sum(obs_len_comp)
  obs_counts_scaled = round(obs_prop * ess)
  
  stan_data = list(
    n_len             = length(len_mids),
    n_ages            = nages,
    len_mids          = len_mids,
    obs_counts        = as.integer(obs_counts_scaled), 
    M                 = M,
    k                 = vb_params$k,
    Linf              = vb_params$linf,
    CV_L              = 0.1, 
    maturity_at_age   = maturity,
    weight_at_age     = waa,
    fixed_selectivity = selectivity
  )
  
  return(stan_data)
}