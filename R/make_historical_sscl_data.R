make_historical_sscl_data = function(sim, burn_in_length, vb_params, 
                                     waa, maturity, nages, M, sr_params, bin_size = 1) {
  
  # subset historical period
  years = 1:burn_in_length
  obs_catch = as.numeric(sim$catch[years])
  obs_len_mat = sim$length_mat[years, ]
  
  # determine bin midpoints
  len_mids = as.numeric(colnames(sim$length_mat))
  if(any(is.na(len_mids))) {
    max_l = vb_params$linf * 1.5
    len_mids = seq(bin_size/2, max_l, by = bin_size)[1:ncol(sim$length_mat)]
  }
  
  # scale length compositions to ESS
  ess = 100
  obs_len_counts = matrix(0, nrow = burn_in_length, ncol = length(len_mids))
  
  for (y in 1:burn_in_length) {
    if (sum(obs_len_mat[y, ]) > 0) {
      prop = obs_len_mat[y, ] / sum(obs_len_mat[y, ])
      obs_len_counts[y, ] = as.integer(round(prop * ess))
    }
  }
  
  stan_data = list(
    nyears         = burn_in_length,
    n_ages         = nages,
    n_len          = length(len_mids),
    catches        = obs_catch,
    obs_len_counts = obs_len_counts,
    len_mids       = len_mids,
    
    M              = M,
    k              = vb_params$k,
    Linf           = vb_params$linf,
    t0             = vb_params$t0,
    CV_L           = 0.1,
    waa            = waa,
    mat            = maturity,
    
    h              = sr_params$h_high, # assuming high steepness as base
    R0             = sr_params$R0,
    SSB0           = sr_params$SSB0,
    sigma_R        = sr_params$sigma_rec,
    sigma_C        = 0.05
  )
  
  return(stan_data)
}