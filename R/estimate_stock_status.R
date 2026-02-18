estimate_stock_status = function(method, model_state = NULL, current_obs = NULL, 
                                 history_data = NULL, mcmc_setup = NULL, model_objects = NULL) {
  
  output = list(
    status = NA,
    updated_state = NULL
  )
  
  if (method == "SCA") {
    
    if (is.null(model_state)) {
      stan_data = make_historical_sca_data(
        sim            = history_data$sim,
        waa            = history_data$waa,
        maturity       = history_data$maturity,
        selectivity    = history_data$selectivity,
        nages          = history_data$nages,
        burn_in_length = history_data$burn_in_length,
        sigma_index    = history_data$sigma_index,
        sigma_catch    = history_data$sigma_catch,
        ess            = history_data$ess
      )
    } else {
      stan_data = get_sca_input(
        stan_data      = model_state,
        true_catch_new = current_obs$catch,
        true_nya_new   = current_obs$nya,
        waa            = current_obs$waa,
        selectivity    = current_obs$selectivity,
        sigma_index    = current_obs$sigma_index,
        sigma_catch    = current_obs$sigma_catch,
        ess            = current_obs$ess
      )
    }
    
    fit = run_sca(stan_data, model_objects$SCA, mcmc_setup)
    
    est_status = tail(fit$biomass$mean, 1) / fit$biomass$mean[1]
    
    output$status        = est_status
    output$updated_state = stan_data
    
  } else if (method == "SPM") {
    
    if (is.null(model_state)) {
      stan_data = make_historical_spm_data(
        sim            = history_data$sim,
        burn_in_length = history_data$burn_in_length,
        waa            = history_data$waa,
        selectivity    = history_data$selectivity,
        sigma_index    = history_data$sigma_index
      )
    } else {
      stan_data = get_spm_input(
        stan_data      = model_state,
        true_catch_new = current_obs$catch,
        true_nya_new   = current_obs$nya,
        waa            = current_obs$waa,
        selectivity    = current_obs$selectivity,
        sigma_index    = current_obs$sigma_index
      )
    }
    
    init_list = list(
      log_r = log(0.3),
      log_k = log(150000),
      log_q = log(0.001),
      sigma_process = 0.1,
      sigma_obs = 0.1,
      log_B_devs = rep(0, stan_data$nyears)
    )
    
    fit = rstan::sampling(
      object  = model_objects$SPM,
      data    = stan_data,
      init    = function() init_list,
      chains  = mcmc_setup$chains,
      iter    = mcmc_setup$niter,
      warmup  = mcmc_setup$nwarmup,
      thin    = mcmc_setup$thin,
      refresh = mcmc_setup$verbose,
      control = list(adapt_delta = mcmc_setup$adapt_delta, max_treedepth = mcmc_setup$max_treedepth)
    )
    
    post = rstan::extract(fit)
    
    output$status        = median(post$B_B0) 
    output$updated_state = stan_data
    
  } else if (method == "LBSCA") {
    
    output$status        = 1.0 
    output$updated_state = model_state
    
  } else if (method == "LBSPR") {
    
    output$status        = 1.0 
    output$updated_state = model_state
    
  } else {
    stop("Unknown estimation method")
  }
  
  return(output)
}