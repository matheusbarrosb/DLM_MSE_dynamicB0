run_sca = function(stan_input, model_obj, mcmc_setup) {
  
  init_list = list(
    log_R0   = 7,
    h        = 0.7,
    q        = log(0.001),
    rec_devs = rep(0, stan_input$nyears - 1) 
  )
  
  if (mcmc_setup$chains > 1) {
    init_vals = lapply(1:mcmc_setup$chains, function(x) init_list)
  } else {
    init_vals = list(init_list)
  }
  
  fit = rstan::sampling(
    object  = model_obj,
    data    = stan_input,
    init    = init_vals,    
    chains  = mcmc_setup$chains,
    iter    = mcmc_setup$niter,
    warmup  = mcmc_setup$nwarmup,
    thin    = mcmc_setup$thin,
    refresh = mcmc_setup$verbose,
    control = list(adapt_delta = mcmc_setup$adapt_delta, max_treedepth = mcmc_setup$max_treedepth)
  )
  
  post = rstan::extract(fit)
  
  df_biomass = data.frame(
    year = 1:stan_input$nyears,
    mean = apply(post$SSB, 2, mean),
    sd   = apply(post$SSB, 2, sd)
  )
  
  post_N = exp(post$log_N) 
  
  mat_N_mean = apply(post_N, c(2,3), mean)
  mat_N_sd   = apply(post_N, c(2,3), sd)
  
  devs_mean = apply(post$rec_devs, 2, mean)
  devs_sd   = apply(post$rec_devs, 2, sd)
  
  df_rec_devs = data.frame(
    year = 2:stan_input$nyears,
    mean = devs_mean,
    sd   = devs_sd
  )
  df_rec_devs = rbind(data.frame(year=1, mean=0, sd=0), df_rec_devs)
  
  df_N_long = as.data.frame(mat_N_mean)
  colnames(df_N_long) = 1:stan_input$nages
  df_N_long$year = 1:stan_input$nyears
  
  df_N_long = tidyr::pivot_longer(
    df_N_long,
    cols      = -year,
    names_to  = "age",
    values_to = "N"
  )
  
  df_N_long$age = as.numeric(df_N_long$age)
  df_N_long$cohort = df_N_long$year - df_N_long$age
  
  # extract convergence information
  rhat_values = summary(fit)$summary[,"Rhat"]
  divergent = sum(rhat_values > 1.1)
  
  
  return(list(
    biomass     = df_biomass,
    N_mean      = mat_N_mean,
    N_sd        = mat_N_sd,
    rec_devs    = df_rec_devs,
    N_long      = df_N_long,
    rhat         = rhat_values,
    divergent    = divergent
  ))
}