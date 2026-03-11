run_mse = function(n_sims, nyears, nages, init_nya,
                   waa, max_harvest_rate, rec_regime_length,
                   rec_type, sr_params, selectivity,
                   survival_mean, survival_sd, mcmc_setup,
                   burn_in_length, hist_harvest_rate, 
                   sca_model_path = NULL, spm_model_path = NULL, 
                   lbspr_model_path = NULL, sscl_model_path = NULL,
                   maturity, thresholds, vb_params, plot = FALSE,
                   estimation = FALSE, est_method = "SCA", bin_size = 2,
                   hcr_type = "absolute_hockey_stick",
                   scenario_name = "Default", 
                   parallel = FALSE) {
  
  log_file = paste0(here::here(), "/mse_progress_log.txt")
  cat(sprintf("\n=== STARTING BATCH: Scenario %s ===\n", scenario_name), file = log_file, append = TRUE)
  
  model_objects = list()
  
  if (estimation == TRUE) {
    if (est_method == "SCA") {
      if (is.null(sca_model_path)) stop("Estimation is TRUE (SCA) but 'sca_model_path' is not provided.")
      cat("Compiling SCA Stan model...\n", file = log_file, append = TRUE)
      message("Compiling SCA Stan model inside run_mse...")
      model_objects$SCA = rstan::stan_model(file = sca_model_path)
    }
    
    if (est_method == "SPM") {
      if (is.null(spm_model_path)) stop("Estimation is TRUE (SPM) but 'spm_model_path' is not provided.")
      cat("Compiling SPM Stan model...\n", file = log_file, append = TRUE)
      message("Compiling SPM Stan model inside run_mse...")
      model_objects$SPM = rstan::stan_model(file = spm_model_path)
    }
    
    if (est_method == "LBSPR") {
      if (is.null(lbspr_model_path)) stop("Estimation is TRUE (LBSPR) but 'lbspr_model_path' is not provided.")
      cat("Compiling LBSPR Stan model...\n", file = log_file, append = TRUE)
      message("Compiling LBSPR Stan model inside run_mse...")
      model_objects$LBSPR = rstan::stan_model(file = lbspr_model_path)
    }
    
    if (est_method == "SS_CL") {
      if (is.null(sscl_model_path)) stop("Estimation is TRUE (SS_CL) but 'sscl_model_path' is not provided.")
      cat("Compiling SS_CL Stan model...\n", file = log_file, append = TRUE)
      message("Compiling SS_CL Stan model inside run_mse...")
      model_objects$SSCL = rstan::stan_model(file = sscl_model_path)
    }
  }
  
  df_list = list()

  if (parallel) {
    
    if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
    if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
    
    cores = parallel::detectCores() - 1 
    message("Starting PARALLEL processing for ", scenario_name, " on ", cores, " cores...")
    cat(sprintf("Spawned %d parallel workers.\n", cores), file = log_file, append = TRUE)
    
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    df_list = foreach::foreach(
      s = 1:n_sims, 
      .packages = c("rstan", "dplyr", "here")
    ) %dopar% {
      
      # force background worker to source R functions
      fun_list = list.files(here::here("R"))
      for (i in 1:length(fun_list)) {
        source(here::here("R", fun_list[i]))
      }
      
      run_simulation(
        nyears            = nyears,
        init_nya          = init_nya,
        waa               = waa,
        nages             = nages,
        rec_regime_length = rec_regime_length,
        max_harvest_rate  = max_harvest_rate,
        survival_mean     = survival_mean,
        survival_sd       = survival_sd,
        sr_params         = sr_params,
        threshold         = thresholds,
        maturity          = maturity,
        selectivity       = selectivity,
        vb_params         = vb_params,
        bin_size          = bin_size,
        burn_in_length    = burn_in_length,
        hist_harvest_rate = hist_harvest_rate,
        estimation        = estimation,
        mcmc_setup        = mcmc_setup,
        est_method        = est_method,
        model_objects     = model_objects,
        rec_type          = rec_type,
        hcr_type          = hcr_type,
        sim_num           = s,
        scenario_name     = scenario_name,
        log_file          = log_file
      )
    }
    
    parallel::stopCluster(cl)
    message("Parallel Scenario ", scenario_name, " complete.")
    cat(sprintf("=== FINISHED PARALLEL BATCH: Scenario %s ===\n", scenario_name), file = log_file, append = TRUE)
    
  } else {

    message("Starting SEQUENTIAL processing for ", scenario_name, "...")
    cat("Running sequentially.\n", file = log_file, append = TRUE)
    
    for (s in 1:n_sims) {
      message("Running simulation ", s, " of ", n_sims)
      
      df_list[[s]] = run_simulation(
        nyears            = nyears,
        init_nya          = init_nya,
        waa               = waa,
        nages             = nages,
        rec_regime_length = rec_regime_length,
        max_harvest_rate  = max_harvest_rate,
        survival_mean     = survival_mean,
        survival_sd       = survival_sd,
        sr_params         = sr_params,
        threshold         = thresholds,
        maturity          = maturity,
        selectivity       = selectivity,
        vb_params         = vb_params,
        bin_size          = bin_size,
        burn_in_length    = burn_in_length,
        hist_harvest_rate = hist_harvest_rate,
        estimation        = estimation,
        mcmc_setup        = mcmc_setup,
        est_method        = est_method,
        model_objects     = model_objects,
        rec_type          = rec_type,
        hcr_type          = hcr_type,
        sim_num           = s,
        scenario_name     = scenario_name,
        log_file          = log_file
      )
    }
    
    message("Sequential Scenario ", scenario_name, " complete.")
    cat(sprintf("=== FINISHED SEQUENTIAL BATCH: Scenario %s ===\n", scenario_name), file = log_file, append = TRUE)
  }

  total_years = burn_in_length + nyears
  abs_biomass_mat = matrix(NA, nrow=n_sims, ncol=total_years)
  est_biomass_mat = matrix(NA, nrow=n_sims, ncol=nyears)
  catch_mat       = matrix(NA, nrow=n_sims, ncol=total_years)
  
  for(s in 1:n_sims) {
    abs_biomass_mat[s, ] = df_list[[s]]$biomass
    est_biomass_mat[s, ] = df_list[[s]]$est_biomass
    catch_mat[s, ]       = df_list[[s]]$catch
  }
  
  df_out = data.frame(
    year = 1:total_years,
    abs_biomass_mean = colMeans(abs_biomass_mat),
    abs_biomass_sd   = apply(abs_biomass_mat, 2, sd),
    catch_mean       = colMeans(catch_mat),
    catch_sd         = apply(catch_mat, 2, sd)
  )
  
  df_est = data.frame(
    year = (burn_in_length+1):total_years,
    est_biomass_mean = colMeans(est_biomass_mat),
    est_biomass_sd   = apply(est_biomass_mat, 2, sd)
  )
  
  df_out = merge(df_out, df_est, by="year", all.x=TRUE)
  
  if(plot) {
    true    = df_out$abs_biomass_mean / df_out$abs_biomass_mean[1]
    true_sd = df_out$abs_biomass_sd / df_out$abs_biomass_mean[1]
    
    par(mfrow = c(2, 1), mar = c(4,4,1,1))
    plot(1:nrow(df_out), true, type = "l", col = "blue", ylim = c(0, max(true + true_sd*2, na.rm=T)), ylab = "Relative Biomass", xlab = "Year")
    lines(1:nrow(df_out), true + 2*true_sd, col = "blue", lty = 2)
    lines(1:nrow(df_out), true - 2*true_sd, col = "blue", lty = 2)
    if(estimation) {
      lines((burn_in_length+1):total_years, df_out$est_biomass_mean, col="red")
    }
    abline(v = burn_in_length, lty = 2, col = "black")
    
    plot(1:nrow(df_out), df_out$catch_mean, type = "l", col = "darkgreen", ylab = "Catch", xlab = "Year", ylim = c(0, max(df_out$catch_mean + df_out$catch_sd*2, na.rm=T)))
    lines(1:nrow(df_out), df_out$catch_mean + 2*df_out$catch_sd, col = "darkgreen", lty = 2)
    lines(1:nrow(df_out), df_out$catch_mean - 2*df_out$catch_sd, col = "darkgreen", lty = 2)
    abline(v = burn_in_length, lty = 2, col = "black")
  }
  
  return(list(raw_sims = df_list, df = df_out))
}