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
  total_time = total_years
  
  # make containers
  abs_biomass_mat  = matrix(NA, nrow = n_sims, ncol = total_years)
  rel_biomass_mat  = matrix(NA, nrow = n_sims, ncol = total_years)
  est_biomass_mat  = matrix(NA, nrow = n_sims, ncol = nyears)
  rec_mat          = matrix(NA, nrow = n_sims, ncol = total_years)
  catch_mat        = matrix(NA, nrow = n_sims, ncol = total_years)
  mean_lengths_mat = matrix(NA, nrow = n_sims, ncol = total_years)
  
  for(s in 1:n_sims) {
    abs_biomass_mat[s, ]  = df_list[[s]]$biomass
    rel_biomass_mat[s, ]  = df_list[[s]]$biomass / df_list[[s]]$ghost_biomass
    est_biomass_mat[s, ]  = df_list[[s]]$est_biomass
    catch_mat[s, ]        = df_list[[s]]$catch
    rec_mat[s, ]          = df_list[[s]]$nya_mat[,1]
    
    # Calculate mean lengths per year
    temp_lengths = numeric(total_years)
    for(y in 1:total_years) {
      temp_lengths[y] = sum(df_list[[s]]$length_mat[y, ] * as.numeric(colnames(df_list[[s]]$length_mat))) / sum(df_list[[s]]$length_mat[y, ])
    }
    mean_lengths_mat[s, ] = temp_lengths
  }
  
  # summarize outputs
  abs_biomass_mean = colMeans(abs_biomass_mat, na.rm = TRUE)
  abs_biomass_sd   = apply(abs_biomass_mat, 2, sd, na.rm = TRUE)
  rel_biomass_mean = colMeans(rel_biomass_mat, na.rm = TRUE)
  rel_biomass_sd   = apply(rel_biomass_mat, 2, sd, na.rm = TRUE)
  
  est_biomass_mean = colMeans(est_biomass_mat, na.rm = TRUE)
  est_biomass_mean = c(rep(NA, burn_in_length), est_biomass_mean)
  est_biomass_sd   = apply(est_biomass_mat, 2, sd, na.rm = TRUE)
  est_biomass_sd   = c(rep(NA, burn_in_length), est_biomass_sd)
  
  catch_mean       = colMeans(catch_mat, na.rm = TRUE)
  catch_sd         = apply(catch_mat, 2, sd, na.rm = TRUE)
  recruitment_mean = colMeans(rec_mat, na.rm = TRUE)
  recruitment_sd   = apply(rec_mat, 2, sd, na.rm = TRUE)
  mean_lenghts_mu  = colMeans(mean_lengths_mat, na.rm = TRUE)
  mean_lengths_sd  = apply(mean_lengths_mat, 2, sd, na.rm = TRUE)
  
  # plotting
  if (plot == TRUE) {
    
    par(mfrow = c(2,3), 
        mar = c(3,3,2,1), 
        oma = c(1,1,1,1), 
        mgp = c(2,0.7,0), 
        tcl = -0.3
    )
    
    plot(abs_biomass_mean, type = "l", xlab = "Year", ylab = "Spawning biomass", ylim = c(0, max( (abs_biomass_mean + abs_biomass_sd)*1.2, na.rm=TRUE)))
    lines(abs_biomass_mean + abs_biomass_sd*1.96, lty = 2)
    lines(abs_biomass_mean - abs_biomass_sd*1.96, lty = 2)
    for(i in 1:min(10, n_sims)) {
      lines(abs_biomass_mat[i,], col = rgb(0,0,0,0.2))
    }
    abline(v = burn_in_length, lty = 2, col = "blue")
    text(x = burn_in_length/2, y = max( (abs_biomass_mean + abs_biomass_sd)*1.2, na.rm=TRUE)*0.8, labels = "Burn-in", col = "blue")
    text(x = burn_in_length + (total_time - burn_in_length)/2, y = max( (abs_biomass_mean + abs_biomass_sd)*1.2, na.rm=TRUE)*0.8, labels = "Projection", col = "blue")
    box(lwd = 2)
    
    plot(rel_biomass_mean, type = "l", xlab = "Year", ylab = expression(B/B[unfished]), ylim = c(0, max( (rel_biomass_mean + rel_biomass_sd)*1.2, na.rm=TRUE)))
    lines(rel_biomass_mean + rel_biomass_sd*1.96, lty = 2)
    lines(rel_biomass_mean - rel_biomass_sd*1.96, lty = 2)
    for(i in 1:min(10, n_sims)) {
      lines(rel_biomass_mat[i,], col = rgb(0,0,0,0.2))
    }
    abline(h = 0.5, lty = 2, col = "red")
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)
    
    plot(recruitment_mean, type = "l", xlab = "Year", ylab = "Recruitment (N)", ylim = c(0, max( (recruitment_mean + recruitment_sd)*1.2, na.rm=TRUE)))
    lines(recruitment_mean + recruitment_sd*1.96, lty = 2)
    lines(recruitment_mean - recruitment_sd*1.96, lty = 2)
    for(i in 1:min(5, n_sims)) {
      lines(rec_mat[i,], col = rgb(0,0,0,0.2))
    }
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)
    
    plot(catch_mean, type = "p", pch = 20, cex = 0.8,
         xlab = "Year", ylab = "Catch (t)", ylim = c(0, max( (catch_mean + catch_sd)*1.2, na.rm=TRUE)))
    lines(catch_mean + catch_sd*1.96, lty = 2)
    lines(catch_mean - catch_sd*1.96, lty = 2)
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)
    
    plot(mean_lenghts_mu, type = "l", pch = 20, cex = 1,
         xlab = "Age", ylab = "Mean length (cm)", ylim = c(min((mean_lenghts_mu - mean_lengths_sd)*0.8, na.rm=TRUE),
                                                           max((mean_lenghts_mu + mean_lengths_sd)*1.2, na.rm=TRUE)))
    lines(mean_lenghts_mu + mean_lengths_sd*1.96, lty = 2)
    lines(mean_lenghts_mu - mean_lengths_sd*1.96, lty = 2)
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)
    
  }
  
  df_out = data.frame(
    year             = 1:total_years,
    abs_biomass_mean = abs_biomass_mean,
    abs_biomass_sd   = abs_biomass_sd,
    rel_biomass_mean = rel_biomass_mean,
    rel_biomass_sd   = rel_biomass_sd,
    est_biomass_mean = est_biomass_mean,
    est_biomass_sd   = est_biomass_sd,
    catch_mean       = catch_mean,
    catch_sd         = catch_sd,
    mean_lenghts_mu  = mean_lenghts_mu,
    mean_lengths_sd  = mean_lengths_sd
  )
  
  output = list(
    raw_sims = df_list,
    df = df_out
  )
  
  return(output)
}