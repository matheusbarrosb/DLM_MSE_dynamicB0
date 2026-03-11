run_simulation = function(nyears, init_nya, waa, nages,
                          rec_regime_length, max_harvest_rate,
                          survival_mean, survival_sd, sr_params,
                          threshold, maturity, selectivity, vb_params, bin_size,
                          burn_in_length, hist_harvest_rate, estimation = FALSE,
                          mcmc_setup, est_method = "SCA", model_objects = NULL,
                          hcr_type = "absolute_hockey_stick",
                          rec_type = c("decoupled", "BV"),
                          sim_num = 1, scenario_name = "Default", log_file = "mse_log.txt") {
  
  total_years = burn_in_length + nyears
  
  nya_mat       = matrix(NA, nrow = total_years, ncol = nages)
  nya_ghost_mat = matrix(NA, nrow = total_years, ncol = nages)
  catch         = numeric(total_years)
  seed_vec      = numeric(total_years)
  est_biomass   = numeric(nyears)
  
  k    = rnorm(1, vb_params$k, vb_params$k_sd)
  linf = rnorm(1, vb_params$linf, vb_params$linf_sd)
  
  max_len  = linf * 1.5 
  len_bins = seq(0, max_len, by = bin_size)
  len_mids = len_bins[-length(len_bins)] + bin_size/2
  n_bins   = length(len_mids)
  
  prob_mat = matrix(0, nrow = nages, ncol = n_bins)
  colnames(prob_mat) = round(len_mids, 2)
  
  for(a in 1:nages) {
    mean_l  = linf * (1 - exp(-k * (a - vb_params$t0)))
    sigma_l = mean_l * 0.1 
    probs = pnorm(len_bins[-1], mean_l, sigma_l) - pnorm(len_bins[-(n_bins+1)], mean_l, sigma_l)
    prob_mat[a, ] = probs / sum(probs)
  }
  
  length_mat = matrix(NA, nrow = total_years, ncol = n_bins)
  colnames(length_mat) = colnames(prob_mat)
  
  est_model_state   = NULL
  est_spawn_biomass = NULL
  est_history       = c()
  
  for (y in 1:total_years) {
    
    seed_vec[y] = sample(1:1e6, 1)
    set.seed(seed_vec[y])
    
    if (y == 1) {
      curr_nya = init_nya
      curr_ghost_nya = init_nya
    } else {
      curr_nya = real_projection$next_nya
      curr_ghost_nya = ghost_projection$next_nya
    }
    
    curr_ssb = sum(curr_nya * waa * maturity)
    
    curr_recruitment = compute_recruitment(
      year              = y,
      sr_params         = sr_params,
      rec_regime_length = rec_regime_length,
      type              = rec_type,
      ssb               = curr_ssb
    )
    
    curr_survival = rnorm(1, mean = survival_mean, sd = survival_sd)
    
    if (y <= burn_in_length) {
      current_thresholds = c(0, 0)
      current_max_hr     = hist_harvest_rate
      est_spawn_biomass  = NULL    
      run_estimation     = FALSE
      hcr_biomass_input  = NULL
      
    } else {
      current_thresholds = threshold
      current_max_hr     = max_harvest_rate
      run_estimation     = estimation 
      
      if (run_estimation == TRUE && y > burn_in_length) {
        if(length(est_history) > 0) {
          hcr_biomass_input = tail(est_history, 1) 
        } else {
          hcr_biomass_input = curr_ssb/sr_params$SSB0
        }
      } else {
        hcr_biomass_input = curr_ssb/sr_params$SSB0
      }
    }
    
    real_projection = project_pop(
      threshold         = current_thresholds, 
      nages             = nages,
      waa               = waa,
      selectivity       = selectivity,
      curr_nya          = curr_nya,
      maturity          = maturity,
      recruitment       = curr_recruitment,
      max_harvest_rate  = current_max_hr,   
      survival          = curr_survival,
      estimation        = run_estimation,
      est_spawn_biomass = if(y > burn_in_length) hcr_biomass_input else NULL,
      recent_catch      = if(y > burn_in_length) catch[y-1] else 0,
      anchor_catch      = if(y > burn_in_length) mean(catch[1:burn_in_length]) else 0,
      hcr_type          = hcr_type
    )
    
    ghost_projection = project_pop(
      threshold         = Inf, 
      nages             = nages,
      waa               = waa,
      selectivity       = selectivity,
      curr_nya          = curr_ghost_nya,
      maturity          =  maturity,
      recruitment       = curr_recruitment,
      max_harvest_rate  = max_harvest_rate,
      survival          = curr_survival,
      estimation        = FALSE, 
      est_spawn_biomass = NULL,
      recent_catch      = 0,
      anchor_catch      = 0,
      hcr_type          = hcr_type
    )
    
    nya_mat[y, ]      = curr_nya
    nya_ghost_mat[y,] = curr_ghost_nya
    catch[y]          = real_projection$total_catch
    
    current_lengths   = real_projection$sample_numbers %*% prob_mat
    length_mat[y, ]   = current_lengths
    
    # -------------------------------------------------------------------------
    # BURN-IN ESTIMATION INITIALIZATION
    # -------------------------------------------------------------------------
    if (y == burn_in_length && estimation == TRUE) {
      
      history_data = list(
        sim            = list(catch = catch, nya_mat = nya_mat, length_mat = length_mat),
        waa            = waa,
        maturity       = maturity,
        selectivity    = selectivity,
        nages          = nages,
        burn_in_length = burn_in_length,
        sigma_index    = 0.05,
        sigma_catch    = 0,
        ess            = 250,
        vb_params      = vb_params,
        len_bins       = len_bins,
        M              = -log(survival_mean)
      )
      
      # Safe Execution Wrapper
      assessment_result = tryCatch({
        estimate_stock_status(
          method        = est_method,
          model_state   = NULL,
          history_data  = history_data,
          mcmc_setup    = mcmc_setup,
          model_objects = model_objects
        )
      }, error = function(e) {
        err_msg = sprintf("[Scenario: %s | Sim: %03d | Year: %02d] ERROR IN ASSESSMENT: %s", scenario_name, sim_num, y, e$message)
        cat(err_msg, "\n", file = log_file, append = TRUE)
        return(NULL)
      })
      
      if (!is.null(assessment_result)) {
        est_model_state   = assessment_result$updated_state
        est_spawn_biomass = assessment_result$status
      } else {
        # Fallback if initialization fails
        est_model_state   = NULL
        est_spawn_biomass = 1.0 
      }
      
      est_history = c(est_history, est_spawn_biomass)
      
      log_msg = sprintf("[Scenario: %s | Sim: %03d | Year: %02d] Initialized Method: %s | Initial B/B0: %.2f", 
                        scenario_name, sim_num, y, est_method, est_spawn_biomass)
      cat(log_msg, "\n", file = log_file, append = TRUE)
    }
    
    # -------------------------------------------------------------------------
    # PROJECTION ESTIMATION LOOP
    # -------------------------------------------------------------------------
    if (y > burn_in_length && estimation == TRUE) {
      
      current_obs = list(
        catch       = catch[y],
        nya         = nya_mat[y, ],
        len_comp    = length_mat[y, ],
        waa         = waa,
        selectivity = selectivity,
        sigma_index = 0.05,
        sigma_catch = 0,
        ess         = 250
      )
      
      # Safe Execution Wrapper
      assessment_result = tryCatch({
        estimate_stock_status(
          method        = est_method,
          model_state   = est_model_state,
          current_obs   = current_obs,
          mcmc_setup    = mcmc_setup,
          model_objects = model_objects
        )
      }, error = function(e) {
        err_msg = sprintf("[Scenario: %s | Sim: %03d | Year: %02d] ERROR IN ASSESSMENT: %s", scenario_name, sim_num, y, e$message)
        cat(err_msg, "\n", file = log_file, append = TRUE)
        return(NULL)
      })
      
      if (!is.null(assessment_result)) {
        est_model_state   = assessment_result$updated_state
        est_spawn_biomass = assessment_result$status
      } else {
        # Fallback to the previous year's estimate to keep the loop alive
        est_spawn_biomass = tail(est_history, 1)
      }
      
      est_biomass[y - burn_in_length] = est_spawn_biomass
      est_history                     = c(est_history, est_spawn_biomass)
      
      log_msg = sprintf("[Scenario: %s | Sim: %03d | Year: %02d] Method: %s | Catch = %d | True B/B0 = %.3f | Est B/B0 = %.3f", 
                        scenario_name, sim_num, y, est_method, round(catch[y]), curr_ssb/sr_params$SSB0, est_spawn_biomass)
      cat(log_msg, "\n", file = log_file, append = TRUE)
    }
  } 
  
  biomass       = rowSums(nya_mat * matrix(waa * maturity, nrow = total_years, ncol = nages, byrow = TRUE))
  ghost_biomass = rowSums(nya_ghost_mat * matrix(waa * maturity, nrow = total_years, ncol = nages, byrow = TRUE))
  
  biomass[biomass < 0] = 0
  nya_mat[nya_mat < 0] = 0
  length_mat[length_mat < 0] = 0
  
  dimnames(nya_mat) = list(paste0("Year_", 1:total_years), paste0("Age_", 1:nages))
  dimnames(nya_ghost_mat) = dimnames(nya_mat)
  
  return(
    list(
      nya_mat       = nya_mat,
      ghost_nya_mat = nya_ghost_mat,
      length_mat    = length_mat,
      biomass       = biomass,
      est_biomass   = est_biomass,
      ghost_biomass = ghost_biomass,
      catch         = catch,
      em_fit        = if(estimation && !is.null(assessment_result)) assessment_result else NULL,
      seeds         = seed_vec,
      burn_in_end   = burn_in_length 
    )
  )
}