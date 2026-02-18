# run_simulation = function(nyears, init_nya, waa, nages,
#                           rec_regime_length, max_harvest_rate,
#                           survival_mean, survival_sd, sr_params,
#                           threshold, maturity, selectivity, vb_params, bin_size,
#                           burn_in_length, hist_harvest_rate, estimation = FALSE,
#                           mcmc_setup, sca_model = NULL, est_method = "SCA", model_objects = NULL,
#                           rec_type = c("decoupled", "BV")) {
#   
#   total_years = burn_in_length + nyears
#   
#   nya_mat       = matrix(NA, nrow = total_years, ncol = nages)
#   nya_ghost_mat = matrix(NA, nrow = total_years, ncol = nages)
#   catch         = numeric(total_years)
#   seed_vec      = numeric(total_years)
#   est_biomass   = numeric(nyears)
#   
#   k    = rnorm(1, vb_params$k, vb_params$k_sd)
#   linf = rnorm(1, vb_params$linf, vb_params$linf_sd)
#   
#   trans_mat = make_transition_matrix(
#     k        = k,
#     linf     = linf,
#     sigma    = 1,
#     bin_size = bin_size
#   )
#   
#   mapping_mat = matrix(0, nrow = nages, ncol = ncol(trans_mat))
#   obs_bins = seq(0, linf, length.out = nages + 1)
#   for (i in 1:nages) {
#     lower = obs_bins[i]
#     upper = obs_bins[i+1]
#     indices = which(as.numeric(rownames(trans_mat)) >= lower &
#                       as.numeric(rownames(trans_mat)) < upper)
#     mapping_mat[i, indices] = 1
#   }
#   
#   collapsed = mapping_mat %*% trans_mat %*% t(mapping_mat)
#   prob_mat = sweep(collapsed, 1,  rowSums(collapsed), "/")
#   
#   length_mat = matrix(NA, nrow = total_years, ncol = ncol(prob_mat))
#   colnames(length_mat) = as.character(round(vb_params$linf * (1 - exp(-vb_params$k*((1:nages) - vb_params$t0)))))
#   
#   est_model_state   = NULL
#   est_spawn_biomass = NULL 
#   model_objects     = list(SCA = sca_model, SPM)
#   
#   for (y in 1:total_years) {
#     
#     seed_vec[y] = sample(1:1e6, 1)
#     set.seed(seed_vec[y])
#     
#     if (y == 1) {
#       curr_nya = init_nya
#       curr_ghost_nya = init_nya
#     } else {
#       curr_nya = real_projection$next_nya
#       curr_ghost_nya = ghost_projection$next_nya
#     }
#     
#     curr_ssb = sum(curr_nya * waa * maturity)
#     
#     curr_recruitment = compute_recruitment(
#       year              = y,
#       sr_params         = sr_params,
#       rec_regime_length = rec_regime_length,
#       type              = rec_type,
#       ssb               = curr_ssb
#     )
#     
#     curr_survival = rnorm(1, mean = survival_mean, sd = survival_sd)
#     
#     if (y <= burn_in_length) {
#       current_thresholds = c(0, 0)
#       current_max_hr     = hist_harvest_rate
#       est_spawn_biomass  = NULL    
#       run_estimation     = FALSE
#       hcr_biomass_input  = NULL
#       
#     } else {
#       current_thresholds = threshold
#       current_max_hr     = max_harvest_rate
#       run_estimation     = estimation 
#       
#       if (run_estimation == TRUE) {
#         if(!is.null(est_spawn_biomass)) {
#           hcr_biomass_input = est_spawn_biomass
#         } else {
#           hcr_biomass_input = curr_ssb/sr_params$SSB0
#         }
#       } else {
#         hcr_biomass_input = curr_ssb/sr_params$SSB0
#       }
#     }
#     
#     real_projection = project_pop(
#       threshold         = current_thresholds, 
#       nages             = nages,
#       waa               = waa,
#       selectivity       = selectivity,
#       curr_nya          = curr_nya,
#       maturity          = maturity,
#       recruitment       = curr_recruitment,
#       max_harvest_rate  = current_max_hr,   
#       survival          = curr_survival,
#       estimation        = run_estimation,
#       est_spawn_biomass = if(y > burn_in_length) hcr_biomass_input else NULL
#     )
#     
#     ghost_projection = project_pop(
#       threshold        = Inf, 
#       nages            = nages,
#       waa              = waa,
#       selectivity      = selectivity,
#       curr_nya         = curr_ghost_nya,
#       maturity         = maturity,
#       recruitment      = curr_recruitment,
#       max_harvest_rate = max_harvest_rate,
#       survival         = curr_survival,
#       estimation       = FALSE, 
#       est_spawn_biomass = NULL
#     )
#     
#     nya_mat[y, ]      = curr_nya
#     nya_ghost_mat[y,] = curr_ghost_nya
#     catch[y]          = real_projection$total_catch
#     current_lengths   = real_projection$next_nya %*% prob_mat
#     length_mat[y, ]   = current_lengths
#     
#     # estimation initialization
#     if (y == burn_in_length && estimation == TRUE) {
#       
#       history_data = list(
#         sim            = list(catch = catch, nya_mat = nya_mat),
#         waa            = waa,
#         maturity       = maturity,
#         selectivity    = selectivity,
#         nages          = nages,
#         burn_in_length = burn_in_length,
#         sigma_index    = 0.05,
#         sigma_catch    = 0,
#         ess            = 250
#       )
#       
#       assessment_result = estimate_stock_status(
#         method        = est_method,
#         model_state   = NULL,
#         history_data  = history_data,
#         mcmc_setup    = mcmc_setup,
#         model_objects = model_objects
#       )
#       
#       est_model_state   = assessment_result$updated_state
#       est_spawn_biomass = assessment_result$status
#       
#       message("Initialized Method: ", est_method, " | Initial B/B0: ", round(est_spawn_biomass, 2))
#     }
#     
#     # estimation projection
#     if (y > burn_in_length && estimation == TRUE) {
#       
#       current_obs = list(
#         catch       = catch[y],
#         nya         = nya_mat[y, ],
#         waa         = waa,
#         selectivity = selectivity,
#         sigma_index = 0.05,
#         sigma_catch = 0,
#         ess         = 250
#       )
#       
#       assessment_result = estimate_stock_status(
#         method        = est_method,
#         model_state   = est_model_state,
#         current_obs   = current_obs,
#         mcmc_setup    = mcmc_setup,
#         model_objects = model_objects
#       )
#       
#       est_model_state                 = assessment_result$updated_state
#       est_spawn_biomass               = assessment_result$status
#       est_biomass[y - burn_in_length] = est_spawn_biomass
#       
#       message(
#         "Method: ", est_method, 
#         " | Catch =", round(catch[y], 0), 
#         " | True B/B0 =", round(curr_ssb/sr_params$SSB0, 3), 
#         " | Est B/B0 =", round(est_spawn_biomass, 3)
#       )
#     }
#     
#   } 
#   
#   biomass       = rowSums(nya_mat * matrix(waa * maturity, nrow = total_years, ncol = nages, byrow = TRUE))
#   ghost_biomass = rowSums(nya_ghost_mat * matrix(waa * maturity, nrow = total_years, ncol = nages, byrow = TRUE))
#   
#   biomass[biomass < 0] = 0
#   nya_mat[nya_mat < 0] = 0
#   length_mat[length_mat < 0] = 0
#   
#   dimnames(nya_mat) = list(paste0("Year_", 1:total_years), paste0("Age_", 1:nages))
#   dimnames(nya_ghost_mat) = dimnames(nya_mat)
#   
#   return(
#     list(
#       nya_mat       = nya_mat,
#       ghost_nya_mat = nya_ghost_mat,
#       length_mat    = length_mat,
#       biomass       = biomass,
#       est_biomass   = est_biomass,
#       ghost_biomass = ghost_biomass,
#       catch         = catch,
#       seeds         = seed_vec,
#       burn_in_end   = burn_in_length 
#     )
#   )
# }

run_simulation = function(nyears, init_nya, waa, nages,
                          rec_regime_length, max_harvest_rate,
                          survival_mean, survival_sd, sr_params,
                          threshold, maturity, selectivity, vb_params, bin_size,
                          burn_in_length, hist_harvest_rate, estimation = FALSE,
                          mcmc_setup, est_method = "SCA", model_objects = NULL,
                          rec_type = c("decoupled", "BV")) {
  
  total_years = burn_in_length + nyears
  
  nya_mat       = matrix(NA, nrow = total_years, ncol = nages)
  nya_ghost_mat = matrix(NA, nrow = total_years, ncol = nages)
  catch         = numeric(total_years)
  seed_vec      = numeric(total_years)
  est_biomass   = numeric(nyears)
  
  k    = rnorm(1, vb_params$k, vb_params$k_sd)
  linf = rnorm(1, vb_params$linf, vb_params$linf_sd)
  
  trans_mat = make_transition_matrix(
    k        = k,
    linf     = linf,
    sigma    = 1,
    bin_size = bin_size
  )
  
  mapping_mat = matrix(0, nrow = nages, ncol = ncol(trans_mat))
  obs_bins = seq(0, linf, length.out = nages + 1)
  for (i in 1:nages) {
    lower = obs_bins[i]
    upper = obs_bins[i+1]
    indices = which(as.numeric(rownames(trans_mat)) >= lower &
                      as.numeric(rownames(trans_mat)) < upper)
    mapping_mat[i, indices] = 1
  }
  
  collapsed = mapping_mat %*% trans_mat %*% t(mapping_mat)
  prob_mat = sweep(collapsed, 1,  rowSums(collapsed), "/")
  
  length_mat = matrix(NA, nrow = total_years, ncol = ncol(prob_mat))
  colnames(length_mat) = as.character(round(vb_params$linf * (1 - exp(-vb_params$k*((1:nages) - vb_params$t0)))))
  
  est_model_state   = NULL
  est_spawn_biomass = NULL 
  
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
      
      if (run_estimation == TRUE) {
        if(!is.null(est_spawn_biomass)) {
          hcr_biomass_input = est_spawn_biomass
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
      est_spawn_biomass = if(y > burn_in_length) hcr_biomass_input else NULL
    )
    
    ghost_projection = project_pop(
      threshold        = Inf, 
      nages            = nages,
      waa              = waa,
      selectivity      = selectivity,
      curr_nya         = curr_ghost_nya,
      maturity         = maturity,
      recruitment      = curr_recruitment,
      max_harvest_rate = max_harvest_rate,
      survival         = curr_survival,
      estimation       = FALSE, 
      est_spawn_biomass = NULL
    )
    
    nya_mat[y, ]      = curr_nya
    nya_ghost_mat[y,] = curr_ghost_nya
    catch[y]          = real_projection$total_catch
    current_lengths   = real_projection$next_nya %*% prob_mat
    length_mat[y, ]   = current_lengths
    
    if (y == burn_in_length && estimation == TRUE) {
      
      history_data = list(
        sim            = list(catch = catch, nya_mat = nya_mat),
        waa            = waa,
        maturity       = maturity,
        selectivity    = selectivity,
        nages          = nages,
        burn_in_length = burn_in_length,
        sigma_index    = 0.05,
        sigma_catch    = 0,
        ess            = 250
      )
      
      assessment_result = estimate_stock_status(
        method        = est_method,
        model_state   = NULL,
        history_data  = history_data,
        mcmc_setup    = mcmc_setup,
        model_objects = model_objects
      )
      
      est_model_state   = assessment_result$updated_state
      est_spawn_biomass = assessment_result$status
      
      message("Initialized Method: ", est_method, " | Initial B/B0: ", round(est_spawn_biomass, 2))
    }
    
    if (y > burn_in_length && estimation == TRUE) {
      
      current_obs = list(
        catch       = catch[y],
        nya         = nya_mat[y, ],
        waa         = waa,
        selectivity = selectivity,
        sigma_index = 0.05,
        sigma_catch = 0,
        ess         = 250
      )
      
      assessment_result = estimate_stock_status(
        method        = est_method,
        model_state   = est_model_state,
        current_obs   = current_obs,
        mcmc_setup    = mcmc_setup,
        model_objects = model_objects
      )
      
      est_model_state                 = assessment_result$updated_state
      est_spawn_biomass               = assessment_result$status
      est_biomass[y - burn_in_length] = est_spawn_biomass
      
      message(
        "Method: ", est_method, 
        " | Catch =", round(catch[y], 0), 
        " | True B/B0 =", round(curr_ssb/sr_params$SSB0, 3), 
        " | Est B/B0 =", round(est_spawn_biomass, 3)
      )
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
      seeds         = seed_vec,
      burn_in_end   = burn_in_length 
    )
  )
}