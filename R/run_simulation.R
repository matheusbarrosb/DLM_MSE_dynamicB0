run_simulation = function(nyears, init_nya, waa, nages,
                          rec_regime_length, harvest_rate,
                          survival_mean, survival_sd, sr_params,
                          threshold, maturity, vb_params, bin_size,
                          rec_type = c("decoupled", "BV")) {
  
  nya_mat       = matrix(NA, nrow = nyears, ncol = nages)
  nya_ghost_mat = matrix(NA, nrow = nyears, ncol = nages)
  catch         = numeric(nyears)
  seed_vec      = numeric(nyears)
  
  for (y in 1:nyears) {
    
    # random seed
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

    # alternate recruitment between high and low productivity regimes every rec_regime_length years
    # if (((y - 1) %/% rec_regime_length) %% 2 == 0) {
    #   curr_recruiment = rnorm(1, mean = recruitment_mean * 1.5, sd = recruitment_sd)
    # } else {
    #   curr_recruiment = rnorm(1, mean = recruitment_mean * 0.25, sd = recruitment_sd)
    # }
    
    curr_recruitment = compute_recruitment(
      
        year              = y,
        sr_params         = sr_params,
        rec_regime_length = rec_regime_length,
        type              = rec_type,
        ssb               = curr_ssb
      )
    
    # generate survival
    curr_survival = rnorm(1, mean = survival_mean, sd = survival_sd)
    
    # realized projection
    real_projection = project_pop(
      
      threshold    = threshold,
      nages        = nages,
      waa          = waa,
      curr_nya     = curr_nya,
      maturity     = maturity,
      recruitment  = curr_recruitment,
      harvest_rate = harvest_rate,
      survival     = curr_survival
    )
    
    # ghost projection to track dynamic B0
    ghost_projection = project_pop(
      
      threshold    = Inf, # to avoid catch
      nages        = nages,
      waa          = waa,
      curr_nya     = curr_ghost_nya,
      maturity     = maturity,
      recruitment  = curr_recruitment,
      harvest_rate = 0,
      survival     = curr_survival
    )
    
    nya_mat[y, ]      = real_projection$next_nya
    nya_ghost_mat[y,] = ghost_projection$next_nya
    catch[y]          = real_projection$total_catch  
    
  }
  
  biomass       = rowSums(nya_mat * matrix(waa * maturity, nrow = nyears, ncol = nages, byrow = TRUE))
  ghost_biomass = rowSums(nya_ghost_mat * matrix(waa * maturity, nrow = nyears, ncol = nages, byrow = TRUE))
  
  # calculate lengths from transition matrix
  k    = rnorm(1, vb_params$k, vb_params$k_sd)
  linf = rnorm(1, vb_params$linf, vb_params$linf_sd)
  
  trans_mat = 
    make_transition_matrix(
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
  
  length_mat = nya_mat %*% prob_mat
  colnames(length_mat) = as.character(round(vb_params$linf * (1 - exp(-vb_params$k*((1:nages) - vb_params$t0)))))
  
  # ensure non-negative biomass, nya, and lengths
  biomass[biomass < 0] = 0
  nya_mat[nya_mat < 0] = 0
  length_mat[length_mat < 0] = 0

  dimnames(nya_mat) = list(
    paste0("Year_", 1:nyears),
    paste0("Age_", 1:nages)
  );dimnames(nya_ghost_mat) = dimnames(nya_mat)
  
  return(
    list(
      nya_mat       = nya_mat,
      ghost_nya_mat = nya_ghost_mat,
      length_mat    = length_mat,
      biomass       = biomass,
      ghost_biomass = ghost_biomass,
      catch         = catch,
      seeds         = seed_vec
    )
  )
  
}
