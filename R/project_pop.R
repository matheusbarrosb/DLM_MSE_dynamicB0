project_pop = function(nages, waa, selectivity, curr_nya, maturity,
                       recruitment, survival, thresholds, est_spawn_biomass = NULL,
                       max_harvest_rate, recruit_type = NULL, estimation = FALSE) {
  
  # convert instantaneous mortality to discrete survival
  M = -log(survival)
  
  # determine harvest rate h based on start of year SSB
  true_start_ssb = sum(curr_nya * waa * maturity)
  if (estimation == TRUE && !is.null(est_spawn_biomass)) {
    h = control_rule(spawn_biomass    = est_spawn_biomass,
                     thresholds       = thresholds,
                     max_harvest_rate = max_harvest_rate)
  } else {
    h = control_rule(spawn_biomass    = true_start_ssb,
                     thresholds       = thresholds,
                     max_harvest_rate = max_harvest_rate)
  }
  
  # calculate the biomass available to the fishery after applying natural mortality
  nya_post_M = curr_nya * exp(-M)
  vuln_bio_post_M = sum(nya_post_M * waa * selectivity)
  
  # realized harvest rate
  target_catch = true_start_ssb * h 
  
  U = target_catch / (vuln_bio_post_M + 1e-10)
  if (U > 0.99) U = 0.99
  
  # calculate survivors - Baranov Catch Equation
  # Survivors = N_start * exp(-M) * (1 - U * sel)
  survivors = curr_nya * exp(-M) * (1 - U * selectivity)
  
  # catch at age
  catch_at_age = curr_nya * exp(-M) * (U * selectivity) * waa

  next_nya    = numeric(nages)
  next_nya[1] = recruitment
  
  for (a in 1:(nages-1)) {
    next_nya[a+1] = survivors[a]
  }

  # add to plus group
  next_nya[nages] = next_nya[nages] + survivors[nages]
  
  total_yield = sum(catch_at_age)
  spawn_biomass_after = sum(next_nya * waa * maturity)
  
  return(list(next_nya       = next_nya,
              spawn_biomass  = spawn_biomass_after,
              catch_at_age   = catch_at_age,
              total_catch    = total_yield))
}