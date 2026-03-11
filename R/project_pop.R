project_pop = function(nages, waa, selectivity, curr_nya, maturity,
                       recruitment, survival, thresholds, est_spawn_biomass = NULL,
                       max_harvest_rate, recruit_type = NULL, estimation = FALSE,
                       recent_catch = NULL, anchor_catch = NULL,
                       hcr_type = "absolute_hockey_stick") { 
  
  M = -log(survival)
  true_start_ssb = sum(curr_nya * waa * maturity)
  
  nya_post_M      = curr_nya * exp(-M)
  vuln_bio_post_M = sum(nya_post_M * waa * selectivity)
  
  if (is.infinite(thresholds[1])) {
    target_catch = 0
  } else {
    
    if (estimation == TRUE && !is.null(est_spawn_biomass)) {
      
      is_absolute = grepl("absolute", hcr_type)
      
      rule_output = control_rule(indicator_history = est_spawn_biomass,
                                 indicator_type    = ifelse(is_absolute, "absolute", "relative"),
                                 hcr_type          = hcr_type,
                                 thresholds        = thresholds,
                                 max_val           = max_harvest_rate,
                                 recent_catch      = recent_catch,
                                 anchor_catch      = anchor_catch)
      
      if (rule_output$type == "harvest_rate") {
        target_catch = true_start_ssb * rule_output$value
      } else {
        target_catch = rule_output$value
      }
      
    } else {
      rule_output = control_rule(indicator_history = true_start_ssb,
                                 indicator_type    = "absolute",
                                 hcr_type          = "absolute_hockey_stick",
                                 thresholds        = thresholds,
                                 max_val           = max_harvest_rate)
      
      target_catch = true_start_ssb * rule_output$value 
    }
  }
  
  U = target_catch / (vuln_bio_post_M + 1e-10)
  if (U > 0.99) U = 0.99
  
  survivors             = curr_nya * exp(-M) * (1 - U * selectivity)
  catch_yield_at_age    = curr_nya * exp(-M) * (U * selectivity) * waa
  catch_numbers_at_age  = curr_nya * exp(-M) * (U * selectivity)
  sample_numbers_at_age = curr_nya * exp(-M) * selectivity
  
  next_nya    = numeric(nages)
  next_nya[1] = recruitment
  for (a in 1:(nages-1)) {
    next_nya[a+1] = survivors[a]
  }
  next_nya[nages] = next_nya[nages] + survivors[nages]
  
  total_yield         = sum(catch_yield_at_age)
  spawn_biomass_after = sum(next_nya * waa * maturity)
  
  return(list(next_nya       = next_nya,
              spawn_biomass  = spawn_biomass_after,
              catch_at_age   = catch_yield_at_age,   
              catch_numbers  = catch_numbers_at_age, 
              sample_numbers = sample_numbers_at_age, 
              total_catch    = total_yield))
}