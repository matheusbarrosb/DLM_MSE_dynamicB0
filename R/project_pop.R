project_pop = function(nages, waa, curr_nya, maturity,
                       recruitment, survival, threshold,
                       harvest_rate, recruit_type = NULL) {
  
  # ensure survival is between 0.1 and 0.9
  if (survival < 0.1) {
    survival = 0.1
  } else if (survival > 0.9) {
    survival = 0.9
  }
  
  next_nya    = numeric(nages)
  next_nya[1] = recruitment
  
  # ages 1 - (nages - 1)
  for (a in 2:(nages-1)) {
    next_nya[a] = curr_nya[a-1] * survival
  }
  
  # plus group
  next_nya[nages] = curr_nya[nages - 1] * survival + curr_nya[nages] * survival
  
  # ensure non negative
  next_nya[next_nya < 0] = 0
  
  next_nya_ghost = next_nya
  
  # calculate biomass
  spawn_biomass = sum(waa * maturity * next_nya)
  
  # define and apply catch, if biomass > threshold, harvest rate = h, else = 0
  catch_at_age = numeric(nages)
  if (spawn_biomass > threshold) {
    h = harvest_rate
    catch = spawn_biomass * h
    
    spawn_biomass_before = spawn_biomass
    
    # apply catch to population
    catch_at_age = (waa * maturity * next_nya) / spawn_biomass * catch
    next_nya = ifelse(waa > 0 & maturity > 0, 
                      next_nya - (catch_at_age / (waa * maturity)), 
                      next_nya)
    next_nya[next_nya < 0] = 0
    
    # update spawning biomass after catch
    spawn_biomass = sum(waa * maturity * next_nya)
    
    # calculate actual catch removed
    catch = spawn_biomass_before - spawn_biomass
    
  } else {
    h = 0
    catch = 0
  }
  
  return(list(next_nya       = next_nya,
              lengths        = lengths,
              spawn_biomass  = spawn_biomass,
              catch_at_age   = catch_at_age,
              total_catch    = catch))
  
}
