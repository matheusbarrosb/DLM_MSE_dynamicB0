control_rule = function(spawn_biomass, thresholds, max_harvest_rate) {
  
  lower = thresholds[1]
  upper = thresholds[2]
  
  h = ifelse(
    spawn_biomass < lower, 0,
    ifelse(
      spawn_biomass > upper, max_harvest_rate,
      max_harvest_rate * (spawn_biomass - lower) / (upper - lower)
    )
  )
  
  return(h)
}


