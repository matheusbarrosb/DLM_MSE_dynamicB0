control_rule = function(indicator_history, indicator_type = c("absolute", "relative"), 
                        hcr_type = c("absolute_hockey_stick", "multiplier_hockey_stick", "proportional_multiplier"),
                        thresholds, max_val, recent_catch = NULL, anchor_catch = NULL) {
  
  indicator_type = match.arg(indicator_type)
  hcr_type       = match.arg(hcr_type)
  
  current_val = tail(indicator_history, 1)
  limit       = thresholds[1]
  target      = thresholds[2]
  
  if (indicator_type == "absolute") {
    
    if (hcr_type != "absolute_hockey_stick") hcr_type = "absolute_hockey_stick"
    
    if (hcr_type == "absolute_hockey_stick") {
      h = ifelse(
        current_val < limit, 0,
        ifelse(
          current_val > target, max_val, 
          max_val * (current_val - limit) / (target - limit)
        )
      )
      return(list(type = "harvest_rate", value = h))
    }
  } 
  
  else if (indicator_type == "relative") {
    
    if (is.null(recent_catch)) stop("recent_catch required for relative HCR")
    
    base_mult = 1.0 
    
    # 1. Multiplier hockey stick 
    if (hcr_type == "multiplier_hockey_stick") {
      base_mult = ifelse(
        current_val < limit, 0,
        ifelse(
          current_val > target, max_val, 
          max_val * (current_val - limit) / (target - limit)
        )
      )
    }
    
    # 2. Proportional multiplier 
    else if (hcr_type == "proportional_multiplier") {
      if (current_val <= limit) {
        base_mult = 0
      } else {
        base_mult = current_val / target
        base_mult = min(base_mult, max_val) 
      }
    }
    
    if (!is.null(anchor_catch) && anchor_catch > 0) {
      if (recent_catch < 0.05 * anchor_catch) {
        recent_catch = 0.05 * anchor_catch
      }
    } else {
      if (recent_catch < 10) recent_catch = 10 
    }
    
    # Final quota
    target_catch = recent_catch * base_mult
    return(list(type = "catch", value = target_catch))
  }
}