# control_rule = function(indicator_history, indicator_type = c("absolute", "relative"), 
#                         thresholds, max_val, recent_catch = NULL, anchor_catch = NULL,
#                         hcr_params = NULL) {
#   
#   indicator_type = match.arg(indicator_type)
#   
#   # default parameters if none are provided
#   if (is.null(hcr_params)) {
#     hcr_params = list(
#       k_down          = 20,
#       k_up            = 15,
#       gamma           = 5,
#       trend_years     = 3,
#       anchor_fraction = 0.05
#     )
#   }
#   
#   # the most recent estimate is the current status
#   current_val = tail(indicator_history, 1)
#   limit       = thresholds[1]
#   target      = thresholds[2]
#   
#   # absolute biomass hcr (sca, spm) ---------------------
#   if (indicator_type == "absolute") {
#     
#     h = ifelse(
#       current_val < limit, 0,
#       ifelse(
#         current_val > target, max_val, 
#         max_val * (current_val - limit) / (target - limit)
#       )
#     )
#     
#     return(list(type = "harvest_rate", value = h))
#     
#   } else if (indicator_type == "relative") {
#     # relative empirical hcr (lbspr) ----------------------
#     
#     if (is.null(recent_catch)) stop("recent_catch required for relative HCR")
#     
#     # base asymmetric logistic multiplier
#     if (current_val <= limit) {
#       base_mult = 0
#     } else if (current_val <= target) {
#       k_down = hcr_params$k_down
#       x_mid  = target - ((target - limit) * 0.25) 
#       S_x    = 1 / (1 + exp(-k_down * (current_val - x_mid)))
#       S_L    = 1 / (1 + exp(-k_down * (limit - x_mid)))
#       S_T    = 1 / (1 + exp(-k_down * (target - x_mid)))
#       base_mult = (S_x - S_L) / (S_T - S_L) 
#     } else {
#       k_up  = hcr_params$k_up
#       x_mid = target + 0.1
#       S_x   = 1 / (1 + exp(-k_up * (current_val - x_mid)))
#       S_T   = 1 / (1 + exp(-k_up * (target - x_mid)))
#       S_inf = 1.0
#       base_mult = 1 + (max_val - 1) * ((S_x - S_T) / (S_inf - S_T))
#     }
#     
#     # trend modifier (slope over recent history)
#     trend_mult = 1.0
#     
#     hist_subset = tail(indicator_history, hcr_params$trend_years)
#     n_years     = length(hist_subset)
#     
#     if (n_years >= hcr_params$trend_years && hcr_params$trend_years > 1) {
#       years = 1:n_years
#       slope = coef(lm(hist_subset ~ years))[2]
#       if (slope < 0) {
#         trend_mult = max(1, 1 + (hcr_params$gamma * slope))
#       }
#     }
#     
#     if (!is.null(anchor_catch) && anchor_catch > 0) {
#       if (recent_catch < hcr_params$anchor_fraction * anchor_catch) {
#         recent_catch = hcr_params$anchor_fraction * anchor_catch
#       }
#     } else {
#       if (recent_catch < 10) recent_catch = 10 
#     }
#     
#     target_catch = recent_catch * base_mult * trend_mult
#     
#     return(list(type = "catch", value = target_catch))
#   }
# }
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