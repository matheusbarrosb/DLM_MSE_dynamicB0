# plot_hcr = function(thresholds, max_val, hcr_params = NULL, 
#                     slopes = c(0.0, -0.01, -0.025, -0.05, -0.15)) {
#   
#   if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
#   
#   if (is.null(hcr_params)) {
#     hcr_params = list(k_down = 20, k_up = 15, gamma = 5, trend_years = 3, anchor_fraction = 0.05)
#   }
#   
#   spr_seq = seq(0, 1, length.out = 150)
#   
#   # wrapper to calculate the multiplier for a given SPR and trajectory
#   calc_multiplier = function(spr, target_slope) {
#     # build a dummy history of length 'trend_years' to achieve the exact target slope
#     n_years = hcr_params$trend_years
#     if (n_years > 1) {
#       years = 1:n_years
#       # y = mx + b. 
#       # b = spr - target_slope * n_years
#       b = spr - target_slope * n_years
#       dummy_history = target_slope * years + b
#     } else {
#       dummy_history = c(spr)
#     }
#     
#     dummy_catch = 1000
#     
#     out = control_rule(indicator_history = dummy_history, 
#                        indicator_type    = "relative", 
#                        thresholds        = thresholds, 
#                        max_val           = max_val, 
#                        recent_catch      = dummy_catch, 
#                        anchor_catch      = 0, # Bypass the anchor for clean visualization
#                        hcr_params        = hcr_params)
#     
#     # Divide output by input to isolate the pure multiplier
#     return(out$value / dummy_catch)
#   }
#   
#   # Dynamically build the dataframe for any number of slopes
#   df_list = lapply(slopes, function(s) {
#     data.frame(
#       SPR        = spr_seq,
#       Slope      = s,
#       Multiplier = sapply(spr_seq, calc_multiplier, target_slope = s)
#     )
#   })
#   
#   # Combine the list of dataframes into one
#   df = do.call(rbind, df_list)
#   
#   # Plot
#   ggplot2::ggplot(df, ggplot2::aes(x = SPR, y = Multiplier, color = Slope, group = as.factor(Slope))) +
#     ggplot2::geom_line(size = 1.2) +
#     ggplot2::geom_vline(xintercept = thresholds[1], linetype = "dashed", color = "red") +
#     ggplot2::geom_vline(xintercept = thresholds[2], linetype = "dashed", color = "black") +
#     ggplot2::geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5) +
#     ggplot2::scale_color_continuous(type = "viridis") +
#     ggplot2::theme_minimal() +
#     ggplot2::labs(
#       subtitle = paste0("k_down: ", hcr_params$k_down, " | k_up: ", hcr_params$k_up, " | gamma: ", hcr_params$gamma),
#       x = "SPR", 
#       y = "Multiplier",
#       color = "Slope"
#     )
# }
plot_hcr = function(thresholds, max_val) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  
  spr_seq = seq(0, 1, length.out = 150)
  rules   = c("multiplier_hockey_stick", "proportional_multiplier")
  
  calc_multiplier = function(spr, r_type) {
    dummy_catch = 1000
    
    out = control_rule(indicator_history = c(spr), 
                       indicator_type    = "relative",
                       hcr_type          = r_type,
                       thresholds        = thresholds, 
                       max_val           = max_val, 
                       recent_catch      = dummy_catch, 
                       anchor_catch      = 0)
    
    return(out$value / dummy_catch)
  }
  
  df_list = list()
  for (r in rules) {
    df_list[[length(df_list) + 1]] = data.frame(
      SPR        = spr_seq,
      Rule       = r,
      Multiplier = sapply(spr_seq, calc_multiplier, r_type = r)
    )
  }
  
  df = do.call(rbind, df_list)
  
  ggplot2::ggplot(df, ggplot2::aes(x = SPR, y = Multiplier, color = Rule)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::geom_vline(xintercept = thresholds[1], linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = thresholds[2], linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5) +
    ggplot2::scale_color_manual(values = c("multiplier_hockey_stick" = "#1f77b4", 
                                           "proportional_multiplier" = "#ff7f0e")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Comparison of Relative HCR Shapes (Linear)",
      x = "Estimated Indicator (SPR or B/B0)", 
      y = "Catch Multiplier",
      color = "Rule Type"
    )
}