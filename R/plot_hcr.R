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
    custom_theme() +
    ggplot2::labs(
      x = "Estimated Indicator (SPR or B/B0)", 
      y = "Catch Multiplier",
      color = "HCR type"
    )
}