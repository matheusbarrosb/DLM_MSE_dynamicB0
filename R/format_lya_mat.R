format_lya_mat = function(lya_mat, vb_params) {
  
  obs_bins = seq(0, vb_params$linf, length.out = ncol(lya_mat)+1)
  obs_mids = (obs_bins[-1] + obs_bins[-length(obs_bins)])/2
  
  mean_lengths = apply(lya_mat, 1, function(row) {
    weighted.mean(obs_mids, w = row)
  })
  
  sd_lengths = apply(lya_mat, 1, function (row) {
    mu = weighted.mean(obs_mids, w = row)
    var = sum(row * (obs_mids - mu)^2) / sum(row)
    return(sqrt(var))
  })
  
  out = 
    data.frame(
      year = 1:nrow(lya_mat),
      mean_lengths = mean_lengths,
      sd_lengths = sd_lengths
    )
  
  rownames(out) = NULL
  
  out
  
}
