make_transition_matrix = function(k, linf, sigma,
                                  bin_size) {
  
  bins = seq(0, linf, bin_size)
  mids = bins[-length(bins)] + bin_size / 2
  n_bins = length(mids)
  
  transition = matrix(0, nrow = n_bins, ncol = n_bins)
  dimnames(transition) = list(mids, mids)
  
  # fill in
  for (i in 1:n_bins) {
    
    Lcurr = mids[i]
    Lexp  = Lcurr + (linf - Lcurr) * (1 - exp(-k))
    Lexp  = max(Lcurr, Lexp) 
    
    probs = diff(
      pnorm(bins, mean = Lexp, sd = sigma)
    )
   
    transition[i, ] = probs / sum(probs)
     
  }
  
  return(transition)
  
}
