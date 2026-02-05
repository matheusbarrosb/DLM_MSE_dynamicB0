compute_recruitment = function(ssb, sr_params, rec_regime_length, year, type = c("decoupled", "BV")) {
  
  if (rec_regime_length < 5) {
    stop("rec_regime_length has to be larger than 5 years.")
  }
  # check parameter names, has to be: h, R0, SSB0, sigma_rec
  if (!all(c("h_low", "h_high", "R0", "SSB0", "sigma_rec") %in% names(sr_params))) {
    stop("sr_params has to contain the following parameters: h, R0, SSB0, sigma_rec")
  }
  if (sr_params$sigma_rec < 0) {
    stop("sigma_rec has to be non-negative.")
  }
  if (sr_params$R0 < 0) {
    stop("R0 has to be non-negative.")
  }
  if (sr_params$SSB0 < 0) {
    stop("SSB0 has to be non-negative.")
  }

  # parameters
  h_high           = sr_params$h_high
  h_low            = sr_params$h_low
  R0               = sr_params$R0
  SSB0             = sr_params$SSB0
  sigma_rec        = sr_params$sigma_rec
  recruitment_mean = sr_params$recruitment_mean
  recruitment_sd   = sr_params$recruitment_sd

  if (type == "decoupled") {
    if (((year - 1) %/% rec_regime_length)%% 2 == 0) {
      recruitment = rnorm(1, mean = recruitment_mean * 1.5, sd = recruitment_sd)
    } else {
      recruitment = rnorm(1, mean = recruitment_mean * 0.5, sd = recruitment_sd)
    }
  } else if (type == "BV") {
    if (((year - 1) %/% rec_regime_length)%% 2 == 0) {
      h = h_high
    } else {
      h = h_low
    }
    # compute recruitment
    recruitment = (4 * h * R0 * ssb) / (SSB0 * (1 - h) + ssb * (5 * h - 1))
  } else {
    stop("type has to be either 'declouped' or 'BV'.")
  }
  
  return(recruitment)
  
}











