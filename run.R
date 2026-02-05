rm(list = ls())
options(error = NULL)
library(ggplot2)
library(TMB)
library(tmbstan)
library(dplyr)
library(bayesplot)
library(here)

# list and source functions
fun_list = list.files(here("R"))
for (i in 1:length(fun_list)) {
  source(here("R", fun_list[i]))
}

# set parameters
n_sims       = 50
nyears       = 60
nages        = 10

init_nya = c(492, 492, 80, 45, 118, 156, 86, 13, 15, 68)/2

maturity     = c(0,0,0,0,0.2,0.5,0.8,1,1,1)
waa          = c(0, 0, 0, 65, 86, 100, 117, 130, 150, 152)
recruitment  = init_nya[1]
threshold    = 5000
harvest_rate = 0.15


vb_params = 
  list(
    k = 0.05, k_sd = 0.005,
    linf = 200, linf_sd = 30,
    t0 = -1
  )

sr_params = list(
  h_low            = 0.3,
  h_high           = 0.7,
  R0               = init_nya[1],
  SSB0             = sum(init_nya * waa * maturity)*4,
  sigma_rec        = 0.5,
  recruitment_mean = init_nya[1]/2,
  recruitment_sd   = (init_nya[1]/2)*0.2
  
)

# run MSE
mse_output = 
  run_mse(n_sims            = n_sims,
          nyears            = nyears,
          nages             = nages,
          init_nya          = init_nya/2,
          waa               = waa,
          rec_regime_length = 15,
          rec_type          = "BV",
          harvest_rate      = harvest_rate,
          maturity          = maturity,
          threshold         = threshold,
          vb_params         = vb_params,
          sr_params         = sr_params,
          plot              = TRUE)










