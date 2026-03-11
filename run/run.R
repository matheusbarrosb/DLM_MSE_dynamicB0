rm(list = ls())
set.seed(444)
options(error = NULL)
library(ggplot2)
library(TMB)
library(tmbstan)
library(dplyr)
library(bayesplot)
library(here)
library(rstan)
library(doParallel)
library(foreach)

# list and source functions
fun_list = list.files(here("R"))
for (i in 1:length(fun_list)) {
  source(here("R", fun_list[i]))
}

# PARAMETERS

# biology ---------------------------------------------
nages        = 10
maturity     = c(0,0,0,0,0.2,0.5,0.8,1,1,1)
waa          = c(0, 0, 0, 65, 86, 100, 117, 130, 150, 152)
selectivity  = c(0, 0.2, 0.4, 0.6, 0.9, 1, 1, 1, 1, 1)

vb_params = 
  list(
    k       = 0.15, 
    k_sd    = 0.01,
    linf    = 30,
    linf_sd = 2,
    t0      = -1
  )

M        = 0.2
surv     = exp(-M)
r0_value = 1000 

# initial stable age distribution 
stable_nya    = numeric(nages)
stable_nya[1] = r0_value
for(a in 2:(nages-1)) stable_nya[a] = stable_nya[a-1] * surv
stable_nya[nages] = (stable_nya[nages-1] * surv) / (1 - surv) # + group

# expected SSB0
ssb0_value = sum(stable_nya * waa * maturity)

sr_params = list(
  h_low            = 0.3,
  h_high           = 0.7,
  R0               = r0_value,
  SSB0             = ssb0_value, 
  sigma_rec        = 0.25,
  recruitment_mean = r0_value / 3,
  recruitment_sd   = (r0_value / 2) * 0.2
)


# management procedure --------------------------------
thresholds = c(0.25, 0.8)

hcr_params = list(
  k_down          = 50, # steepness of the catch reduction when SPR drops below the target. # higher = more aggressive braking (cliff-like drop)
  k_up            = 100, # steepness of the catch increase when SPR goes above the target. lower = more cautious, slower quota increases during good years
  gamma           = -1, # sensitivity to stock decline. higher = massive quota penalties if the stock trend (slope) is negative
  trend_years     = 5, 
  anchor_fraction = 0.05, # minimum allowed catch floor (as a fraction of historical average catch).
  thresholds      = thresholds
)

plot_hcr(
  thresholds    = thresholds, 
  max_val       = 1.05
)

# assessment model ------------------------------------
mcmc_setup = 
  list(
    chains        = 1,
    niter         = 100,
    nwarmup       = 50,
    thin          = 1,
    adapt_delta   = 0.99,
    max_treedepth = 15,
    verbose       = 0
  )


# run MSE ---------------------------------------------
mse_output = 
  run_mse(n_sims            = 100,
          nyears            = 60,
          burn_in_length    = 40,
          hist_harvest_rate = 0.05,
          nages             = 10,
          init_nya          = stable_nya,
          waa               = waa,
          selectivity       = selectivity,
          rec_regime_length = 5,
          rec_type          = "BV",
          survival_mean     = surv,
          survival_sd       = 0.05,
          max_harvest_rate  = 1.05,
          maturity          = maturity,
          threshold         = thresholds,
          vb_params         = vb_params,
          sr_params         = sr_params,
          plot              = TRUE,
          sca_model_path    = here("estimation", "SCA", "SCA_log.stan"),
          spm_model_path    = here("estimation", "SPM", "SPM_log.stan"),
          lbspr_model_path  = here("estimation", "LBSPR", "LBSPR.stan"),
          sscl_model_path   = here("estimation", "SS_CL", "SS_CL.stan"),
          mcmc_setup        = mcmc_setup,
          estimation        = TRUE,
          est_method        = "SPM",
          hcr_type          = "absolute_hockey_stick",
          parallel          = TRUE) 
# 
# 
# nages = 10
# nyears = 25
# nya = matrix(NA, nrow = nyears, ncol = nages);nya
# 
# 
# lwa = 0.0001
# lwb = 3.12
# age = 1:10
# waa = lwa*age^lwb * 1200
# 
# h = 0.6
# M = 0.2
# R0 = 500
# init_nya = numeric(nages)
# init_nya[1] = R0/2
# for (a in 2:(nages-1)) init_nya[a] = init_nya[a-1] * exp(-M)
# init_nya[nages] = (init_nya[a-1] * exp(-M))/(1 - exp(-M))
# 
# nya[1,] = init_nya
# 
# sigma = 0.1
# for (t in 2:nyears) {
#   # calculate recruitment
#   nya[t,1] = (4 * h * R0 * sum(nya[t-1,4:10]*waa[4:10])) / (SSB0 * (1 - h) + sum(nya[t-1,4:10]*waa[4:10]) * (5 * h - 1))
#   for (a in 2:nages) {
#     if (a < nages) {
#     nya[t, a] = nya[t-1, a-1] * exp(- (M + rnorm(1,0,sigma))) 
#     } else {
#       nya[t, a] = nya[t-1, a-1] * exp(- (M + rnorm(1,0,sigma))) + nya[t-1, a] * exp(- (M + rnorm(1,0,sigma))) 
#     }
#   }
# }
# 
# ssb = rowSums(nya[,4:10])
# 
# plot(ssb, type = "l")

# maybe I can use parameters of a well studied species. maybe multiple species?





















