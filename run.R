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

init_nya = c(492, 492, 80, 45, 118, 156, 86, 13, 15, 68)

maturity     = c(0,0,0,0,0.2,0.5,0.8,1,1,1)
waa          = c(0, 0, 0, 65, 86, 100, 117, 130, 150, 152)
selectivity  = c(0, 0, 0, 0.5, 0.75, 0.8, 0.9, 1, 1, 1)
recruitment  = init_nya[1]
threshold    = 10000
harvest_rate = 0.1


vb_params = 
  list(
    k       = 0.15, 
    k_sd    = 0.005,
    linf    = 200,
    linf_sd = 30,
    t0      = -1
  )

M = 0.2
surv = 0.8
r0_value = 1000 

# Build a stable age distribution (N at age)
stable_nya = numeric(nages)
stable_nya[1] = r0_value
for(a in 2:(nages-1)) stable_nya[a] = stable_nya[a-1] * surv
stable_nya[nages] = (stable_nya[nages-1] * surv) / (1 - surv) # Plus group

# expected SSB0
ssb0_value = sum(stable_nya * waa * maturity)


sr_params = list(
  h_low            = 0.25,
  h_high           = 0.75,
  R0               = r0_value,
  SSB0             = ssb0_value, 
  sigma_rec        = 0.1,
  recruitment_mean = r0_value / 3,
  recruitment_sd   = (r0_value / 2) * 0.2
)

# run MSE
mse_output = 
  run_mse(n_sims            = n_sims,
          nyears            = nyears+15,
          nages             = nages,
          init_nya          = stable_nya,
          waa               = waa,
          selectivity =      selectivity,
          rec_regime_length = 15,
          rec_type          = "BV",
          harvest_rate      = harvest_rate,
          maturity          = maturity,
          threshold         = threshold,
          vb_params         = vb_params,
          sr_params         = sr_params,
          plot              = TRUE)


# calculate reference points
# 1. life history equilibrium (phi_0)
surv_m = 0.8
phi_0 = sum(surv_m^(0:8) * waa[1:9] * maturity[1:9]) + 
  (surv_m^9 / (1 - surv_m)) * waa[10] * maturity[10]

get_equil_stats = function(H, h_val, R0, SSB0, phi_0, waa, maturity, nages, surv_m) {
  # phi_f: spawning biomass per recruit under fishing
  s_f = surv_m * (1 - H)
  phi_f = sum(s_f^(0:8) * waa[1:9] * maturity[1:9]) + 
    (s_f^9 / (1 - s_f)) * waa[10] * maturity[10]
  
  # R_eq: equilibrium recruitment
  # if the term in brackets is <= 0, the population cannot persist at this H
  num = 4 * h_val * (phi_f / phi_0) - (1 - h_val)
  den = (5 * h_val - 1) * (phi_f / phi_0)
  
  if (num <= 0) return(list(yield = 0, ssb = 0))
  
  R_eq = R0 * (num / den)
  
  # Yield per recruit (ypr)
  ypr = sum(s_f^(0:8) * waa[1:9] * maturity[1:9] * H) + 
    (s_f^9 / (1 - s_f)) * waa[10] * maturity[10] * H
  
  return(list(yield = R_eq * ypr, ssb = R_eq * phi_f))
}

# 2. grid search for MSY
regimes = list(Low = 0.25, High = 0.75)
ref_results = data.frame()

for(name in names(regimes)) {
  h_val = regimes[[name]]
  h_grid = seq(0, 0.5, length.out = 500)
  
  stats = lapply(h_grid, get_equil_stats, h_val=h_val, R0=sr_params$R0, 
                 SSB0=sr_params$SSB0, phi_0=phi_0, waa=waa, 
                 maturity=maturity, nages=nages, surv_m=0.8)
  
  yields = sapply(stats, function(x) x$yield)
  ssbs   = sapply(stats, function(x) x$ssb)
  
  max_idx = which.max(yields)
  
  ref_results = rbind(ref_results, data.frame(
    Regime = name,
    MSY    = round(yields[max_idx], 2),
    Bmsy   = round(ssbs[max_idx], 2),
    Hmsy   = round(h_grid[max_idx], 3),
    B0     = round(sr_params$SSB0, 2)
  ))
}

print(ref_results)