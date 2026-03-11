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
threshold_list = list(
  c(0.2, 0.4), c(0.2, 0.4), c(0.2, 0.4),
  c(0.3, 0.7)
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


# define scenarios ------------------------------------
scenarios = data.frame(
  scenario_name = c("SCA_Absolute", "SPM_Absolute", "SSCL_Absolute", "LBSPR_Proportional"),
  est_method    = c("SCA", "SPM", "SS_CL", "LBSPR"),
  hcr_type      = c("absolute_hockey_stick", "absolute_hockey_stick", "absolute_hockey_stick", "proportional_multiplier"),
  max_hr        = c(0.2, 0.2, 0.2, 1.05) # 0.2 biological rate for absolute, 1.05 catch multiplier for relative
)

# run -------------------------------------------------
all_mse_results = list()
for (i in 1:nrow(scenarios)) {
  
  current_scenario = scenarios[i, ]
  
  message("\n=======================================================")
  message("Executing Scenario: ", current_scenario$scenario_name)
  message("Model: ", current_scenario$est_method, " | HCR: ", current_scenario$hcr_type)
  message("=======================================================\n")
  
  mse_output = 
    run_mse(n_sims            = 100,
            nyears            = 60,
            burn_in_length    = 40,
            hist_harvest_rate = 0.05,
            nages             = nages,
            init_nya          = stable_nya,
            waa               = waa,
            selectivity       = selectivity,
            rec_regime_length = 5,
            rec_type          = "BV",
            survival_mean     = surv,
            survival_sd       = 0.05,
            max_harvest_rate  = current_scenario$max_hr,
            maturity          = maturity,
            threshold         = threshold_list[[i]],
            vb_params         = vb_params,
            sr_params         = sr_params,
            plot              = FALSE,
            sca_model_path    = here("estimation", "SCA", "SCA_log.stan"),
            spm_model_path    = here("estimation", "SPM", "SPM_log.stan"),
            lbspr_model_path  = here("estimation", "LBSPR", "LBSPR.stan"),
            sscl_model_path   = here("estimation", "SS_CL", "SS_CL.stan"),
            mcmc_setup        = mcmc_setup,
            estimation        = TRUE,
            est_method        = current_scenario$est_method,
            hcr_type          = current_scenario$hcr_type,
            scenario_name     = current_scenario$scenario_name,
            parallel          = TRUE) 
  
  save_path = here(paste0(current_scenario$scenario_name, "_results.rds"))
  saveRDS(mse_output, file = save_path)
  message("--> Checkpoint saved: ", save_path)
  
  all_mse_results[[current_scenario$scenario_name]] = mse_output
}

message("\nAll Scenarios Completed Successfully!")