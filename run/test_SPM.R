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

# PARAMETERS ------------------------------------------

# biology 
nages        = 10
maturity     = c(0,0,0,0,0.2,0.5,0.8,1,1,1)
waa          = c(5, 18, 35, 65, 86, 100, 117, 130, 150, 152) # Non-zero weights
selectivity  = maturity

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
stable_nya[nages] = (stable_nya[nages-1] * surv) / (1 - surv) 

ssb0_value = sum(stable_nya * waa * maturity)

sr_params = list(
  h_low            = 0.3,
  h_high           = 0.7,
  R0               = r0_value,
  SSB0             = ssb0_value, 
  sigma_rec        = 0.2,
  recruitment_mean = r0_value / 3,
  recruitment_sd   = (r0_value / 2) * 0.2
)

# configure assessment model --------------------------
mcmc_setup = 
  list(
    chains        = 1,
    niter         = 1000,
    nwarmup       = 500,
    thin          = 1,
    adapt_delta   = 0.99,
    max_treedepth = 15,
    verbose       = 0
  )

# run MSE (Generating Historical Data) ----------------
mse_output = 
  run_mse(n_sims            = 1,
          nyears            = 40,
          burn_in_length    = 40,
          hist_harvest_rate = 0.2,    # Strong enough to generate contrast
          nages             = nages,
          init_nya          = stable_nya,
          waa               = waa,
          selectivity       = selectivity,
          rec_regime_length = 15,     # Balanced regime
          rec_type          = "BV",
          survival_mean     = surv,
          survival_sd       = 0.05,
          max_harvest_rate  = 0.15,
          maturity          = maturity,
          thresholds        = c(0.2, 0.7),
          vb_params         = vb_params,
          sr_params         = sr_params,
          plot              = FALSE,
          estimation        = FALSE)

# test SPM model --------------------------------------
model_dir = here("estimation", "SPM")
stanc(file.path(model_dir, "SPM_log.stan")) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
spm_model = stan_model(file = file.path(model_dir, "SPM_log.stan"))

# grab data
sim_1 = mse_output$raw_sims[[1]]

stan_input = make_historical_spm_data(
  sim            = sim_1, 
  burn_in_length = 40,    
  waa            = waa,
  selectivity    = selectivity,
  sigma_index    = 0.1 
);str(stan_input)

init_fun = function() {
  list(
    log_r = log(0.1),                             
    log_k = log(ssb0_value),           
    log_q = log(0.001)           
  )
}

fit = sampling(
  object  = spm_model,
  data    = stan_input,        
  init    = init_fun,        
  chains  = 1,                
  iter    = 200,            
  warmup  = 50,            
  thin    = 1,                
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

post = rstan::extract(fit)

abs_B_matrix = post$B

est_b_median = apply(abs_B_matrix, 2, median)
est_b_lower  = apply(abs_B_matrix, 2, quantile, probs = 0.025)
est_b_upper  = apply(abs_B_matrix, 2, quantile, probs = 0.975)

true_b = numeric(stan_input$nyears)
for(y in 1:stan_input$nyears) {
  true_b[y] = sum(sim_1$nya_mat[y, ] * waa * selectivity)
}

df_plot = data.frame(
  Year      = 1:stan_input$nyears,
  True_B    = true_b/true_b[1],
  Est_Med   = est_b_median/est_b_median[1],
  Est_Low   = est_b_lower/est_b_lower[1],
  Est_High  = est_b_upper/est_b_upper[1]
)

ggplot(df_plot, aes(x = Year)) +
  geom_ribbon(aes(ymin = Est_Low, ymax = Est_High, fill = "95% CI"), alpha = 0.3) +
  geom_line(aes(y = Est_Med, color = "Estimated Absolute"), size = 1) +
  geom_point(aes(y = True_B, color = "True"), size = 2, shape = 21) +
  scale_color_manual(values = c("Estimated" = "black", "True" = "red")) +
  scale_fill_manual(values = c("95% CI" = "grey30")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1.4)
  ) +
  ylim(0, max(c(df_plot$Est_High, df_plot$True_B), na.rm = TRUE) * 1.2)
