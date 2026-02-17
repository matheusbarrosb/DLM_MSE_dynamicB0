#rm(list = ls())
set.seed(444)
options(error = NULL)
library(ggplot2)
library(TMB)
library(tmbstan)
library(dplyr)
library(bayesplot)
library(here)
library(rstan)

# list and source functions
fun_list = list.files(here("R"))
for (i in 1:length(fun_list)) {
  source(here("R", fun_list[i]))
}

# set parameters --------------------------------------
n_sims       = 1000
nyears       = 60
nages        = 10

init_nya = c(492, 492, 80, 45, 118, 156, 86, 13, 15, 68)

maturity     = c(0,0,0,0,0.2,0.5,0.8,1,1,1)
waa          = c(0, 0, 0, 65, 86, 100, 117, 130, 150, 152)
selectivity  = c(0, 0, 0, 0.5, 0.75, 0.8, 0.9, 1, 1, 1)
selectivity  = maturity
recruitment  = init_nya[1]
thresholds    = c(0.2, 0.7)
max_harvest_rate = 0.15


vb_params = 
  list(
    k       = 0.15, 
    k_sd    = 0.005,
    linf    = 200,
    linf_sd = 30,
    t0      = -1
  )

M = 0.2
surv = exp(-M)
r0_value = 1000 

# Build a stable age distribution (N at age)
stable_nya = numeric(nages)
stable_nya[1] = r0_value
for(a in 2:(nages-1)) stable_nya[a] = stable_nya[a-1] * surv
stable_nya[nages] = (stable_nya[nages-1] * surv) / (1 - surv) # Plus group

# expected SSB0
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
    chains = 1,
    niter   = 30,
    nwarmup = 10,
    thin   = 1,
    adapt_delta = 0.99,
    max_treedepth = 15,
    verbose = 0
  )


# run MSE ---------------------------------------------
mse_output = 
  run_mse(n_sims            = 100,
          nyears            = 40,
          burn_in_length    = 40,
          hist_harvest_rate = 0.1,
          nages             = nages,
          init_nya          = stable_nya,
          waa               = waa,
          selectivity       = selectivity,
          rec_regime_length = 15,
          rec_type          = "BV",
          survival_mean     = surv,
          survival_sd       = 0.05,
          max_harvest_rate  = max_harvest_rate,
          maturity          = maturity,
          threshold         = thresholds,
          vb_params         = vb_params,
          sr_params         = sr_params,
          plot              = TRUE,
          sca_model_path    = here("estimation", "SCA", "SCA_log.stan"),
          mcmc_setup        = mcmc_setup,
          estimation        = FALSE)

# test SCA model --------------------------------------
model_dir = here("estimation", "SCA")
stanc(file.path(model_dir, "SCA_log.stan"))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
sca_model = stan_model(file = file.path(model_dir, "SCA_log.stan"))

# grab data
sim_1 = mse_output$sims[[1]]
stan_input = make_historical_sca_data(
  sim            = sim_1, 
  burn_in_length = 40,    
  waa            = waa,
  maturity       = maturity,
  selectivity    = selectivity,
  nages          = nages,
  M              = 0.2,
  sigma_index    = 0.05,   
  ess            = 100    
);str(stan_input)

init_fun = function() {
  list(
    log_R0   = 7,                       
    h        = 0.5,
    q        = log(0.001),              
    rec_devs = rep(0,40-1)
  )
}

fit = sampling(
  object = sca_model,
  data   = stan_input,       
  init   = init_fun,        
  chains = 1,               
  iter   = 250,            
  warmup = 100,            
  thin   = 1,               
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

post = rstan::extract(fit)


est_ssb_median = apply(post$SSB, 2, median)
est_ssb_lower  = apply(post$SSB, 2, quantile, probs = 0.025)
est_ssb_upper  = apply(post$SSB, 2, quantile, probs = 0.975)

true_ssb = sim_1$biomass[1:stan_input$nyears]

df_plot = data.frame(
  Year      = 1:stan_input$nyears,
  True_SSB  = true_ssb,
  Est_Med   = est_ssb_median,
  Est_Low   = est_ssb_lower,
  Est_High  = est_ssb_upper
)

ggplot(df_plot, aes(x = Year)) +
  geom_ribbon(aes(ymin = Est_Low, ymax = Est_High, fill = "95% CI"), alpha = 0.3) +
  geom_line(aes(y = Est_Med, color = "Estimated"), size = 1) +
  geom_point(aes(y = True_SSB, color = "True"), size = 2, shape = 21) +
  scale_color_manual(values = c("Estimated" = "black", "True" = "red")) +
  scale_fill_manual(values = c("95% CI" = "grey30")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1.4),
    
  ) +
  ylim(0, max(df_plot$Est_High) * 1.2)







