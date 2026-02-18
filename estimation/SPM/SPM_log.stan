data {
  int<lower=1> nyears;
  int<lower=1> nages;
  vector[nyears] catches;
  vector[nyears] index;
}
parameters {
  real log_r;
  real log_k;
  real log_q;
  real<lower=0> sigma_process;
  real<lower=0> sigma_obs;
  vector[nyears] log_B_devs;            
}
transformed parameters {
  real r = exp(log_r);
  real K = exp(log_k);
  real q = exp(log_q);
  vector[nyears] B;
  B[1] = K;
  for (t in 2:nyears) {
    real surplus_prod = r * B[t-1] * (1 - B[t-1]/K);
    real Btemp        = B[t-1] + surplus_prod - catches[t-1];
    if(Btemp < 1e-1) Btemp = 1e-1;
    B[t] = Btemp * exp(log_B_devs[t-1]);
  }
}
model {
  // priors
  log_r ~ normal(-1.5, 0.5);
  log_k ~ normal(12.5, 2);
  log_q ~ normal(-7, 2);
  
  sigma_process ~ normal(0, 0.2);
  sigma_obs     ~ normal(0, 0.2);

  log_B_devs ~ std_normal();
  
  // likelihood
  index ~ lognormal(log(B*q), sigma_obs);
}
generated quantities {
  real MSY  = (r*K)/4;
  real Bmsy = K/2;
  real Fmsy = r/2;
  real B_B0 = B[nyears]/K;
}