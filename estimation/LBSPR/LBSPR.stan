data {
  int<lower=1> n_len;
  int<lower=1> n_ages;
  vector[n_len] len_mids;
  int<lower=0> obs_counts[n_len]; 
  
  // fixed parameters
  real<lower=0> M;
  real<lower=0> k;
  real<lower=0> Linf;
  real<lower=0> CV_L;
  vector[n_ages] maturity_at_age;
  vector[n_ages] weight_at_age;
  
  // fixed selectivity from OM
  vector[n_ages] fixed_selectivity;
}

transformed data {
  matrix[n_ages, n_len] P_la;
  
  // pre-calculate growth kernel
  for (a in 1:n_ages) {
    real mu_L = Linf * (1.0 - exp(-k * a));
    real sigma_L = mu_L * CV_L;
    for (l in 1:n_len) {
      P_la[a, l] = exp(normal_lpdf(len_mids[l] | mu_L, sigma_L));
    }
    P_la[a, ] = P_la[a, ] / sum(P_la[a, ]);
  }
}

parameters {
  real<lower=0, upper=5> FM_ratio; 
}

transformed parameters {
  vector[n_ages] N_eq;
  vector[n_ages] C_age;
  vector[n_len]  pred_len_prop;
  
  real F_rate = M * FM_ratio;
  
  // ----------------------------------------------------
  // 1. Equilibrium Dynamics WITH PLUS GROUP
  // ----------------------------------------------------
  N_eq[1] = 1.0; 
  
  // Normal ages (up to n_ages - 1)
  for (a in 2:(n_ages-1)) {
    real Z_prev = M + F_rate * fixed_selectivity[a-1];
    N_eq[a] = N_eq[a-1] * exp(-Z_prev);
  }
  
  // Plus Group (Age n_ages)
  // N_plus = N_incoming / (1 - exp(-Z_plus))
  {
    real Z_prev = M + F_rate * fixed_selectivity[n_ages-1];
    real Z_plus = M + F_rate * fixed_selectivity[n_ages];
    
    real N_incoming = N_eq[n_ages-1] * exp(-Z_prev);
    
    // Geometric series sum for the tail
    N_eq[n_ages] = N_incoming / (1.0 - exp(-Z_plus));
  }
  
  // ----------------------------------------------------
  // 2. Catch at Age
  // ----------------------------------------------------
  for (a in 1:n_ages) {
    real Z = M + F_rate * fixed_selectivity[a];
    
    // Baranov catch equation
    if (Z > 1e-5) 
      C_age[a] = (F_rate * fixed_selectivity[a] / Z) * N_eq[a] * (1.0 - exp(-Z));
    else 
      C_age[a] = 0;
  }
  
  // predicted length comps
  pred_len_prop = P_la' * C_age;
  pred_len_prop = (pred_len_prop + 1e-10) / sum(pred_len_prop + 1e-10);
}

model {
  // priors
  FM_ratio ~ lognormal(log(1), 0.5); 
  
  // likelihood
  obs_counts ~ multinomial(pred_len_prop);
}

generated quantities {
  real SPR;
  {
    real SSB_F0 = 0;
    real SSB_Fcurr = 0;
    
    vector[n_ages] N_0;
  
    N_0[1] = 1.0;
    for(a in 2:(n_ages-1)) {
      N_0[a] = N_0[a-1] * exp(-M);
    }
    
    // Plus group for unfished state
    N_0[n_ages] = (N_0[n_ages-1] * exp(-M)) / (1.0 - exp(-M));
    
    for(a in 1:n_ages) {
      SSB_F0    += N_0[a] * maturity_at_age[a] * weight_at_age[a];
      SSB_Fcurr += N_eq[a] * maturity_at_age[a] * weight_at_age[a];
    }
    SPR = SSB_Fcurr / SSB_F0;
  }
}
