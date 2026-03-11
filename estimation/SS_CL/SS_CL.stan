data {
  int<lower=1> nyears;
  int<lower=1> n_ages;
  int<lower=1> n_len;
  
  vector[nyears] catches;
  int<lower=0> obs_len_counts[nyears, n_len]; 
  vector[n_len] len_mids;
  
  // fixed biology
  real<lower=0> M;
  real<lower=0> k;
  real<lower=0> Linf;
  real t0;
  real<lower=0> CV_L;
  vector[n_ages] waa;
  vector[n_ages] mat;
  
  // stock-recruit
  real<lower=0> h;
  real<lower=0> R0;
  real<lower=0> SSB0;
  
  // error terms
  real<lower=0> sigma_R;
  real<lower=0> sigma_C;
}

transformed data {
  matrix[n_ages, n_len] P_la;
  vector[n_ages] N0;
  
  // pre-calculate age-length key
  for (a in 1:n_ages) {
    real mu_L = Linf * (1.0 - exp(-k * (a - t0)));
    real sigma_L = mu_L * CV_L;
    for (l in 1:n_len) {
      P_la[a, l] = exp(normal_lpdf(len_mids[l] | mu_L, sigma_L));
    }
    P_la[a, ] = P_la[a, ] / sum(P_la[a, ]);
  }
  
  // initialize unfished state
  N0[1] = R0;
  for (a in 2:(n_ages-1)) {
    N0[a] = N0[a-1] * exp(-M);
  }
  N0[n_ages] = (N0[n_ages-1] * exp(-M)) / (1.0 - exp(-M));
}

parameters {
  vector<lower=-10, upper=2>[nyears] log_F;
  vector[nyears-1] rec_devs;
  
  // length-based selectivity
  real<lower=0> SL50;
  real<lower=0.1> SL_diff;
}

transformed parameters {
  vector[n_len] S_len;
  vector[n_ages] S_age;
  matrix[nyears, n_ages] N;
  matrix[nyears, n_ages] C_num;
  vector[nyears] pred_catch;
  matrix[nyears, n_len] pred_len_prop;
  vector[nyears] SSB;
  vector[nyears] F_rate = exp(log_F);
  
  // calculate selectivity
  for (l in 1:n_len) {
    S_len[l] = 1.0 / (1.0 + exp(-2.944 * (len_mids[l] - SL50) / SL_diff));
  }
  S_age = P_la * S_len;
  
  // population dynamics
  N[1, ] = N0'; 
  
  for (t in 1:nyears) {
    
    // calculate SSB
    SSB[t] = sum(N[t, ]' .* waa .* mat);
    
    // recruitment for next year
    if (t < nyears) {
      real exp_R = (4 * h * R0 * SSB[t]) / (SSB0 * (1 - h) + SSB[t] * (5 * h - 1));
      N[t+1, 1]  = exp_R * exp(rec_devs[t] - 0.5 * sigma_R^2);
    }
    
    // mortality and catch
    for (a in 1:n_ages) {
      real Z = M + F_rate[t] * S_age[a];
      
      // survivors to next year
      if (t < nyears) {
        if (a < n_ages) {
          N[t+1, a+1] = N[t, a] * exp(-Z);
        } else {
          N[t+1, n_ages] += N[t, a] * exp(-Z); // plus group
        }
      }
      
      // baranov catch equation
      if (Z > 1e-5) {
        C_num[t, a] = (F_rate[t] * S_age[a] / Z) * N[t, a] * (1.0 - exp(-Z));
      } else {
        C_num[t, a] = 0.0;
      }
    }
    
    // aggregate catch biomass
    pred_catch[t] = sum(C_num[t, ]' .* waa);
    
    // calculate length proportions
    pred_len_prop[t, ] = (P_la' * C_num[t, ]')';
    pred_len_prop[t, ] = (pred_len_prop[t, ] + 1e-10) / sum(pred_len_prop[t, ] + 1e-10);
  }
}

model {
  // priors
  log_F    ~ normal(-2, 2);
  rec_devs ~ normal(0, sigma_R);
  SL50     ~ normal(0.5 * Linf, 0.2 * Linf);
  SL_diff  ~ normal(0.1 * Linf, 0.05 * Linf);
  
  // likelihoods
  for (t in 1:nyears) {
    // catch penalty (lognormal)
    catches[t] ~ lognormal(log(pred_catch[t]), sigma_C);
    
    // length composition (multinomial)
    obs_len_counts[t, ] ~ multinomial(pred_len_prop[t, ]');
  }
}

generated quantities {
  real curr_status = SSB[nyears] / SSB0;
}