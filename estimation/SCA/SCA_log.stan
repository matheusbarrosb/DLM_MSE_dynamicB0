data {
  int<lower=1> nyears;
  int<lower=1> nages;
  
  vector[nages] maturity;
  matrix[nyears, nages] waa; 
  real<lower=0> M;    
  vector[nyears] catches;
  vector[nages] selectivity; 
  
  vector[nyears] index;
  int<lower=0> nya_obs[nyears, nages];
  
  real<lower=0> catch_sigma;
  real<lower=0> index_sigma;
}

parameters {
  real<lower=0> log_R0;             
  real<lower=0.2, upper=1> h;       
  vector[nyears-1] rec_devs;        
  real q;                           
}

transformed parameters {
  matrix[nyears, nages] log_N;      
  vector[nyears] SSB;               
  vector[nyears] pred_catch;        
  vector[nyears] pred_index;        

  real R0 = exp(log_R0);
  real SSB0;
  vector[nages] log_npr_virgin;     
  
  // INITIALIZATION
  log_npr_virgin[1] = 0; 
  for (a in 2:(nages-1)) {
    log_npr_virgin[a] = log_npr_virgin[a-1] - M;
  }
  log_npr_virgin[nages] = log_npr_virgin[nages-1] - M - log(1 - exp(-M));
  
  SSB0 = sum(exp(log_npr_virgin) .* to_vector(waa[1]) .* maturity) * R0;
  
  for (a in 1:nages) {
    log_N[1, a] = log(R0) + log_npr_virgin[a];
  }
  
  for (y in 1:nyears) {
    
    // 1. calculate SSB (start of year)
    SSB[y] = sum(exp(row(log_N, y)) .* row(waa, y) .* to_row_vector(maturity));
    
    // 2. predict index (Start of year)
    pred_index[y] = exp(q) * sum(exp(row(log_N, y)) .* row(waa, y) .* to_row_vector(selectivity));
    
    // 3. recruitment for next year
    if (y < nyears) {
       real num = 4 * h * R0 * SSB[y];
       real den = SSB0 * (1 - h) + SSB[y] * (5 * h - 1);
       //log_N[y+1, 1] = log(num / den) + rec_devs[y] - 0.5 * 0.1^2; 
       log_N[y+1, 1] = log(num / den) + rec_devs[y] - 0.5 * 0.2^2;
    }

    // 4. calculate harvest rate U on post-M biomass
    //(survival first, then catch)
    
    real b_available_post_M = 0;
    
    for(a in 1:nages) {
      // Calculate biomass AFTER natural mortality
      b_available_post_M += exp(log_N[y,a]) * exp(-M) * waa[y,a] * selectivity[a];
    }
    
    // calculate U
    // U = catch / biomass available after M
    real U = catches[y] / (b_available_post_M + 1e-10);
    
    if (U > 0.99) U = 0.99;
    
    pred_catch[y] = catches[y]; 
    
    // 5. Project Forward
    if (y < nyears) {
      for (a in 1:nages) {
        
        // step A: apply mortality first
        real n_post_M = exp(log_N[y,a]) * exp(-M);
        
        // Step B: then catch 
        // N_next = N_post_M * (1 - U * Selectivity)
        real n_end = n_post_M * (1 - U * selectivity[a]);
        
        if (a < nages) {
           if(n_end < 1e-10) n_end = 1e-10;
           
           if (a < nages - 1) {
             log_N[y+1, a+1] = log(n_end);
           } else {
             // Plus Group
             real n_plus_post_M = exp(log_N[y,nages]) * exp(-M);
             real n_plus_end    = n_plus_post_M * (1 - U * selectivity[nages]);
             log_N[y+1, nages] = log(n_end + n_plus_end); 
           }
        }
      }
    }
  }
}

model {
  log_R0   ~ normal(5, 2);
  h        ~ beta(5, 5);
  rec_devs ~ normal(0, 0.2);
  q        ~ normal(0, 5); 
  
  index ~ lognormal(log(pred_index), index_sigma);
  
  for (y in 1:nyears) {
    vector[nages] p_raw;
    p_raw = exp(row(log_N, y))' .* selectivity; 
    nya_obs[y] ~ multinomial((p_raw + 1e-10) / sum(p_raw + 1e-10));
  }
}

generated quantities {
  real est_SSB_terminal = SSB[nyears];
  real est_B0 = SSB0;
  vector[nyears] est_recruitment = exp(col(log_N, 1));
}
