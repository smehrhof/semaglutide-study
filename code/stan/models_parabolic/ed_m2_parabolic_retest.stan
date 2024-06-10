// #include /pre/license.stan

data {
  int<lower=1> N; // number of subjects
  int<lower=1> T; // max number of trials per subject
  int<lower=1> N_time; // number of testing sessions
  array[N, N_time] int<lower=1, upper=T> T_subj; // number of trials per subject per session
  array[N, N_time, T] real<lower=0> effort_a; // effort of offer
  array[N, N_time, T] real<lower=0> amount_a; // reward of offer
  array[N, N_time, T] real<lower=0> effort_b; // effort of reject option
  array[N, N_time, T] real<lower=0> amount_b; // reward of reject option
  array[N, N_time, T] int<lower=0, upper=1> choice; // 0 for rejected offer, 1 for accepted offer
}

transformed data {
}

parameters {

// Group-level correlation matrix (cholesky factor for faster computation)
cholesky_factor_corr[N_time] L_R_kE;
cholesky_factor_corr[N_time] L_R_kR;

// Group-level parameter
// Means
vector[N_time] kE_mean;
vector[N_time] kR_mean;

// SD
vector<lower=0>[N_time] kE_sd;
vector<lower=0>[N_time] kR_sd;

// Subject-level parameters
matrix[N_time, N] kE_pr; 
matrix[N_time, N] kR_pr; 

}

transformed parameters {
  // Individual-level parameter off-sets (for non-centered parameterization)
  matrix[N_time, N] kE_tilde;
  matrix[N_time, N] kR_tilde;

  // Individual-level parameters
  matrix<lower=0, upper=1>  [N, N_time] kE;
  matrix<lower=0, upper=6>  [N, N_time] kR;

  // Construct individual offsets (for non-centered parameterization)
  kE_tilde  = diag_pre_multiply(kE_sd, L_R_kE)    * kE_pr;
  kR_tilde  = diag_pre_multiply(kR_sd, L_R_kR)    * kR_pr; 

  // Compute individual-level parameters from non-centered parameterization
  for (time in 1:N_time){
    for (i in 1:N) {
      kE[i, time]   = Phi_approx(kE_mean[time] + kE_tilde[time,i]) * 1;
      kR[i, time]   = Phi_approx(kR_mean[time] + kR_tilde[time,i]) * 6;
    }
  }
}

model {

  // Prior on cholesky factor of correlation matrix
  L_R_kE   ~ lkj_corr_cholesky(1);
  L_R_kR   ~ lkj_corr_cholesky(1); 
  
  // Group-level priors (hyperparameters)
  kE_mean   ~ normal(0, 1);
  kR_mean   ~ normal(0, 1);

  kE_sd     ~ normal(0, 0.2);
  kR_sd     ~ normal(0, 0.2);

  // Subject-level priors
  to_vector(kE_pr)   ~ normal(0, 1);
  to_vector(kR_pr)   ~ normal(0, 1);

  for(time in 1:N_time){
    for(i in 1:N){
        // Define values
        real sv_a;
        real sv_b;

        for(t in 1:(T_subj[i, time])){
            sv_a   = (kR[i, time] * amount_a[i, time, t])  - (kE[i, time] * (effort_a[i, time, t]^2));
            sv_b   = (kR[i, time] * amount_b[i, time, t])  - (kE[i, time] * (effort_b[i, time, t]^2));
            choice[i, time, t] ~ bernoulli_logit(1 * (0 + (sv_a - sv_b)));
        }
    }
  }
}

generated quantities {
  // test-retest correlations
  corr_matrix[N_time] R_kE;
  corr_matrix[N_time] R_kR;

  // For log-likelihood calculation
  array [N, N_time, T] real post_pred;
  array [N, N_time, T] real log_lik;

  // Reconstruct correlation matrices from cholesky factor
  R_kE     = L_R_kE    * L_R_kE';
  R_kR  = L_R_kR * L_R_kR';

  // initialize LL and post_pred arrays to -1
  for (i in 1:N) {
    post_pred[i,,] = rep_array(-1, N_time, T);
    log_lik[i,,] = rep_array(-1.0, N_time, T);
  }

  { // local section, this saves time and space
  for (time in 1:N_time){
    for (i in 1:N) {
      // Define values
      real sv_a;
      real sv_b;

      for (t in 1:(T_subj[i, time])) {
        log_lik[i, time, t] = 0;

        sv_a   = (kR[i, time] * amount_a[i, time, t])  - (kE[i, time] * (effort_a[i, time, t]^2));
        sv_b   = (kR[i, time] * amount_b[i, time, t])  - (kE[i, time] * (effort_b[i, time, t]^2));

        log_lik[i, time, t] += bernoulli_logit_lpmf(choice[i, time, t] | 1 * (0 + (sv_a - sv_b)));

        post_pred[i, time, t] = bernoulli_rng(inv_logit(1 * (0 + (sv_a - sv_b))));
      }
    }
  }
  }

}
