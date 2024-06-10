// #include /pre/license.stan

data {
  int<lower=1> N; // number of subjects
  int<lower=1> T; // max number of trials per subject
  array[N] int<lower=1, upper=T> Tsubj; // number of trials per subject
  array[N, T] real<lower=0> effort_a; // effort of offer
  array[N, T] real<lower=0> amount_a; // reward of offer
  array[N, T] real<lower=0> effort_b; // effort of reject option
  array[N, T] real<lower=0> amount_b; // reward of reject option
  array[N, T] int<lower=0, upper=1> choice; // 0 for rejected offer, 1 for accepted offer
}

transformed data {
}

parameters {
// Declare all parameters as vectors for vectorizing

  // Group-level parameters
  vector[2] mu_pr;
  vector<lower=0>[2] sigma;

  // Subject-level parameters 
  vector[N] kE_pr;
  vector[N] a_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0, upper=1>[N] kE;
  vector<lower=-10, upper=10>[N] a;

  for (i in 1:N) {
    kE[i]    = (Phi_approx(mu_pr[1] + sigma[1] * kE_pr[i]));
    a[i]    = (Phi_approx(mu_pr[2] + sigma[2] * a_pr[i]) * 20) - 10;
  }
}

model {
// Priors
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // individual parameters
  kE_pr    ~ normal(0, 1);
  a_pr    ~ normal(0, 1);

  for (i in 1:N) {
    // Define values
    real sv_a;
    real sv_b;

    for (t in 1:(Tsubj[i])) {
      sv_a   = amount_a[i, t]  - (kE[i] * effort_a[i, t]^2);
      sv_b   = amount_b[i, t]  - (kE[i] * effort_b[i, t]^2);
      choice[i, t] ~ bernoulli_logit(a[i] + (sv_a - sv_b));
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_kE;
  real<lower=-10, upper=10> mu_a;

  // For log likelihood calculation
  array [N] real log_lik;

  // For posterior predictive check
  array [N, T] real y_pred;

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }

  mu_kE    = (Phi_approx(mu_pr[1]));
  mu_a    = (Phi_approx(mu_pr[2]) * 20) - 10;

  { // local section, this saves time and space
    for (i in 1:N) {
      // Define values
      real sv_a;
      real sv_b;

      log_lik[i] = 0;

      for (t in 1:(Tsubj[i])) {
        sv_a   = amount_a[i, t]  - (kE[i] * effort_a[i, t]^2);
        sv_b   = amount_b[i, t]  - (kE[i] * effort_b[i, t]^2);
        log_lik[i] += bernoulli_logit_lpmf(choice[i, t] | a[i] + (sv_a - sv_b));

        // generate posterior prediction for current trial
        y_pred[i, t] = bernoulli_rng(inv_logit(a[i] + (sv_a - sv_b)));
      }
    }
  }
}
