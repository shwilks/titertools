
functions {
#include /functions/normal_int_censored_likelihood.stan
}

data {
  int<lower=1> N;
  int<lower=1> N_ags;
  int<lower=1> N_srs;
  vector[N] upper_lims;
  vector[N] lower_lims;
  int ag[N];
  int sr[N];
  real mu_prior_mu;
  real mu_prior_sigma;
  real sigma_prior_alpha;
  real sigma_prior_beta;
  real sigma_prior_ag_effect;
  real sigma_prior_sr_effect;
}

parameters {
  real mu;
  real<lower=0> sigma;
  vector[N_ags] ag_effects;
  vector[N_srs] sr_effects;
}

model {

  // Calculate likelihood of priors
  mu ~ normal(mu_prior_mu, mu_prior_sigma);
  sigma ~ inv_gamma(sigma_prior_alpha, sigma_prior_beta);
  ag_effects ~ normal(0, sigma_prior_ag_effect);
  sr_effects ~ normal(0, sigma_prior_sr_effect);

  // Work out likelihood of each titer
  for (i in 1:N) {

    // Add in mixed effect for each antigen and serum
    real logtiter = mu + ag_effects[ag[i]] + sr_effects[sr[i]];

    // Work out titer likelihood
    target += normal_int_censored_likelihood(
      lower_lims[i],
      upper_lims[i],
      logtiter, sigma
    );

  }

}
