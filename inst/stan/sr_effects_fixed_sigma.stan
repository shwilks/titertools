
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
  real ag_mean_prior_mean;
  real ag_mean_prior_sigma;
  real<lower=0.01> sigma;
  real sigma_prior_sr_effect;
}

parameters {
  vector[N_ags] ag_means;
  vector[N_srs] sr_effects;
}

model {

  // Calculate likelihood of priors
  sr_effects ~ normal(0, sigma_prior_sr_effect);
  ag_means ~ normal(ag_mean_prior_mean, ag_mean_prior_sigma);

  // Work out likelihood of each titer
  for (i in 1:N) {

    // Add in mixed effect for each antigen and serum
    real logtiter = ag_means[ag[i]] + sr_effects[sr[i]];

    // Work out titer likelihood
    target += normal_int_censored_likelihood(
      lower_lims[i],
      upper_lims[i],
      logtiter, sigma
    );

  }

}
