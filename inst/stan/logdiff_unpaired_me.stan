
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
  int ag12[N];
  int sr[N];
  real ags1_mu_prior_mu;
  real ags1_mu_prior_sigma;
  real sigma_prior_alpha;
  real sigma_prior_beta;
  real sigma_prior_ag_effect;
  real sr_logdiffs_mu;
  real sr_logdiffs_sigma;
}

parameters {
  real ags1_mu;
  real<lower=0> sigma;
  vector[N_ags] ag_effects;
  vector[N_srs] sr_logdiffs;
}

model {

  // Calculate likelihood of priors
  ags1_mu ~ normal(ags1_mu_prior_mu, ags1_mu_prior_sigma);
  sigma ~ inv_gamma(sigma_prior_alpha, sigma_prior_beta);
  ag_effects ~ normal(0, sigma_prior_ag_effect);
  sr_logdiffs ~ normal(sr_logdiffs_mu, sr_logdiffs_sigma);

  // Work out likelihood of each titer
  for (i in 1:N) {

    // Estimate log titer
    real logtiter;
    if (ag12[i] == 1) {
      logtiter = ags1_mu + ag_effects[ag[i]];
    } else {
      logtiter = ags1_mu + ag_effects[ag[i]] + sr_logdiffs[sr[i]];
    }

    // Work out titer likelihood
    target += normal_int_censored_likelihood(
      lower_lims[i],
      upper_lims[i],
      logtiter, sigma
    );

  }

}
