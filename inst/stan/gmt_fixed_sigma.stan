
functions {
#include /functions/normal_int_censored_likelihood.stan
}

data {
  int<lower=1> N;
  vector[N] upper_lims;
  vector[N] lower_lims;
  real mu_prior_mu;
  real mu_prior_sigma;
  real sigma_prior_alpha;
  real sigma_prior_beta;
  real<lower=0> sigma;
}

parameters {
  real mu;
}

model {

  // Calculate likelihood of priors
  mu ~ normal(mu_prior_mu, mu_prior_sigma);

  // Work out likelihood of each titer
  for (i in 1:N) {

    target += normal_int_censored_likelihood(
      lower_lims[i],
      upper_lims[i],
      mu, sigma
    );

  }

}
