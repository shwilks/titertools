
functions {
#include /functions/normal_int_censored_likelihood.stan
}

data {
  real mu_prior_mu;
  real mu_prior_sigma;
  real sigma_prior_alpha;
  real sigma_prior_beta;
}

parameters {
  real mu;
  real<lower=0.01> sigma;
}

model {

  // Calculate likelihood of priors
  mu ~ normal(mu_prior_mu, mu_prior_sigma);
  sigma ~ inv_gamma(sigma_prior_alpha, sigma_prior_beta);

}
