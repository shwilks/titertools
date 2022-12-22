
functions {
#include /functions/normal_int_censored_likelihood.stan
}

data {
  int<lower=1> N;
  vector[N] upper_lims1;
  vector[N] lower_lims1;
  vector[N] upper_lims2;
  vector[N] lower_lims2;
  real intercept_prior_mu;
  real intercept_prior_sigma;
  real slope_prior_mu;
  real slope_prior_sigma;
  real sigma_prior_alpha;
  real sigma_prior_beta;
}

parameters {
  real intercept;
  real slope;
  real<lower=0.01> sigma;
  real logtiter1mean;
  real<lower=0.01> logtiter1sd;
  vector[N] log2titers1;
}

model {

  // Calculate likelihood of priors
  intercept ~ normal(intercept_prior_mu, intercept_prior_sigma);
  slope ~ normal(slope_prior_mu, slope_prior_sigma);
  sigma ~ inv_gamma(sigma_prior_alpha, sigma_prior_beta);
  log2titers1 ~ normal(logtiter1mean, logtiter1sd);

  for (i in 1:N) {

    // Calculate likelihood of the measurements 1 given log2titers1
    target += normal_int_censored_likelihood(
      lower_lims1[i],
      upper_lims1[i],
      log2titers1[i],
      sigma1
    );

    // Calculate likelihood of the measurements 2 given log2titers1 and slope and intercept
    target += normal_int_censored_likelihood(
      lower_lims2[i],
      upper_lims2[i],
      log2titers1[i]*slope + intercept,
      sigma2
    );

  }

}
