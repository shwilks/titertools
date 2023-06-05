
functions {
#include /functions/normal_int_censored_likelihood.stan
}

data {
  int<lower=0> N;
  int<lower=0> N_ags;
  int<lower=0> N_srs;
  int<lower=0> N_sr_groups;
  vector[N] upper_logtiter_lims;
  vector[N] lower_logtiter_lims;
  vector[N_srs] upper_homologous_logtiter_lims;
  vector[N_srs] lower_homologous_logtiter_lims;
  int ags[N];
  int srs[N];
  int sr_groups[N];
}

parameters {
  vector<lower=0.01, upper=2>[N_sr_groups - 1] sr_group_slopes;
  vector[N_ags] ag_folddrops;
  real<lower=0> logtiter_error_sigma;
}

model {

  // Set variables
  real folddrop;
  real sr_group_slope;

  // Add a very flat prior to antigen folddrops to help with convergence
  ag_folddrops ~ normal(0, 1000);

  // Work out likelihood of each titer
  for (i in 1:N) {

    // Fetch antigen number and serum group number
    int ag = ags[i];
    int sr = srs[i];
    int sr_group = sr_groups[i];

    // Fetch sr slope parameter
    if (sr_group == 1) sr_group_slope = 1;
    else               sr_group_slope = sr_group_slopes[sr_group - 1];

    // Predict titer
    folddrop = ag_folddrops[ag]*sr_group_slope;

    // Accumulate likelihood of predicted titer
    target += normal_int_censored_likelihood(
      lower_logtiter_lims[i] - upper_homologous_logtiter_lims[sr],
      upper_logtiter_lims[i] - lower_homologous_logtiter_lims[sr],
      folddrop,
      logtiter_error_sigma
    );

  }

}
