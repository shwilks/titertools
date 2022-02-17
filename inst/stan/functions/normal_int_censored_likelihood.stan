
real normal_int_censored_likelihood(real lower_lim, real upper_lim, real mu, real sigma) {

  real result;

  if (is_inf(lower_lim) && is_inf(upper_lim)) {

    result = 0;

  } else if (lower_lim == upper_lim) {

    result = normal_lpdf(
      lower_lim | mu, sigma
    );

  } else if (!is_inf(lower_lim) && !is_inf(upper_lim)) {

    result = log_diff_exp(
      normal_lcdf(upper_lim | mu, sigma),
      normal_lcdf(lower_lim | mu, sigma)
    );

  } else if (!is_inf(lower_lim) && is_inf(upper_lim)) {

    result = normal_lccdf(
      lower_lim | mu, sigma
    );

  } else if (is_inf(lower_lim) && !is_inf(upper_lim)) {

    result = normal_lcdf(
      upper_lim | mu, sigma
    );

  }

  return result;

}
