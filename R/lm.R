
#' Calculate the linear relationship between two sets of titers
#'
#' This function is useful for example for calculating differences where you have a set of
#' pre-vaccination and post-vaccination samples and you would like to know the mean response
#' size, accounting for non-detectable values.
#'
#' @param titers1 The first titerset from which to calculate the mean difference
#' @param titers2 The second titerset to which to calculate the mean difference
#' @param ci_method The method to use when calculating the confidence intervals
#' @param ci_level The confidence level to use when calculating confidence intervals
#' @param dilution_stepsize The dilution stepsize. Generally 1 for titers of
#' the form '20', '40', '80', and 0 for titers that are read out continuously.
#' @param mu_prior_mu Prior for the mean of the normal distribution to estimate the mean
#' @param mu_prior_sigma Prior for the standard deviation of the normal distribution to estimate the mean
#' @param sigma_prior_alpha Prior for the alpha parameter of the inverse gamma distribution to estimate the standard deviation
#' @param sigma_prior_beta Prior for the beta parameter of the inverse gamma distribution to estimate the standard deviation
#' @param options Options for the sampler
#'
#' @export
titerlm <- function(
    titers1,
    titers2,
    ci_method = "HDI",
    ci_level = 0.95,
    dilution_stepsize = NA,
    intercept_prior_mu = 0,
    intercept_prior_sigma = 10000,
    slope_prior_mu = 0,
    slope_prior_sigma = 10000,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 5,
    return_all_pars = FALSE,
    options = list()
) {

  # Check input
  if (!is.vector(titers1) || !is.vector(titers2) || length(titers1) != length(titers2)) {
    "Titers must be input as vectors of matching length"
  }

  # Infer dilution stepsize if necessary
  if (is.na(dilution_stepsize)) dilution_stepsize <- infer_dilution_stepsize(titers1)

  # Remove na titers
  na_titers <- is_na_titer(titers1) | is_na_titer(titers2)
  titers1 <- titers1[!na_titers]
  titers2 <- titers2[!na_titers]

  if (length(titers1) == 0) {
    return(
      matrix(
        NA_real_, 2, 3,
        dimnames = list(
          c("intercept", "slope"),
          c("estimate", "lower", "upper")
        )
      )
    )
  }

  # Calculate titer limits
  titerlims1 <- calc_titer_lims(titers1, dilution_stepsize)
  titerlims2 <- calc_titer_lims(titers2, dilution_stepsize)

  # Setup data
  standata <- list(
    upper_lims1 = as.array(titerlims1$max_titers),
    lower_lims1 = as.array(titerlims1$min_titers),
    upper_lims2 = as.array(titerlims2$max_titers),
    lower_lims2 = as.array(titerlims2$min_titers),
    N = sum(!na_titers),
    intercept_prior_mu = intercept_prior_mu,
    intercept_prior_sigma = intercept_prior_sigma,
    slope_prior_mu = slope_prior_mu,
    slope_prior_sigma = slope_prior_sigma,
    sigma_prior_alpha = sigma_prior_alpha,
    sigma_prior_beta = sigma_prior_beta
  )

  # Set initial conditions based on a basic linear regression
  lmfit <- lm(titerlims2$log_titers ~ titerlims1$log_titers)
  initdata <- list(
    intercept = unname(lmfit$coefficients[1]),
    slope = unname(lmfit$coefficients[2]),
    sigma = sd(lmfit$residuals),
    logtiter1mean = mean(titerlims1$log_titers),
    logtiter1sd = sd(titerlims1$log_titers),
    log2titers1 = titerlims1$log_titers
  )

  # Optimize parameters
  result <- rstan::optimizing(
    stanmodels$lm,
    data = standata,
    init = initdata,
    hessian = TRUE
  )

  # Decide which parameters to return
  if (return_all_pars) pars <- names(result$par)
  else                 pars <- names(result$par)[1:3]

  # Calculate output
  result <- calc_confint(
    result = result,
    standata = standata,
    model = stanmodels$gmt,
    pars = pars,
    method = ci_method,
    level = ci_level,
    options = options
  )

  # Add attributes
  attr(result, "dilution_stepsize") <- dilution_stepsize

  # Return the result
  pars[pars == "sigma"] <- "sd"
  result[pars,]

}
