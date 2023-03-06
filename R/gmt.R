
#' Calculate the geometric mean titer of a set of titers.
#'
#' @param titers A vector of titers
#' @param ci_method The method to use when calculating the confidence intervals
#' @param ci_level The confidence level to use when calculating confidence intervals
#' @param dilution_stepsize The dilution stepsize. Generally 1 for titers of
#' the form '20', '40', '80', and 0 for titers that are read out continuously.
#' @param mu_prior_mu Prior for the mean of the normal distribution to estimate the mean
#' @param mu_prior_sigma Prior for the standard deviation of the normal distribution to estimate the mean
#' @param sigma A fixed parameter for sigma, default is NA and a prior inverse gamma distribution distribution will be assumed based on sigma_prior_alpha and sigma_prior_beta
#' @param sigma_prior_alpha Prior for the alpha parameter of the inverse gamma distribution to estimate the standard deviation
#' @param sigma_prior_beta Prior for the beta parameter of the inverse gamma distribution to estimate the standard deviation
#' @param options Options for the sampler
#'
#' @export
gmt <- function(
  titers,
  ci_method = "HDI",
  ci_level = 0.95,
  dilution_stepsize = NA,
  mu_prior_mu = 0,
  mu_prior_sigma = 100,
  sigma = NA,
  sigma_prior_alpha = 2,
  sigma_prior_beta = 5,
  options = list()
  ) {

  # Remove na titers
  titers <- titers[!is_na_titer(titers)]
  if (length(titers) == 0) {
    return(
      matrix(
        NA_real_, 2, 3,
        dimnames = list(
          c("mean", "sd"),
          c("estimate", "lower", "upper")
        )
      )
    )
  }

  # Infer dilution stepsize if necessary
  if (is.na(dilution_stepsize)) dilution_stepsize <- infer_dilution_stepsize(titers)

  # Calculate titer limits
  titerlims <- calc_titer_lims(titers, dilution_stepsize)

  # Setup data
  standata <- list(
    upper_lims = as.array(titerlims$max_titers),
    lower_lims = as.array(titerlims$min_titers),
    N = length(titers),
    mu_prior_mu = mu_prior_mu,
    mu_prior_sigma = mu_prior_sigma,
    sigma_prior_alpha = sigma_prior_alpha,
    sigma_prior_beta = sigma_prior_beta
  )

  # Set initial conditions
  initdata <- list(
    mu = mean(titerlims$log_titers)
  )

  # Optimize parameters
  if (is.na(sigma)) {

    # Where sigma is not fixed
    initdata$sigma <- pmax(sd(titerlims$log_titers), 0.1, na.rm = T)

    result <- rstan::optimizing(
      stanmodels$gmt,
      data = standata,
      init = initdata,
      hessian = TRUE
    )

  } else {

    # Where sigma is fixed
    standata$sigma <- sigma

    result <- rstan::optimizing(
      stanmodels$gmt_fixed_sigma,
      data = standata,
      init = initdata,
      hessian = TRUE
    )

  }

  calc_confint(
    result = result,
    standata = standata,
    model = stanmodels$gmt,
    pars = names(result$par),
    method = ci_method,
    level = ci_level,
    options = options
  )

}


#' Calculate the geometric mean titer of a titertable, taking into account
#' antigen reactivity biases (ag_effect) and reactivity biases for individual
#' sera (sr_effect).
#'
#' @param titers A matrix of titers
#' @param ci_method The method to use when calculating the confidence intervals
#' @param ci_level The confidence level to use when calculating confidence intervals
#' @param dilution_stepsize The dilution stepsize. Generally 1 for titers of
#' the form '20', '40', '80', and 0 for titers that are read out continuously.
#' @param mu_prior_mu Prior for the mean of the normal distribution to estimate the mean
#' @param mu_prior_sigma Prior for the standard deviation of the normal distribution to estimate the mean
#' @param sigma_prior_alpha Prior for the alpha parameter of the inverse gamma distribution to estimate the standard deviation
#' @param sigma_prior_beta Prior for the beta parameter of the inverse gamma distribution to estimate the standard deviation
#' @param sigma_prior_ag_effect Prior for the antigen effect
#' @param sigma_prior_sr_effect Prior for the serum effect
#' @param options Options for the sampler
#'
#' @export
gmt_me <- function(
  titers,
  ci_method = "HDI",
  ci_level = 0.95,
  dilution_stepsize = 0,
  mu_prior_mu = 0,
  mu_prior_sigma = 100,
  sigma_prior_alpha = 2,
  sigma_prior_beta = 5,
  sigma_prior_ag_effect = 0.7,
  sigma_prior_sr_effect = 2,
  options = list()
) {

  ag_num <- matrix(seq_len(ncol(titers)), nrow(titers), ncol(titers), byrow = T)
  sr_num <- matrix(seq_len(nrow(titers)), nrow(titers), ncol(titers), byrow = F)
  na_titers <- is_na_titer(titers)
  titerlims <- calc_titer_lims(titers, dilution_stepsize)

  if (sum(!na_titers) == 0) {
    return(
      matrix(
        NA_real_, 2 + ncol(titers) + nrow(titers), 3,
        dimnames = list(
          c("mean", "sd", sprintf("ag_effects[%s]", seq_len(ncol(titers))), sprintf("sr_effects[%s]", seq_len(nrow(titers)))),
          c("estimate", "lower", "upper")
        )
      )
    )
  }

  # Setup data
  standata <- list(
    N                     = sum(!na_titers),
    N_ags                 = ncol(titers),
    N_srs                 = nrow(titers),
    upper_lims            = as.array(titerlims$max_titers[!na_titers]),
    lower_lims            = as.array(titerlims$min_titers[!na_titers]),
    ag                    = as.array(as.vector(ag_num)[!na_titers]),
    sr                    = as.array(as.vector(sr_num)[!na_titers]),
    mu_prior_mu           = mu_prior_mu,
    mu_prior_sigma        = mu_prior_sigma,
    sigma_prior_alpha     = sigma_prior_alpha,
    sigma_prior_beta      = sigma_prior_beta,
    sigma_prior_ag_effect = sigma_prior_ag_effect,
    sigma_prior_sr_effect = sigma_prior_sr_effect
  )

  # Set initial conditions
  initdata <- list(
    mu         = mean(titerlims$log_titers[!na_titers]),
    sigma      = pmax(sd(titerlims$log_titers[!na_titers]), 0.1, na.rm = T),
    ag_effects = as.array(rep(0, ncol(titers))),
    sr_effects = as.array(rep(0, nrow(titers)))
  )

  # Optimize parameters
  result <- rstan::optimizing(
    stanmodels$gmt_me,
    data = standata,
    init = initdata,
    hessian = TRUE
  )

  calc_confint(
    result = result,
    standata = standata,
    model = stanmodels$gmt_me,
    pars = names(result$par),
    method = ci_method,
    level = ci_level,
    options = options
  )

}


#' Estimate sera effects
#'
#' @param titers A matrix of titers
#' @param ci_method The method to use when calculating the confidence intervals
#' @param ci_level The confidence level to use when calculating confidence intervals
#' @param dilution_stepsize The dilution stepsize. Generally 1 for titers of
#' the form '20', '40', '80', and 0 for titers that are read out continuously.
#' @param ag_mean_prior_mean Prior for the mean parameter of the antigen mean values
#' @param ag_mean_prior_sigma Prior for the standard deviation parameter of the antigen mean values
#' @param sigma A fixed parameter for sigma, default is NA and a prior inverse gamma distribution distribution will be assumed based on sigma_prior_alpha and sigma_prior_beta
#' @param sigma_prior_alpha Prior for the alpha parameter of the inverse gamma distribution to estimate the standard deviation
#' @param sigma_prior_beta Prior for the beta parameter of the inverse gamma distribution to estimate the standard deviation
#' @param sigma_prior_sr_effect Prior for the serum effect
#' @param options Options for the sampler
#'
#' @export
estimate_sr_effects <- function(
    titers,
    ci_method = "HDI",
    ci_level = 0.95,
    dilution_stepsize = 0,
    ag_mean_prior_mean = 0,
    ag_mean_prior_sigma = 1000,
    sigma = NA,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 5,
    sigma_prior_sr_effect = 4,
    options = list()
) {

  # Determine NA ags and sera
  na_ags <- rowSums(titers != "*") == 0
  na_srs <- colSums(titers != "*") == 0

  # Remove NA ags and seras
  titers <- titers[!na_ags, !na_srs, drop = FALSE]

  num_ags <- nrow(titers)
  num_srs <- ncol(titers)

  ag_num <- matrix(seq_len(num_ags), nrow(titers), ncol(titers), byrow = F)
  sr_num <- matrix(seq_len(num_srs), nrow(titers), ncol(titers), byrow = T)
  na_titers <- is_na_titer(titers)
  titerlims <- calc_titer_lims(titers, dilution_stepsize)

  if (sum(!na_titers) == 0) {
    return(
      matrix(
        NA_real_, 1 + num_ags + num_srs, 3,
        dimnames = list(
          c("sd", sprintf("ag_means[%s]", seq_len(num_ags)), sprintf("sr_effects[%s]", seq_len(num_srs))),
          c("estimate", "lower", "upper")
        )
      )
    )
  }


  # Setup data
  standata <- list(
    N                     = sum(!na_titers),
    N_ags                 = num_ags,
    N_srs                 = num_srs,
    upper_lims            = as.array(titerlims$max_titers[!na_titers]),
    lower_lims            = as.array(titerlims$min_titers[!na_titers]),
    ag                    = as.array(as.vector(ag_num)[!na_titers]),
    sr                    = as.array(as.vector(sr_num)[!na_titers]),
    ag_mean_prior_mean    = ag_mean_prior_mean,
    ag_mean_prior_sigma   = ag_mean_prior_sigma,
    sigma_prior_alpha     = sigma_prior_alpha,
    sigma_prior_beta      = sigma_prior_beta,
    sigma_prior_sr_effect = sigma_prior_sr_effect
  )

  # Set initial conditions
  initdata <- list(
    sr_effects = as.array(rep(0, num_srs))
  )

  # Optimize parameters
  if (is.na(sigma)) {

    # Sigma not fixed
    initdata$sigma <- pmax(sd(titerlims$log_titers[!na_titers]), 0.1, na.rm = T)

    result <- rstan::optimizing(
      stanmodels$sr_effects,
      data = standata,
      init = initdata,
      hessian = TRUE
    )

  } else {

    # Sigma fixed
    standata$sigma <- sigma

    result <- rstan::optimizing(
      stanmodels$sr_effects_fixed_sigma,
      data = standata,
      init = initdata,
      hessian = TRUE
    )

  }

  # Fetch results
  results <- calc_confint(
    result = result,
    standata = standata,
    model = stanmodels$sr_effects,
    pars = names(result$par),
    method = ci_method,
    level = ci_level,
    options = options
  )

  # Replace NA ags and sera
  ag_means_results_subset  <- results[sprintf("ag_means[%s]", seq_len(num_ags)), , drop = FALSE]
  sr_effect_results_subset <- results[sprintf("sr_effects[%s]", seq_len(num_srs)), , drop = FALSE]

  ag_means_results <- matrix(
    data = NA,
    nrow = length(na_ags),
    ncol = 3,
    dimnames = list(
      sprintf("ag_means[%s]", seq_along(na_ags)),
      colnames(ag_means_results_subset)
    )
  )

  sr_effect_results <- matrix(
    data = NA,
    nrow = length(na_srs),
    ncol = 3,
    dimnames = list(
      sprintf("sr_effects[%s]", seq_along(na_srs)),
      colnames(sr_effect_results_subset)
    )
  )

  ag_means_results[!na_ags,] <- ag_means_results_subset
  sr_effect_results[!na_srs,] <- sr_effect_results_subset

  # Return the results
  output <- rbind(
    ag_means_results,
    sr_effect_results
  )

  if (is.na(sigma)) {
    output <- rbind(
      results["sd", , drop = FALSE],
      output
    )
  }

  output

}


