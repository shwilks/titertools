
#' Calculate the mean difference between two paired sets of titers
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
#' @param sigma A fixed parameter for sigma, default is NA and a prior inverse gamma distribution distribution will be assumed based on sigma_prior_alpha and sigma_prior_beta
#' @param sigma_prior_alpha Prior for the alpha parameter of the inverse gamma distribution to estimate the standard deviation
#' @param sigma_prior_beta Prior for the beta parameter of the inverse gamma distribution to estimate the standard deviation
#' @param options Options for the sampler
#'
#' @export
log2diff <- function(
  titers1,
  titers2,
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

  # Check input
  if (!is.character(titers1) ||
      !is.character(titers2) ||
      !is.vector(titers1) ||
      !is.vector(titers2)
    ) {
    stop("Titers must be input as character vectors")
  }

  if (length(titers1) != length(titers2)) {
    stop("Titers inputs must be matching length")
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
          c("mean", "sd"),
          c("estimate", "lower", "upper")
        )
      )
    )
  }

  # Calculate titer limits
  titerdifflims <- calc_titer_diff_lims(
    titers1 = titers1,
    titers2 = titers2,
    dilution_stepsize = dilution_stepsize
  )

  # Setup data
  standata <- list(
    upper_lims = as.array(titerdifflims$max_diffs),
    lower_lims = as.array(titerdifflims$min_diffs),
    N = sum(!na_titers),
    mu_prior_mu = mu_prior_mu,
    mu_prior_sigma = mu_prior_sigma,
    sigma_prior_alpha = sigma_prior_alpha,
    sigma_prior_beta = sigma_prior_beta
  )

  # Set initial conditions
  initdata <- list(
    mu = mean(titerdifflims$logtiter_diffs)
  )

  # Optimize parameters
  if (is.na(sigma)) {

    # Where sigma is not fixed
    initdata$sigma <- pmax(sd(titerdifflims$logtiter_diffs), 0.1, na.rm = T)

    result <- tryCatch(
      rstan::optimizing(
        stanmodels$gmt,
        data = standata,
        init = initdata,
        hessian = TRUE
      ),
      error = function(e) {
        list(
          par = vapply(initdata, \(x) NA, numeric(1))
        )
      }
    )

    result <- calc_confint(
      result = result,
      standata = standata,
      model = stanmodels$gmt,
      pars = names(result$par),
      method = ci_method,
      level = ci_level,
      options = options
    )

  } else {

    # Where sigma is fixed
    standata$sigma <- sigma

    result <- tryCatch(
      rstan::optimizing(
        stanmodels$gmt_fixed_sigma,
        data = standata,
        init = initdata,
        hessian = TRUE
      ),
      error = function(e) {
        list(
          par = vapply(initdata, \(x) NA, numeric(1))
        )
      }
    )

    result <- calc_confint(
      result = result,
      standata = standata,
      model = stanmodels$gmt_fixed_sigma,
      pars = names(result$par),
      method = ci_method,
      level = ci_level,
      options = options
    )

  }

  # Add attributes
  attr(result, "dilution_stepsize") <- dilution_stepsize

  # Return the result
  result

}


#' Calculate the mean difference between two paired sets of titers,
#' where each set is a matrix of titers.
#'
#' This function is useful for example for calculating differences where you have a set of
#' pre-vaccination and post-vaccination samples and you would like to know the mean response
#' size, accounting for non-detectable values. The calculations take into account
#' antigen reactivity biases (ag_effect) and reactivity biases for individual
#' sera (sr_effect).
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
#' @param sigma_prior_ag_effect Prior for the antigen effect
#' @param sigma_prior_sr_effect Prior for the serum effect
#' @param options Options for the sampler
#'
#' @export
log2diff_me <- function(
  titers1,
  titers2,
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

  # Check input
  if (!is.matrix(titers1) || !is.matrix(titers2) || dim(titers1) != dim(titers2)) {
    "Titers must be input as a matrix of matching dimensions"
  }

  # Set variables
  ag_num <- matrix(seq_len(ncol(titers1)), nrow(titers1), ncol(titers1), byrow = T)
  sr_num <- matrix(seq_len(nrow(titers1)), nrow(titers1), ncol(titers1), byrow = F)
  na_titers <- is_na_titer(titers1) | is_na_titer(titers2)

  # Calculate titer limits
  titerdifflims <- calc_titer_diff_lims(
    titers1 = titers1,
    titers2 = titers2,
    dilution_stepsize = dilution_stepsize
  )

  # Setup data
  standata <- list(
    N                     = sum(!na_titers),
    N_ags                 = ncol(titers1),
    N_srs                 = nrow(titers1),
    upper_lims            = as.array(titerdifflims$max_diffs[!na_titers]),
    lower_lims            = as.array(titerdifflims$min_diffs[!na_titers]),
    ag                    = as.array(as.vector(ag_num))[!na_titers],
    sr                    = as.array(as.vector(sr_num))[!na_titers],
    mu_prior_mu           = mu_prior_mu,
    mu_prior_sigma        = mu_prior_sigma,
    sigma_prior_alpha     = sigma_prior_alpha,
    sigma_prior_beta      = sigma_prior_beta,
    sigma_prior_ag_effect = sigma_prior_ag_effect,
    sigma_prior_sr_effect = sigma_prior_sr_effect
  )

  # Set initial conditions
  initdata <- list(
    mu         = mean(titerdifflims$logtiter_diffs[!na_titers]),
    sigma      = pmax(sd(titerdifflims$logtiter_diffs[!na_titers]), 0.1, na.rm = T),
    ag_effects = as.array(rep(0, ncol(titers1))),
    sr_effects = as.array(rep(0, nrow(titers1)))
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


#' Calculate the mean difference between two unpaired matrices of titers
#'
#' This function is useful for example for calculating differences where you have a set of
#' pre-vaccination and post-vaccination samples and you would like to know the mean response
#' size, accounting for non-detectable values. The calculations take into account
#' antigen reactivity biases (ag_effect) and reactivity biases for individual
#' sera (sr_effect).
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
#' @param sigma_prior_ag_effect Prior for the antigen effect
#' @param sigma_prior_sr_effect Prior for the serum effect
#' @param options Options for the sampler
#'
#' @export
log2diff_unpaired_me <- function(
  titers1,
  titers2,
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

  # Check input
  if (!is.matrix(titers1) || !is.matrix(titers2) || dim(titers1) != dim(titers2)) {
    "Titers must be input as a matrix of matching dimensions"
  }

  # Set variables
  titers <- cbind(titers1, titers2)
  ag12 <- cbind(
    matrix(1, nrow(titers1), ncol(titers1)),
    matrix(2, nrow(titers2), ncol(titers2))
  )
  ag_num <- matrix(seq_len(ncol(titers)), nrow(titers), ncol(titers), byrow = T)
  sr_num <- matrix(seq_len(nrow(titers)), nrow(titers), ncol(titers), byrow = F)
  na_titers <- is_na_titer(titers)

  # Calculate titer limits
  titerlims <- calc_titer_lims(
    titers = titers,
    dilution_stepsize = dilution_stepsize
  )

  # Setup data
  standata <- list(
    N                     = sum(!na_titers),
    N_ags                 = ncol(titers),
    N_srs                 = nrow(titers),
    upper_lims            = as.array(titerlims$max_titers[!na_titers]),
    lower_lims            = as.array(titerlims$min_titers[!na_titers]),
    ag                    = as.array(as.vector(ag_num)[!na_titers]),
    sr                    = as.array(as.vector(sr_num)[!na_titers]),
    ag12                  = as.array(as.vector(ag12)[!na_titers]),
    mu_prior_mu           = mu_prior_mu,
    mu_prior_sigma        = mu_prior_sigma,
    sigma_prior_alpha     = sigma_prior_alpha,
    sigma_prior_beta      = sigma_prior_beta,
    sigma_prior_ag_effect = sigma_prior_ag_effect,
    sigma_prior_sr_effect = sigma_prior_sr_effect,
    sr_logdiffs_mu        = 0,
    sr_logdiffs_sigma     = 100
  )

  # Set initial conditions
  ag1_mu_init <- mean(titerlims$log_titers[!na_titers & ag12 == 1])
  ag2_mu_init <- mean(titerlims$log_titers[!na_titers & ag12 == 2])

  initdata <- list(
    ags1_mu     = ag1_mu_init,
    sigma       = pmax(sd(titerlims$log_titers[!na_titers]), 0.1, na.rm = T),
    ag_effects  = as.array(rep(0, ncol(titers))),
    sr_logdiffs = as.array(rep(ag2_mu_init - ag1_mu_init, nrow(titers)))
  )

  # Optimize parameters
  result <- rstan::optimizing(
    stanmodels$logdiff_unpaired_me,
    data = standata,
    init = initdata,
    hessian = TRUE
  )

  calc_confint(
    result = result,
    standata = standata,
    model = stanmodels$logdiff_unpaired_me,
    pars = names(result$par),
    method = ci_method,
    level = ci_level,
    options = options
  )

}



