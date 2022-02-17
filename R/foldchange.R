
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
#' @param options Options for the sampler
#'
#' @export
foldchange <- function(
  titers1,
  titers2,
  ci_method = "quap",
  ci_level = 0.95,
  options = list()
  ) {

  # Check input
  if (!is.vector(titers1) || !is.vector(titers2) || length(titers1) != length(titers2)) {
    "Titers must be input as vectors of matching length"
  }

  # Remove na titers
  na_titers <- is_na_titer(titers1) | is_na_titer(titers2)
  titers1 <- titers1[!na_titers]
  titers2 <- titers2[!na_titers]

  # Calculate titer limits
  titerdifflims <- calc_titer_diff_lims(
    titers1 = titers1,
    titers2 = titers2,
    dilution_stepsize = 0
  )

  # Setup data
  standata <- list(
    upper_lims = titerdifflims$max_diffs,
    lower_lims = titerdifflims$min_diffs,
    N = sum(!na_titers),
    mu_prior_mu = 0,
    mu_prior_sigma = 100,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 5
  )

  # Set initial conditions
  initdata <- list(
    mu = mean(titerdifflims$logtiter_diffs),
    sigma = pmax(sd(titerdifflims$logtiter_diffs), 0.1)
  )

  # Optimize parameters
  result <- rstan::optimizing(
    stanmodels$gmt,
    data = standata,
    init = initdata,
    hessian = TRUE
  )

  # Calculate output
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
#' @param options Options for the sampler
#'
#' @export
foldchange_me <- function(
  titers1,
  titers2,
  ci_method = "quap",
  ci_level = 0.95,
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
    dilution_stepsize = 0
  )

  # Setup data
  standata <- list(
    N                     = sum(!na_titers),
    N_ags                 = ncol(titers1),
    N_srs                 = nrow(titers1),
    upper_lims            = titerdifflims$max_diffs[!na_titers],
    lower_lims            = titerdifflims$min_diffs[!na_titers],
    ag                    = as.vector(ag_num)[!na_titers],
    sr                    = as.vector(sr_num)[!na_titers],
    mu_prior_mu           = 0,
    mu_prior_sigma        = 100,
    sigma_prior_alpha     = 2,
    sigma_prior_beta      = 5,
    sigma_prior_ag_effect = 0.7,
    sigma_prior_sr_effect = 2
  )

  # Set initial conditions
  initdata <- list(
    mu         = mean(titerdifflims$logtiter_diffs[!na_titers]),
    sigma      = sd(titerdifflims$logtiter_diffs[!na_titers]),
    ag_effects = rep(0, ncol(titers1)),
    sr_effects = rep(0, nrow(titers1))
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

