
#' Calculate the geometric mean titer of a set of titers
#'
#' @param titers A vector of titers
#' @param ci_method The method to use when calculating the confidence intervals
#' @param ci_level The confidence level to use when calculating confidence intervals
#' @param options Options for the sampler
#'
#' @export
gmt <- function(
  titers,
  ci_method = "quap",
  ci_level = 0.95,
  options = list()
  ) {

  # Remove na titers
  titers <- titers[!is_na_titer(titers)]

  # Calculate titer limits
  titerlims <- calc_titer_lims(titers, 0)

  # Setup data
  standata <- list(
    upper_lims = titerlims$max_titers,
    lower_lims = titerlims$min_titers,
    N = length(titers),
    mu_prior_mu = 0,
    mu_prior_sigma = 100,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 5
  )

  # Set initial conditions
  initdata <- list(
    mu = mean(titerlims$log_titers),
    sigma = sd(titerlims$log_titers)
  )

  # Optimize parameters
  result <- rstan::optimizing(
    stanmodels$gmt,
    data = standata,
    init = initdata,
    hessian = TRUE
  )

  # Save output
  output <- list(
    mean = result$par["mu"],
    sd = result$par["sigma"],
    mean_lower = NA,
    mean_upper = NA,
    sd_lower = NA,
    sd_upper = NA
  )

  calc_confint(
    result = result,
    model = stanmodels$gmt,
    pars = names(result$par),
    method = ci_method,
    level = ci_level,
    options = options
  )

}


#' Calculate the geometric mean titer of a set of titers
#'
#' @param titers A matrix of titers
#' @param ci_method The method to use when calculating the confidence intervals
#' @param ci_level The confidence level to use when calculating confidence intervals
#' @param options Options for the sampler
#'
#' @export
gmt_me <- function(
  titers,
  ci_method = "quap",
  ci_level = 0.95,
  options = list()
) {

  ag_num <- matrix(seq_len(ncol(titers)), nrow(titers), ncol(titers), byrow = T)
  sr_num <- matrix(seq_len(nrow(titers)), nrow(titers), ncol(titers), byrow = F)
  na_titers <- is_na_titer(titers)
  titerlims <- calc_titer_lims(titers, 0)

  # Setup data
  standata <- list(
    N                     = sum(!na_titers),
    N_ags                 = ncol(titers),
    N_srs                 = nrow(titers),
    upper_lims            = titerlims$max_titers[!na_titers],
    lower_lims            = titerlims$min_titers[!na_titers],
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
    mu         = mean(titerlims$log_titers[!na_titers]),
    sigma      = sd(titerlims$log_titers[!na_titers]),
    ag_effects = rep(0, ncol(titers)),
    sr_effects = rep(0, nrow(titers))
  )

  # Optimize parameters
  browser()
  result <- rstan::optimizing(
    stanmodels$gmt_me,
    data = standata,
    init = initdata,
    hessian = TRUE
  )

  calc_confint(
    result = result,
    model = stanmodels$gmt_me,
    pars = names(result$par),
    method = ci_method,
    level = ci_level,
    options = options
  )

}


