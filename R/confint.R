
#' Calculate the confidence intervals
#'
#' @param result The result of the optimised model.
#' @param standata The data passed to stan.
#' @param model The stan model.
#' @param pars The parameters that were optimised.
#' @param method The method for calculating the confidence interval.
#' One of 'quap', 'ETI', 'HDI', 'BCI', 'SI', see details.
#' @param level The confidence level to use when calculating confidence intervals
#' @param options Options for the sampler
#'
#' @details
#' `quap`: Confidence interval is computed using quadratic approximation.
#' For ETI, HDI, BCI and SI, the confidence interval is computed across samples from the posterior distribution.
#' `ETI`: Equal-Tailed Interval. The probability of being below this interval is equal to the probability of being above it. A 95% ETI will always have 2.5% of the distribution on either side of its limits. See https://easystats.github.io/bayestestR/reference/eti.html
#' `HDI`: Highest Density Interval. Not equal-tailed. See https://easystats.github.io/bayestestR/reference/hdi.html
#' `BCI`: Bias Corrected and Accelerated Interval. See https://easystats.github.io/bayestestR/reference/bci.html
#' `SI`: Support Interval. 'A support interval contains only the values of the parameter that predict the observed data better than average'. See https://easystats.github.io/bayestestR/reference/si.html
#'
#' @noRd
#'
calc_confint <- function(
  result,
  standata,
  model,
  pars,
  method = "HDI",
  level = 0.95,
  options = list()
) {

  # Setup output
  output <- cbind(result$par, NA, NA)
  colnames(output) <- c("estimate", "lower", "upper")

  # Set options
  options <- do.call(sampler_options, options)

  # Calculate confidence intervals
  if (method == "quap") {

    covariance_matrix <- solve(result$hessian)
    alpha <- 1 - level
    par_sd <- -diag(covariance_matrix)
    names(par_sd) <- names(result$par)

    ci_output <- do.call(
      rbind,
      lapply(pars, function(par) {
        result$par[par] + qnorm(c(0.5, alpha/2, 1 - alpha/2))*par_sd[par]
      })
    )

  } else if (method %in% c("ETI", "HDI", "BCI", "SI")) {

    result_sampling <- rstan::sampling(
      model,
      data = standata,
      # init = list(as.list(result$par)),
      # chains = 1,
      refresh = 0,
      iter = options$iter,
      warmup = options$warmup
    )

    ci_output <- do.call(
      rbind,
      lapply(pars, function(par) {
        par_sample <- rstan::extract(result_sampling, pars = par)[[1]]
        par_ci <- bayestestR::ci(par_sample, ci = level, method = method)
        c(result$par[par], par_ci$CI_low, par_ci$CI_high)
      })
    )

  } else {

    stop("method must be one of 'quap', 'ETI', 'HDI', 'BCI' or 'SI'")

  }

  # Return output
  output[pars,] <- ci_output
  rownames(output)[rownames(output) == "mu"] <- "mean"
  rownames(output)[rownames(output) == "sigma"] <- "sd"
  output

}
