
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
      init = list(as.list(result$par)),
      chains = 1,
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

    stop("ci_type must be one of 'quap', 'ETI', 'HDI', 'BCI' or 'SI'")

  }

  # Return output
  output[pars,] <- ci_output
  rownames(output)[rownames(output) == "mu"] <- "mean"
  rownames(output)[rownames(output) == "sigma"] <- "sd"
  output

}
