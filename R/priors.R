
sample_prior <- function(
    mu_prior_mu = 0,
    mu_prior_sigma = 100,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 5
    ) {

    # Setup data
    standata <- list(
      mu_prior_mu = mu_prior_mu,
      mu_prior_sigma = mu_prior_sigma,
      sigma_prior_alpha = sigma_prior_alpha,
      sigma_prior_beta = sigma_prior_beta
    )

    # Set initial conditions
    initdata <- list(
      mu = mu_prior_mu,
      sigma = 1
    )

    # Run the sampler
    result <- rstan::sampling(
      stanmodels$priors,
      data = standata,
      init = list(as.list(initdata)),
      chains = 1,
      refresh = 0,
      iter = 2000,
      warmup = 1000
    )

    # Return the result
    rstan::extract(result)

}

# output <- sample_prior(
#   sigma_prior_alpha = 2,
#   sigma_prior_beta = 0.5
# )
# hist(output$sigma, xlim = c(0, 5), breaks = seq(from = 0, to = max(output$sigma)+0.2, by = 0.2), ylim = c(0, 400))
# abline(v = quantile(output$sigma, 0.025))
# abline(v = quantile(output$sigma, 0.975))
