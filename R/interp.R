
rcensnorm <- function(mean, sd, lower_lims, upper_lims) {

  mean <- rep_len(mean, length(lower_lims))
  sd   <- rep_len(sd, length(lower_lims))

  bounded_titer <- !is.na(lower_lims) & lower_lims < upper_lims

  min_p <- pnorm(lower_lims[bounded_titer], mean[bounded_titer], sd[bounded_titer])
  max_p <- pnorm(upper_lims[bounded_titer], mean[bounded_titer], sd[bounded_titer])

  samples <- runif(sum(bounded_titer), min_p, max_p)
  output  <- lower_lims
  output[bounded_titer] <- qnorm(samples, mean[bounded_titer], sd[bounded_titer])
  output

}

gmt_interp_logtiters <- function(result, titers) {

  titerlims <- calc_titer_lims(titers, 0)
  rcensnorm(
    mean = result["mean", "estimate"],
    sd = result["mean", "estimate"],
    lower_lims = titerlims$min_titers,
    upper_lims = titerlims$max_titers
  )

}

gmt_interp_titers <- function(result, titers) {

  log_titers <- gmt_interp_logtiters(result, titers)
  output <- as.character(round(2^log_titers*10, 2))
  output[titers == "*"] <- "*"
  output

}
