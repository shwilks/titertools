
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

impute_gmt_logtiters <- function(result, titers) {

  titerlims <- calc_titer_lims(titers, 0)
  result <- rcensnorm(
    mean = result["mean", "estimate"],
    sd = result["sd", "estimate"],
    lower_lims = titerlims$min_titers,
    upper_lims = titerlims$max_titers
  )

  attr(result, "censoring") <- titerlims$censoring
  result

}

impute_gmt_titers <- function(result, titers) {

  log_titers <- impute_gmt_logtiters(result, titers)
  output <- as.character(round(2^log_titers*10, 2))
  output[titers == "*"] <- "*"
  output

}

impute_log2diff <- function(
  result,
  titers1,
  titers2,
  options = list()
  ) {

  titerlims <- calc_titer_diff_lims(
    titers1 = titers1,
    titers2 = titers2,
    dilution_stepsize = attr(result, "dilution_stepsize"),
    options = options
  )

  result <- rcensnorm(
    mean = result["mean", "estimate"],
    sd = result["sd", "estimate"],
    lower_lims = titerlims$min_diffs,
    upper_lims = titerlims$max_diffs
  )

  attr(result, "censoring") <- titerlims$censoring
  result

}
