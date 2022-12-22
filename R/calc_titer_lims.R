
is_na_titer <- function(titers) {

  titers %in% c("*", ".") | is.na(titers)

}



titer_fit_options <- function(
  min_titer_possible = -Inf,
  max_titer_possible = Inf
) {

  list(
    min_titer_possible = min_titer_possible,
    max_titer_possible = max_titer_possible
  )

}



#' Get titer limits
#'
#' Function for getting upper and lower limits of measured titers on the log scale.
#'
#' @param titers A numeric/character vector or matrix of titer measurements.
#' @param dilution_stepsize The dilution stepsize, see details
#' @param options A named list of options to pass to `titer_fit_options()`
#'
#' @return Returns a list of length three with values log_titers, max_titers and min_titers, giving the
#' numeric vectors of the straight log converted, upper and lower bounds of the titers on the log scale.
#'
#' @details HI titers which typically follow a 2-fold dilution series starting
#'   from 1/10, then 1/20, 1/40 etc. This represents a "dilution stepsize" of 1
#'   when converted to the log2 scale. When no inhibition was recorded at the
#'   highest dilution, the value is typically recorded as <10 but the
#'   optimization regime effectively treats this as a <=5, the rationale being
#'   that, had the dilution series been continued to higher concentrations, the
#'   next lowest titer would have been a 5. Over time the method has also been
#'   applied to other neutralization assays that sometimes have a continuous
#'   read out with a lower end, in these cases a <10 really means a <10 since
#'   any other values like 9.8 or 7.62 would also be possible. To indicate these
#'   continuous cases, you can specify the dilution stepsize as 0. Equally, if
#'   the dilution regime followed a different pattern, you can also set that
#'   here.
#'
#' @examples
#' # Calculate the titer limits of a set of HI titers
#' titer_lims <- calc_titer_lims(
#'   titers = c("20", "320", "<10", ">1280"),
#'   dilution_stepsize = 1
#' )
#'
#' # Calculate the titer limits assuming non-default upper and lower bounds for non-detectable
#' # and greater-than titers.
#' titer_lims <- calc_titer_lims(
#'   titers = c("20", "320", "<10", ">1280"),
#'   dilution_stepsize = 1,
#'   options = list(
#'     min_titer_possible = -Inf,
#'     max_titer_possible = 14
#'   )
#' )
#'
#' @noRd
#'
calc_titer_lims <- function(
  titers,
  dilution_stepsize,
  options = list()
) {

  # Get options
  options <- do.call(titer_fit_options, options)

  # Find less than and greater than titers and convert them to a numeric form
  lessthan_titers <- grepl(x = titers, pattern = "<")
  morethan_titers <- grepl(x = titers, pattern = ">")
  na_titers       <- is.na(titers) | titers == "*" | titers == "."

  numeric_titers <- titers
  numeric_titers[na_titers] <- NA
  numeric_titers <- gsub("(<|>)", "", numeric_titers)
  mode(numeric_titers) <- "numeric"

  # Convert titers to the log scale
  log_titers <- log2(numeric_titers / 10)
  log_titers[lessthan_titers] <- log_titers[lessthan_titers] - dilution_stepsize
  log_titers[morethan_titers] <- log_titers[morethan_titers] + dilution_stepsize
  max_titers <- log_titers + dilution_stepsize / 2
  min_titers <- log_titers - dilution_stepsize / 2
  min_titers[lessthan_titers] <- options$min_titer_possible
  max_titers[morethan_titers] <- options$max_titer_possible

  # Create a vector indicating censoring
  censoring <- rep("none", length(titers))
  if (dilution_stepsize > 0) censoring[] <- "interval"
  censoring[lessthan_titers] <- "left"
  censoring[morethan_titers] <- "right"
  censoring <- factor(censoring, levels = c("none", "left", "right", "interval", "full"))

  list(
    log_titers = log_titers,
    max_titers = max_titers,
    min_titers = min_titers,
    censoring  = censoring
  )

}


calc_titer_diff_lims <- function(
  titers1,
  titers2,
  dilution_stepsize,
  options = list()
) {

  # Get limits to titer measurements
  titer1_lims <- calc_titer_lims(titers1, dilution_stepsize, options)
  titer2_lims <- calc_titer_lims(titers2, dilution_stepsize, options)

  # Compute maximum and mininmum differences
  max_diffs <- titer2_lims$max_titers - titer1_lims$min_titers
  min_diffs <- titer2_lims$min_titers - titer1_lims$max_titers

  logtiters1 <- titer1_lims$log_titers
  logtiters2 <- titer2_lims$log_titers

  # Determine censoring
  censoring <- rep("none", length(titers1))
  censoring[titer1_lims$censoring == "none"     & titer2_lims$censoring == "none"]     <- "none"
  censoring[titer1_lims$censoring == "none"     & titer2_lims$censoring == "interval"] <- "interval"
  censoring[titer1_lims$censoring == "none"     & titer2_lims$censoring == "right"]    <- "right"
  censoring[titer1_lims$censoring == "none"     & titer2_lims$censoring == "left"]     <- "left"
  censoring[titer1_lims$censoring == "interval" & titer2_lims$censoring == "none"]     <- "none"
  censoring[titer1_lims$censoring == "interval" & titer2_lims$censoring == "interval"] <- "interval"
  censoring[titer1_lims$censoring == "interval" & titer2_lims$censoring == "right"]    <- "right"
  censoring[titer1_lims$censoring == "interval" & titer2_lims$censoring == "left"]     <- "left"
  censoring[titer1_lims$censoring == "left"     & titer2_lims$censoring == "none"]     <- "right"
  censoring[titer1_lims$censoring == "left"     & titer2_lims$censoring == "interval"] <- "right"
  censoring[titer1_lims$censoring == "left"     & titer2_lims$censoring == "left"]     <- "full"
  censoring[titer1_lims$censoring == "left"     & titer2_lims$censoring == "right"]    <- "right"
  censoring[titer1_lims$censoring == "right"    & titer2_lims$censoring == "none"]     <- "left"
  censoring[titer1_lims$censoring == "right"    & titer2_lims$censoring == "interval"] <- "left"
  censoring[titer1_lims$censoring == "right"    & titer2_lims$censoring == "right"]    <- "full"
  censoring[titer1_lims$censoring == "right"    & titer2_lims$censoring == "left"]     <- "left"
  censoring <- factor(censoring, levels = c("none", "left", "right", "interval", "full"))

  # Return result
  list(
    max_diffs = max_diffs,
    min_diffs = min_diffs,
    logtiter_diffs = logtiters2 - logtiters1,
    censoring = censoring
  )

}


sampler_options <- function(
  iter = 10000,
  warmup = 2000
) {

  list(
    iter = iter,
    warmup = warmup
  )

}


infer_dilution_stepsize <- function(titers) {

  warning("Dilution stepsize inferred as 0, to suppress this message set dilution stepsize explicitly using the 'dilution_stepsize' argument")
  return(0)

}
