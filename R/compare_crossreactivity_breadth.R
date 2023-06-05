
#' Compare cross-reactivity breadth
#'
#' This function takes an input of a table of data and different serum groups
#' and calculates how fold-drops against the different antigens measured in the
#' titer table relate to each other amongst the serum groups, relating
#' response breadth of each serum group by some factor to the response breadth in
#' the reference serum group. The typical use case is to compare the response
#' breadth of different vaccine groups, where the magnitude of fold-drops for
#' each of the antigens could be assumed to follow the same general pattern but
#' be reduced by some factor in some of the serum groups.
#'
#' @param titers A named matrix of titer data, where column names are the serum
#'   names and row names are the antigen names. Censored data should be included
#'   as e.g. "<20"
#' @param sr_groups A character vector of serum group names for each of the sera
#'   in the titer table
#' @param homologous_ag A string giving the name of the homologous antigen for
#'   all the serum groups (should be only 1 antigen name from those included in
#'   the table)
#' @param reference_sr_group A string giving the name of the reference serum group
#' @param dilution_stepsize The log2 dilution stepwise to assume for the titer
#'   data, for example 1 if titer data relates to 2-fold dilutions or 0 if it is
#'   continuous data.
#'
#' @returns Returns a list of parameter estimates from the fit.
#'
#'   "slope_factors" is the factor by which fold-differences for each of the
#'   antigens vs the homologous antigen vary according to the size of the
#'   fold-differences seen in the reference_sr_group. For example a
#'   sr_group_slope estimate of 0.5 says that fold-differences tended to be
#'   about half on the log2 scale for that serum group, so a 4-fold drop for the
#'   reference_sr_group would only be 2-fold for that group.
#'
#'   "ag_folddrops" Is a list of estimates for the baseline antigen
#'   fold-differences seen across the serum groups (on the log2 scale). The
#'   fold-difference for a particular serum group would be the ag folddrop *
#'   serum group slope.
#'
#'   "logtiter_error_sigma" Is the estimate made for the standard deviation of
#'   the fit error on the log2 scale.
#'
#' @export
#'
compare_crossreactivity_breadth <- function(
    titers,
    sr_groups,
    homologous_ag,
    reference_sr_group,
    dilution_stepsize
  ) {

  # Process data
  ags <- rownames(titers)
  srs <- colnames(titers)
  sr_groups <- droplevels(sr_groups)
  sr_groups <- forcats::fct_relevel(sr_groups, reference_sr_group)
  sr_group_levels <- levels(sr_groups)

  # Split out titers according to homologous antigen
  titers_homologous_ag <- titers[rownames(titers) == homologous_ag, , drop = F]
  titers_comparison_ags <- titers[rownames(titers) != homologous_ag, , drop = F]

  # Get homologous titer lims
  homologous_logtiter_lims <- calc_titer_lims(
    titers = titers_homologous_ag,
    dilution_stepsize = dilution_stepsize
  )

  # Get titer lims
  logtiter_lims <- calc_titer_lims(
    titers = titers_comparison_ags,
    dilution_stepsize = dilution_stepsize
  )
  lower_logtiter_lims <- logtiter_lims$min_titers
  upper_logtiter_lims <- logtiter_lims$max_titers

  # Get antigen and sera matrices
  ag_matrix <- matrix(rownames(lower_logtiter_lims), nrow = nrow(lower_logtiter_lims), ncol = ncol(lower_logtiter_lims), byrow = F)
  sr_matrix <- matrix(colnames(lower_logtiter_lims), nrow = nrow(lower_logtiter_lims), ncol = ncol(lower_logtiter_lims), byrow = T)
  sr_groups_matrix <- matrix(sr_groups, nrow = nrow(lower_logtiter_lims), ncol = ncol(lower_logtiter_lims), byrow = T)

  # Remove NA cases
  na_cases <- is.na(lower_logtiter_lims) | is.na(upper_logtiter_lims)
  lower_logtiter_lims <- lower_logtiter_lims[!na_cases]
  upper_logtiter_lims <- upper_logtiter_lims[!na_cases]
  ag_matrix           <- ag_matrix[!na_cases]
  sr_matrix           <- sr_matrix[!na_cases]
  sr_groups_matrix    <- sr_groups_matrix[!na_cases]

  # Get sera levels
  ag_levels <- rownames(titers_comparison_ags)
  sr_levels <- colnames(titers_comparison_ags)

  # Factor antigen and sera cases
  ag_matrix <- factor(ag_matrix, ag_levels)
  sr_matrix <- factor(sr_matrix, sr_levels)
  sr_groups_matrix <- factor(sr_groups_matrix, sr_group_levels)

  # Assemble data list
  standata <- list(
    N = length(lower_logtiter_lims),
    N_ags = nlevels(ag_matrix),
    N_srs = nlevels(sr_matrix),
    N_sr_groups = nlevels(sr_groups_matrix),
    upper_logtiter_lims = upper_logtiter_lims,
    lower_logtiter_lims = lower_logtiter_lims,
    upper_homologous_logtiter_lims = as.vector(homologous_logtiter_lims$max_titers),
    lower_homologous_logtiter_lims = as.vector(homologous_logtiter_lims$min_titers),
    ags = as.numeric(ag_matrix),
    srs = as.numeric(sr_matrix),
    sr_groups = as.numeric(sr_groups_matrix)
  )

  initdata <- list(
    sr_group_slopes = array(rep(0.9, nlevels(sr_groups) - 1)),
    ag_folddrops = array(rep(0, nlevels(ag_matrix))),
    sr_homologous_logtiters = as.vector(homologous_logtiter_lims$log_titers),
    logtiter_error_sigma = 1
  )

  # Optimize model
  result <- rstan::optimizing(
    stanmodels$compare_crossreactivity_breadth,
    data = standata,
    init = initdata,
    hessian = TRUE
  )

  # Calculate confidence intervals
  output <- calc_confint(
    result = result,
    standata = standata,
    model = stanmodels$compare_crossreactivity_breadth,
    pars = names(result$par),
    method = "HDI",
    level = 0.95,
    options = list(
      iter = 5000,
      warmup = 1000
    )
  )

  # Reformat output
  slope_results <- output[grepl("sr_group_slopes", rownames(output)), , drop = F] |>
    tibble::as_tibble(
      rownames = "variable"
    ) |>
    dplyr::mutate(
      sr_group = sr_group_levels[-1]
    )

  ag_folddrop_results <- output[grepl("ag_folddrops", rownames(output)), , drop = F] |>
    tibble::as_tibble(
      rownames = "variable"
    ) |>
    dplyr::mutate(
      ag = levels(ag_matrix)
    )

  logtiter_error_sigma_results <- output[grepl("logtiter_error_sigma", rownames(output)), , drop = F]

  # Return the output
  list(
    slope_factors = slope_results,
    ag_folddrops = ag_folddrop_results,
    logtiter_error_sigma = logtiter_error_sigma_results
  )

}
