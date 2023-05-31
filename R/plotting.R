
# hist_gmt <- function(
#   titers,
#   binwidth = NULL,
#   ...
# ) {
#
#   # Calculate the log2diff
#   gmt_est <- gmt(
#     titers = titers,
#     ...
#   )
#
#   # Impute the log2diff
#   imputed_logtiters <- impute_gmt_logtiters(
#     result = gmt_est,
#     titers = titers
#   )
#
#   # Plot the histogram
#   tibble::tibble(
#     logtiters = as.vector(imputed_logtiters),
#     censoring = attr(imputed_logtiters, "censoring")
#   ) %>%
#     dplyr::filter(
#       !is.na(logtiters)
#     ) %>%
#     ggplot2::ggplot(
#       ggplot2::aes(
#         x = logtiters,
#         fill = censoring
#       )
#     ) +
#     ggplot2::geom_histogram(
#       binwidth = binwidth
#     )
#
# }
#
#
# hist_log2diff <- function(
#   titers1,
#   titers2,
#   binwidth = NULL,
#   ...
#   ) {
#
#   # Calculate the log2diff
#   log2diff_est <- log2diff(
#     titers1 = titers1,
#     titers2 = titers2,
#     ...
#   )
#
#   # Impute the log2diff
#   imputed_log2diff <- impute_log2diff(
#     result = log2diff_est,
#     titers1 = titers1,
#     titers2 = titers2
#   )
#
#   # Plot the histogram
#   tibble::tibble(
#     log2diff = as.vector(imputed_log2diff),
#     censoring = attr(imputed_log2diff, "censoring")
#   ) %>%
#     dplyr::filter(
#       !is.na(log2diff)
#     ) %>%
#     ggplot2::ggplot(
#       ggplot2::aes(
#         x = log2diff,
#         fill = censoring
#       )
#     ) +
#     ggplot2::geom_histogram(
#       binwidth = binwidth
#     )
#
# }
