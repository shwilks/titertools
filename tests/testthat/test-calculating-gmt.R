
rm(list = ls())
library(Racmacs)
# rstantools::rstan_config(); pkgbuild::compile_dll(); devtools::load_all();
map <- read.acmap("~/Dropbox/labbook/publications/in review/duke-sars2-cartography-paper/data/map_ndsubset_no_outliers.ace")
titers <- t(titerTable(map)[c("B.1.351", "P.1"), srGroups(map) == "D614G"])
titers1 <- t(titerTable(map)[c("D614G", "B.1.1.7"), srGroups(map) == "D614G"])
titers2 <- t(titerTable(map)[c("B.1.351"), srGroups(map) == "D614G"])

result <- gmt(titers1[,1], ci_method = "ETI")

# gmt_me_interp_logtiters <- function(result, titers) {
#
#   titerlims <- calc_titer_lims(titers, 0)
#   ag_num <- matrix(seq_len(ncol(titers)), nrow(titers), ncol(titers), byrow = T)
#   sr_num <- matrix(seq_len(nrow(titers)), nrow(titers), ncol(titers), byrow = F)
#   rcensnorm(
#     mean = unname(result["mean", "estimate"] +
#       result[sprintf("ag_effects[%s]", ag_num), "estimate"] +
#       result[sprintf("sr_effects[%s]", sr_num), "estimate"]),
#     sd = result["mean", "estimate"],
#     lower_lims = titerlims$min_titers,
#     upper_lims = titerlims$max_titers
#   )
#
# }
#
# hist(gmt_me_interp_logtiters(result, titers), breaks = -20:10)
# gmt_interp_titers(result, titers2)

