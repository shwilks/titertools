
rm(list = ls())
library(Racmacs)
# rstantools::rstan_config(); pkgbuild::compile_dll(); devtools::load_all();
map <- read.acmap("~/Dropbox/labbook/publications/in review/duke-sars2-cartography-paper/data/map_ndsubset_no_outliers.ace")
titers <- t(titerTable(map)[c("B.1.351", "P.1"), srGroups(map) == "D614G"])
titers1 <- t(titerTable(map)[c("D614G", "B.1.1.7"), srGroups(map) == "D614G", drop = F])
titers2 <- t(titerTable(map)[c("B.1.351", "P.1"), srGroups(map) == "D614G", drop = F])

log2diff("*", "*", ci_method = "HDI")
log2diff("<10", "<10", ci_method = "HDI")
log2diff("10", "*", ci_method = "HDI")

gmt("*")
gmt(c("40", "*"))
gmt(c("<10", "<10", "<10"))

titers1 <- c("562.47", "22.02", "<20", "47.23", "21.44", "133.97", "43.64", "<20")
titers2 <- c("629.29", "<20", "30.58", "28.34", "27.78", "<20", "<20", "21.84")

titers1 <- c("<20", "<20")
titers2 <- c("<20", "<20")

log2diff(titers1, titers2, dilution_stepsize = 0, ci_method = "HDI")
stop()

model_sample <- rstan::sampling(
  stanmodels$gmt,
  data = standata,
  init = list(as.list(staninit)),
  chains = 1,
  refresh = 0,
  iter = 10000,
  warmup = 2000
)


titers_na <- titers1
titers_na[] <- "*"
gmt_me(titers_na)

# result <- gmt(titers1[,1], ci_method = "ETI")
result <- log2diff_unpaired_me(titers1, titers2, ci_method = "ETI")

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

