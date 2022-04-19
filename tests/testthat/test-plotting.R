
rm(list = ls())
devtools::load_all()
library(Racmacs)
map <- read.acmap("~/Dropbox/labbook/publications/in review/duke-sars2-cartography-paper/data/map_ndsubset_no_outliers.ace")
titers <- t(titerTable(map)[c("B.1.351", "P.1"), ])
titers1 <- t(titerTable(map)[c("D614G"), srGroups(map) == "D614G"])
titers2 <- t(titerTable(map)[c("B.1.351"), srGroups(map) == "D614G"])
logtiters1 <- t(logtiterTable(map)[c("D614G"), srGroups(map) == "D614G"])
logtiters2 <- t(logtiterTable(map)[c("B.1.351"), srGroups(map) == "D614G"])

# hist_log2diff(
#   titers1,
#   titers2
# ) %>% plot()

hist_gmt(
  titers2,
  binwidth = 1
) %>% plot()
