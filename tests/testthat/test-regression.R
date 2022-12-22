
# Set seed
# set.seed(100)

# Set variables
ntiters <- 1000
dilution_stepsize <- 0

titerintercept <- -4
titerslope <- 1

titer1error_sd <- 0.2
titer2error_sd <- 0.2

titer1threshold <- -20
titer2threshold <- -20

# Generate titers
# titers1 <- runif(ntiters, min = -3, max = 12)
titers1 <- rnorm(ntiters, mean = 6, sd = 4)
titers2 <- titers1*titerslope + titerintercept

# Add titer noise
measuredtiters1 <- titers1 + rnorm(length(titers1), 0, titer1error_sd)
measuredtiters2 <- titers2 + rnorm(length(titers2), 0, titer2error_sd)

# Censor titers
if (dilution_stepsize == 1) {
  measuredtiters1 <- round(measuredtiters1)
  measuredtiters2 <- round(measuredtiters2)
}

censored1 <- measuredtiters1 <= titer1threshold
censored2 <- measuredtiters2 <= titer2threshold

measuredtiters1[censored1] <- titer1threshold
measuredtiters2[censored2] <- titer2threshold

# Create raw titers
rawmeasuredtiters1 <- as.character(2^measuredtiters1*10)
rawmeasuredtiters2 <- as.character(2^measuredtiters2*10)
rawmeasuredtiters1[censored1] <- paste0("<", 2^(titer1threshold+1)*10)
rawmeasuredtiters2[censored2] <- paste0("<", 2^(titer2threshold+1)*10)

# Do a linear regression
lm_result <- lm(measuredtiters2 ~ measuredtiters1)

# Calculate the fold difference
log2diff_result <- log2diff(
  titers1 = rawmeasuredtiters1,
  titers2 = rawmeasuredtiters2,
  ci_method = "quap",
  dilution_stepsize = 1
)

# # Totally censored
# totcensored <- measuredtiters1 == -1 & measuredtiters2 == -1
# titers1 <- titers1[!totcensored]
# titers2 <- titers2[!totcensored]
# measuredtiters1 <- measuredtiters1[!totcensored]
# measuredtiters2 <- measuredtiters2[!totcensored]
# rawmeasuredtiters1 <- rawmeasuredtiters1[!totcensored]
# rawmeasuredtiters2 <- rawmeasuredtiters2[!totcensored]

# Do an orthogonal regression
titerlm_result <- titerlm(
  titers1 = rawmeasuredtiters1,
  titers2 = rawmeasuredtiters2,
  ci_method = "quap",
  dilution_stepsize = 1,
  return_all_pars = TRUE
)
print(titerlm_result[1:5,])

# plot(titers1, titerlm_result[-(1:5), 1], xlim = c(-5, 15), ylim = c(-5, 15))
# titersubset <- measuredtiters2 == -1 & measuredtiters1 == -1
# points(titers1[titersubset], titerlm_result[-(1:5), 1][titersubset], col = "red")
# abline(0, 1, lty = 2)
# abline(h = -1)
# stop()

# Plot the results
plot(
  x = measuredtiters1,
  y = measuredtiters2,
  xlim = c(-1, 13),
  ylim = c(-1, 13),
  col = "#00000033",
  cex = 0.5
)
abline(0, 1, lty = 3)

abline(titerintercept, titerslope, col = "blue", lwd = 2)
abline(lm_result, col = "red", lwd = 2, lty = 2)
abline(log2diff_result["mean", "estimate"], 1, col = "green", lwd = 2, lty = 2)
abline(titerlm_result["intercept", "estimate"], titerlm_result["slope", "estimate"], col = "purple", lwd = 2, lty = 2)
