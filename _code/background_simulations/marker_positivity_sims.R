# setup: compare sens, spec, and TPR for calling marker gene positivity given 
# expression level (expected counts in the spatial platform), background rate, and 
# expression threshold for calling positivity

library(ggplot2)
rm(list = ls())

# range of expected expression levels:
mu = seq(0.02, 3, by = 0.01)
# background levels:
b = c(0.02, 0.0005)  # 2nd value is 40x lower
# thresholds:
thresh = c(1,2)
# rate of true positivity:
rate = c(.5, .1, 0.01)

# set up data frame of all possible parameters:
df <- expand.grid(mu, b, thresh, rate)
colnames(df) = c("mu", "b", "thresh", "rate")
df$setting = paste0("Background = ", df$b, "; Threshold = ", df$thresh, " counts")
df$setting = gsub("1 counts", "1 count", df$setting)
  
# calc performance metrics:
df$sens <- 1 - ppois(df$thresh - 1, df$mu)
df$spec <- ppois(df$thresh - 1, df$b)
df$tpr <- (df$sens * df$rate) / ((df$sens * df$rate) + ((1 - df$spec) * (1 - df$rate)))

# unstratified:
png("marker detection performance.png", units = "in", width = 8, height = 8, res = 400)
par(mfrow = c(2,2))
par(mar = c(5,6,2,1))
use = TRUE
df$col =  c("blue","red", "dodgerblue2", "orange")[1 + (df$b == b[2]) + 2 * (df$thresh== thresh[2])]

settingcols = c("blue","red", "dodgerblue2", "orange")
names(settingcols) = unique(df$setting)
# sens:
plot(df$mu[use], df$sens[use], col = 0, cex.lab = 1.3,
     xlab = "Expected counts / cell", ylab = "Sensitivity to detect")
use1 = (use & (df$setting == unique(df$setting)[1])) & df$rate == rate[1]
use2 = (use & (df$setting == unique(df$setting)[3])) & df$rate == rate[1]
lines(df$mu[use1], df$sens[use1], col = 1)
lines(df$mu[use2], df$sens[use2], col = 1, lty = 2)
legend("bottomright", lty = c(NA,1:2), legend = c("Count threshold:", thresh))

# tpr:
for (r in rate) {
  use.r = use & (df$rate == r)
  plot(df$mu[use.r], df$tpr[use.r], col = 0, cex.lab = 1.3,
       xlab = "Expected counts / cell", ylab = paste0("True Positive Rate for a cell\ntype with ", r, " prevalence"),
       ylim = c(0, 1)) # pch = c(1 + (df$rate[use.r] == rate[2]) + 2*(df$rate[use] == rate[3])),
  for (setting in unique(df$setting)) {
    use.r.s <- use.r & (df$setting == setting)
    lines(df$mu[use.r.s], df$tpr[use.r.s], col = settingcols[setting], lwd = 1.5)
  }
  if (r == rate[1]) {
    legend("bottomright", lty = 1, col = settingcols, legend = names(settingcols), cex = 0.7)
  }
}
dev.off()
 
