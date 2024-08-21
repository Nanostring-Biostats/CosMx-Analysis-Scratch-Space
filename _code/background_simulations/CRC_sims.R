library(Matrix)

rm(list = ls())
load("CRC excerpt.RData")
ls()

# subset to lower region for concise plotting:
sub = xy[,2] < 4.75
counts = counts[sub, ]
xy = xy[sub, ]
neg = neg[sub]
totcounts = totcounts[sub]
clust = clust[sub]


# simulate a gene with almost no background by down-sampling a high expresser: 
# the mean counts of COL1A1 were 6.34; mean gene had 0.15 mean counts.
# so downsample its counts to the level of an average gene:
set.seed(0)
nobg <- sapply(counts[, "COL1A1"], function(x){rpois(1, x * 0.075 / 6.34)})

# now simulate background atop it, using a negative control gene
bg <- counts[, "negprobe"]
c(mean(bg), mean(neg)) # the negative control gene we use for "background" is representative of the negprobes at large, in fact a bit higher background
# proportion of background counts:
sum(bg) / sum(bg + nobg)

hicols <- colorRampPalette(c("grey80", "darkblue"))(87)
locols <- colorRampPalette(c("grey80", "darkblue"))(4)

png("COL1A1 with downsampling and background.png", width = 8.5, height = 3.75, units = "in", res = 400)
par(mfrow = c(1,4))
par(mar = c(0.1,0.1,3,0.1))
plot(xy, pch = 16, cex = 0.1 + 0.1 * (counts[, "COL1A1"] > 0),
     main = "COL1A1 counts\n(a high expresser)",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = hicols[1 + pmin(length(hicols) - 1, counts[, "COL1A1"])],
     ylim =c(-0.4,4.5), cex.main = 0.9)
plot(xy, pch = 16, cex = 0.1 + 0.1 * (nobg > 0), 
     main = "COL1A1 downsampled, simulating\na low expresser with near-zero background",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = locols[1 + pmin(length(locols) - 1, nobg)],
     ylim =c(-0.4,4.5), cex.main = 0.9)
plot(xy, pch = 16, cex = 0.1 + 0.1 * (bg > 0), 
     main = "Background",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = locols[1 + pmin(length(locols) - 1, bg)],
     ylim = c(-0.4,4.5), cex.main = 0.9)
plot(xy, pch = 16, cex = 0.1  + 0.1 * (nobg + bg > 0), 
     main = "COL1A1 downsampled,\nplus background",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = locols[1 + pmin(length(locols) - 1, nobg + bg)],
     ylim =c(-0.4,4.5), cex.main = 0.9)
dev.off() 

# traces along a narrow band:
xpositions = seq(min(xy[,1]), max(xy[,1]), length.out = 100)
means <- matrix(NA,4,length(xpositions))
ylim = c(2.3,2.4)
iny <- ((xy[,2] > ylim[1]) & (xy[,2] < ylim[2]))
for (i in 1:length(xpositions) ) {
  inds <- ((xy[,1] > xpositions[i] - 0.05) & (xy[,1] < xpositions[i] + 0.05)) & iny
  means[1, i] <- sum(counts[inds, "COL1A1"])
  means[2, i] <- sum(nobg[inds])
  means[3, i] <- sum(bg[inds])
  means[4, i] <- sum((nobg + bg)[inds])
}


png("COL1A1 with downsampling and background - plus traces.png", width = 8.5, height = 4.5, units = "in", res = 400)
layout(t(matrix(1:8, 4)), heights = rep(c(4,.5,4)))
par(mar = c(0.1,0.1,3,0.1))
plot(xy, pch = 16, cex = 0.1 + 0.1 * (counts[, "COL1A1"] > 0),
     main = "COL1A1 counts\n(a high expresser)",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = hicols[1 + pmin(length(hicols) - 1, counts[, "COL1A1"])],
     ylim =c(-0.4,4.5), cex.main = 0.9)
rect(0, ylim[1], 3.05, ylim[2], col = NULL, border = scales::alpha("red", .5))
plot(xy, pch = 16, cex = 0.1 + 0.1 * (nobg > 0), 
     main = "COL1A1 downsampled, simulating\na low expresser with near-zero background",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = locols[1 + pmin(length(locols) - 1, nobg)],
     ylim =c(-0.4,4.5), cex.main = 0.9)
rect(0, ylim[1], 3.05, ylim[2], col = NULL, border = scales::alpha("red", .5))
plot(xy, pch = 16, cex = 0.1 + 0.1 * (bg > 0), 
     main = "Background",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = locols[1 + pmin(length(locols) - 1, bg)],
     ylim = c(-0.4,4.5), cex.main = 0.9)
rect(0, ylim[1], 3.05, ylim[2], col = NULL, border = scales::alpha("red", .5))
plot(xy, pch = 16, cex = 0.1  + 0.1 * (nobg + bg > 0), 
     main = "COL1A1 downsampled,\nplus background (22.4% of all counts)",
     asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = locols[1 + pmin(length(locols) - 1, nobg + bg)],
     ylim =c(-0.4,4.5), cex.main = 0.9)
rect(0, ylim[1], 3.05, ylim[2], col = NULL, border = scales::alpha("red", .5))
par(mar = c(0.1,0.1,0.1,0.1))
for (i in 1:4){
  if (i == 1) {
    plotylim = range(means[1,])
  } else {
    plotylim = range(means[4,])
  }
  plot(means[i,], col = 0,xaxt = "n", yaxt = "n", ylim = plotylim)
  lines(means[i,], col = "red")
}
dev.off() 

# plot low expressers
png("low expressers.png", width = 4.25, height = 3.75, units = "in", res = 400)
par(mfrow = c(1,2))
par(mar = c(0.1,0.1,3,0.1))
for (gene in c("CCL23", "IGF2")) {
  plot(xy, pch = 16, cex = 0.1 + 0.1 * (counts[, gene] > 0), 
       main = paste0(gene, " counts\n", round(mean(counts[,gene])/mean(neg), 2), "-fold above background\n", 
                     round(100*mean(neg) / mean(counts[,gene]), 0), "% counts from background"),
       asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       col = locols[1 + pmin(length(locols) - 1, counts[, gene])],
       ylim =c(-0.4,4.5), cex.main = 0.8)
}
dev.off() 

#### look at DE:

# arbitrary design variable: top vs. lower half:
context <- xy[,2] > median(xy[,2])
use = clust == "CAF"
# ttests:
mod.nobg <- lm((nobg / totcounts)[use] ~ context[use])
mod.bg <- lm(((nobg + bg) / totcounts)[use] ~ context[use])
tab <- signif(rbind(
  c(mod.nobg$coefficients[2], confint(mod.nobg)[2,], summary(mod.nobg)$coef[2,4]),
  c(mod.bg$coefficients[2], confint(mod.bg)[2,], summary(mod.bg)$coef[2,4])
), 3)
colnames(tab) <- c("Estimate", "Lower 95% conf int", "Upper 95% conf int", "p-value")
rownames(tab) <- c("Downsampled so simulate no background", "With simulated background")
tab
write.csv(tab, file = "DE from 18k fibroblasts.csv")

# now downsample to 1k cells:
set.seed(0)
sub = sample(which(use), 1000)
mod.nobg <- lm((nobg / totcounts)[sub] ~ context[sub])
mod.bg <- lm(((nobg + bg) / totcounts)[sub] ~ context[sub])
tab <- signif(rbind(
  c(mod.nobg$coefficients[2], confint(mod.nobg)[2,], summary(mod.nobg)$coef[2,4]),
  c(mod.bg$coefficients[2], confint(mod.bg)[2,], summary(mod.bg)$coef[2,4])
), 3)
colnames(tab) <- c("Estimate", "Lower 95% conf int", "Upper 95% conf int", "p-value")
rownames(tab) <- c("Downsampled so simulate no background", "With simulated background")
tab
write.csv(tab, file = "DE from 1k fibroblasts.csv")
