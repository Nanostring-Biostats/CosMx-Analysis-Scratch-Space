rm(list = ls())
library(InSituType)

# example profiles:
data(ioprofiles)

# load the 6kplex content:
barcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))
genes <- barcodes$Hs_6k$gene
genes <- setdiff(genes, "INS")
ioprofiles <- ioprofiles[is.element(rownames(ioprofiles), genes), ]

# how much signal and background to sim (based on a typical 6k dataset):
totcounts <- 698.039    
#meanneg <- 0.01805              
propactuallybg <- 0.15514892
totrealcounts <- totcounts * (1 - propactuallybg)
meanneg <- totcounts * propactuallybg / length(genes)

# simulate counts with and without background:
set.seed(0)
sim.no.bg <- rmultinom(n = 1, size = totrealcounts, prob = ioprofiles[, "macrophage"])
sim.bg <- rpois(nrow(ioprofiles), lambda = meanneg)
sim.yes.bg <- rmultinom(n = 1, size = totrealcounts, prob = ioprofiles[, "macrophage"]) + sim.bg

# what's the "FDR"?
sum(sim.no.bg) # 589
sum(sim.yes.bg) # 676
sum(sim.yes.bg) - sum(sim.no.bg) # 87
sum(sim.bg) / sum(sim.yes.bg) # 0.129

# barplot of one cell:
png("1cell_barplot.png", width = 8.5, height = 3, units = "in", res = 500)
par(mar = c(4,5,1,0))
bp = barplot(t(cbind(sim.no.bg, sim.bg)), cex.lab = 0.9, ylim = c(0, 15.5),
        col = c('cornflowerblue', "red"), border = c('cornflowerblue', "red"),
        names.arg = rep("", length(sim.bg)), ylab = "Counts in one cell (simulated)", xlab = "6k panel genes")
legend("topright", lty = 1, col = c('cornflowerblue', "red"), legend = c("True counts", "Background counts"))
showgenes <- rownames(sim.no.bg)[order(sim.no.bg, decreasing = T)[1:7]]
text(bp[match(showgenes, rownames(sim.no.bg))], (sim.no.bg + sim.bg)[showgenes, 1], showgenes, cex = 0.6)
dev.off()


# scatterplots vs. truth:
png("1cell_scatterplots.png", width = 8, height = 8, units = "in", res = 500)
par(mfrow = c(2,2))
par(mar = c(5,5,0.1,0.1))

expectedprofile <- ioprofiles[, "macrophage"] / sum(ioprofiles[, "macrophage"]) * totcounts
# no background, linear-scale:
plot(expectedprofile, jitter(sim.no.bg, amount = 0.3),
     col = scales::alpha("darkblue", 0.5), pch = 16, cex.lab = 1.2,
     xlab = "Expected counts", ylab = "Simulated cell, no background")
legend("topleft", legend = paste0("cor = ", round(cor(ioprofiles[, "macrophage"], sim.no.bg), 2)))
# yes background, linear:
plot(expectedprofile, jitter(sim.no.bg+ sim.bg, amount = 0.3),
     col = scales::alpha("darkblue", 0.5), pch = 16, cex.lab = 1.2,
     xlab = "Expected counts", ylab = "Simulated cell with background")
legend("topleft", legend = paste0("cor = ", round(cor(ioprofiles[, "macrophage"], sim.no.bg+ sim.bg), 2)))
# no background, log-scale:
plot(pmax(expectedprofile, 1e-3), jitter(pmax(sim.no.bg, 0.5), amount = 0.1), log = "xy",
     col = scales::alpha("darkblue", 0.5), pch = 16, cex.lab = 1.2,
     xlab = "Expected counts", ylab = "Simulated cell, no background", yaxt = "n")
axis(2, at = c(0.5,1,2,5,10), labels = c(0,1,2,5,10))
legend("topleft", legend = paste0("cor = ", round(cor(
  log2(pmax(ioprofiles[, "macrophage"], 0.5)), 
  log2(pmax(sim.no.bg, 0.5)))
  , 2)))
# yes background, log-scale:
plot(pmax(expectedprofile, 1e-3), jitter(pmax(sim.no.bg+ sim.bg, 0.5), amount = 0.1), log = "xy",
     col = scales::alpha("darkblue", 0.5), pch = 16, cex.lab = 1.2,
     xlab = "Expected counts", ylab = "Simulated cell with background", yaxt = "n")
axis(2, at = c(0.5,1,2,5,10), labels = c(0,1,2,5,10))
legend("topleft", legend = paste0("cor = ", round(cor(
  log2(pmax(ioprofiles[, "macrophage"], 0.5)), 
  log2(pmax(sim.no.bg+ sim.bg, 0.5)))
  , 2)))
dev.off()


# simulate a cell with 20k counts:
set.seed(0)
sim.highcount.real <- rmultinom(n = 1, size = 20000, prob = ioprofiles[, "macrophage"])
sim.highcount.bg <- rpois(nrow(ioprofiles), lambda = meanneg * 20000 / totcounts)
sim.highcount <- sim.highcount.real + sim.highcount.bg

sum(sim.highcount.real)
sum(sim.highcount.bg)
sum(sim.highcount.bg) / sum(sim.highcount)

expectedprofilehigh <- ioprofiles[, "macrophage"] / sum(ioprofiles[, "macrophage"]) * 20000

png("1cell_scatterplots_highcount.png", width = 8, height = 4, units = "in", res = 500)
par(mfrow = c(1,2))
par(mar = c(5,5,0.1,0.1))
# no background, linear-scale:
plot(expectedprofilehigh, jitter(sim.highcount, amount = 0.3),
     col = scales::alpha("darkblue", 0.5), pch = 16, cex.lab = 1.2,
     xlab = "Expected counts", ylab = "Simulated high-count cell")
legend("topleft", legend = paste0("cor = ", round(cor(ioprofiles[, "macrophage"], sim.highcount), 2)))
# no background, log-scale:
plot(pmax(expectedprofilehigh, 1e-3), jitter(pmax(sim.highcount, 0.5), amount = 0.1), log = "xy",
     col = scales::alpha("darkblue", 0.5), pch = 16, cex.lab = 1.2,
     xlab = "Expected counts", ylab = "Simulated high-count cell", yaxt = "n", xaxt = "n")
legend("topleft", legend = paste0("cor = ", round(cor(
  log2(pmax(ioprofiles[, "macrophage"], 0.5)), 
  log2(pmax(sim.highcount, 0.5)))
  , 2)))
axis(2, at = c(0.5,1,2,5,10,20,50,100,200,500), labels = c(0,1,2,5,10,20,50,100,200,500))
axis(1, at = c(1e-3,1e-2,0.1,1,10,100,500), labels = c(1e-3,1e-2,0.1,1,10,100,500))
dev.off()
