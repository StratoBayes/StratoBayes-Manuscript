# Visualise results with artificial data
#
# load package
library(StratoBayes)
#
# Example 1
# 
result2a <- readRDS("results/model-illustration/alignment-2-peaks.rds")
result2b <- readRDS("results/model-illustration/alignment-2-peaks-noPenalty.rds")
result2c <- readRDS("results/model-illustration/alignment-2-second-signal.rds")
#
data2a <- result2a$data
#
# directory for saving images
dir1 <- "figures/validation-figure/"
# figures
width <- 2.1
height <- 4
parOp <- par()
parCustom <- parOp
parCustom$mar <- c(4.1, 4.1, 1.6, 1.1)
parCustom$mgp <- c(2.25, 0.8, 0)
xlim <- range(data2a$signal$value)
ylim <- range(data2a$signal$height)
mainCex <- 1
ptCex <- 2/3
alpha <- 0.5

burnIn = 10000
cl2a <- Cluster(result2a, burnin = burnIn, clusterMethod = "hdbscan", minPts = 100)

figName <- paste0(dir1, "fig2-section1.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 1, mgp = c(2.5, 0.8, 0))
plot(data2a, sites = 1, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, type = "n", yaxt = "n")
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(data2a, sites = 1, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, yaxt = "n")
axis(2, at = -1:7 * pi, labels = c(NA, 0, NA, expression(2*pi), NA, expression(4*pi), NA, expression(6*pi), NA))
dev.off()

figName <- paste0(dir1, "fig2-section2.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 1, mgp = c(2.5, 0.8, 0))
plot(data2a, sites = 2, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, type = "n", yaxt = "n")
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(data2a, sites = 2, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, yaxt = "n")
axis(2, at = -1:7 * pi, labels = c(NA, 0, NA, expression(2*pi), NA, expression(4*pi), NA, expression(6*pi), NA))

dev.off()

figName <- paste0(dir1, "fig2-aligned-section-al1.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 1, mgp = c(2.5, 0.8, 0))
plot(result2a, xlim = xlim, ylim = ylim, main = "alignment 1", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites = F, type = "n", yaxt = "n",
     stratCluster = cl2a)
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(result2a, xlim = xlim, ylim = ylim, main = "alignment 1", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites = F, yaxt = "n",
     stratCluster = cl2a)
axis(2, at = -1:7 * pi, labels = c(NA, 0, NA, expression(2*pi), NA, expression(4*pi), NA, expression(6*pi), NA))

dev.off()

figName <- paste0(dir1, "fig2-aligned-section-al2.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 1, mgp = c(2.5, 0.8, 0))
plot(result2a, xlim = xlim, ylim = ylim, main = "alignment 2", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites= F, alignment = 2, type = "n", yaxt = "n", stratCluster = cl2a)
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(result2a, xlim = xlim, ylim = ylim, main = "alignment 2", overridePar = F,
     cex = ptCex, alpha = alpha, alignment = 2, separateSites = F, yaxt = "n",
      stratCluster = cl2a)
axis(2, at = -1:7 * pi, labels = c(NA, 0, NA, expression(2*pi), NA, expression(4*pi), NA, expression(6*pi), NA))

dev.off()

width <- 2.1
height <- 3.4
figName <- paste0(dir1, "fig2-beanplot-alpha.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 1, mgp = c(2.5, 0.8, 0))

BeanPlot(result2a, parameters = 1, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
         xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 10000, yaxt = "n", burnIn = burnIn, stratCluster = cl2a)

abline(h = c(2 * pi, 4 * pi), lty = 3, lwd = 1.5)
axis(2, at = -1:7 * pi, labels = c(NA, 0, NA, expression(2*pi), expression(3*pi), expression(4*pi), NA, expression(6*pi), NA))

dev.off()

figName <- paste0(dir1, "fig2-beanplot-gamma.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 0, mgp = c(3, 0.8, 0))

BeanPlot(result2a, parameters = 2, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(" "~gamma),
         xaxt = "n", srt = 0, adj = 0.5, ylab = NA, maxIterations = 10000, yaxt = "n", burnIn = burnIn, stratCluster = cl2a)
abline(h = log(2), lty = 3, lwd = 1.5)
axis(2, at = log(c(1.6, 1.7, 1.8, 1.9, 2., 2.1)), labels = c(NA, "ln(1.7)", "ln(1.8)", "ln(1.9)", "ln(2)", NA),
     cex.axis = 0.9)

dev.off()

iterations <- seq(800, 60000, 800)

width <- 4.2
height <- 2.2
figName <- paste0(dir1, "fig2-traceplot-alpha.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 1, mgp = c(2.5, 0.8, 0))

TracePlot(result2a, parameters = 1, iterations = iterations, overridePar = FALSE, ylab = expression(alpha), col = hcl.colors(6, "Plasma", alpha = seq(0.5, 0.65, length.out = 6))[c(1,3,5)],
          srt = 0, adj = 0.5, type = "o", cex = .75, lwd = 1, yaxt = "n")
axis(2, at = -1:7 * pi, labels = c(NA, 0, NA, expression(2*pi), expression(3*pi), expression(4*pi), NA, expression(6*pi), NA))

dev.off()

figName <- paste0(dir1, "fig2-traceplot-gamma.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
par(las = 0, mgp = c(2.5, 0.8, 0))

TracePlot(result2a, parameters = 2, iterations = iterations, overridePar = FALSE, ylab = NA,
          srt = 0, adj = 0.5, type = "o", cex = .75, lwd = 1, yaxt = "n",
          col = hcl.colors(6, "Plasma", alpha = seq(0.5, 0.65, length.out = 6))[c(1,3,5)])
axis(2, at = log(c(1.6, 1.7, 1.8, 1.9, 2., 2.1)), labels = c(NA, "ln(1.7)", "ln(1.8)", "ln(1.9)", "ln(2)", NA),
     cex.axis = 0.9)

dev.off()


### Extra bean plot no overlap penalty

width <- 2.1
height <- 3.4
ylim <- c(5.75, 13.5)

col1 <- rgb(0, 0, 0, 0.25)
col2 <- rgb(0.75, 0.25, 0, 0.33)


b1 <- BeanPlot(result2a, parameters = 1, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)

b2 <- BeanPlot(result2b, parameters = 1, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)


y1 <- b1[[1]][[1]]$y
y2 <- b2[[1]][[1]]$y
x1 <- b1[[1]][[1]]$x
x2 <- b2[[1]][[1]]$x
max12 <- max(y1, y2)
indb1 <- which(y1 > max12/1000)
indb1a <- indb1[1:which.max(diff(indb1))]
indb1b <- indb1[(which.max(diff(indb1))+1):length(indb1)]

indb2 <- which(y2 > max12/1000)
indb2a <- indb2[1:which.max(diff(indb2))]
indb2b <- indb2[(which.max(diff(indb2))+1):length(indb2)]

figName <- paste0(dir1, "fig2-beanplot-alpha-no-overlapPen.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(0, type = "n", xlim = c(-1,1) * max12, ylim = ylim)

abline(h = c(2, 4) * pi, lty = 3, lwd = 1.5)

polygon(c(-y1[indb1a], y1[rev(indb1a)]), x1[c(indb1a, rev(indb1a))], col = col1)
polygon(c(-y1[indb1b], y1[rev(indb1b)]), x1[c(indb1b, rev(indb1b))], col = col1)

polygon(c(-y2[indb2a], y2[rev(indb2a)]), x2[c(indb2a, rev(indb2a))], col = col2)
polygon(c(-y2[indb2b], y2[rev(indb2b)]), x2[c(indb2b, rev(indb2b))], col = col2)
dev.off()

width <- 2.1
height <- 3.4
ylim <- c(0.54, 0.76)

b1 <- BeanPlot(result2a, parameters = 2, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)

b2 <- BeanPlot(result2b, parameters = 2, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)


y1 <- b1[[1]][[1]]$y
y2 <- b2[[1]][[1]]$y
x1 <- b1[[1]][[1]]$x
x2 <- b2[[1]][[1]]$x
max12 <- max(y1, y2)
indb1 <- which(y1 > max12/1000)
indb1a <- indb1[1:which.max(diff(indb1))]
indb1b <- indb1[(which.max(diff(indb1))+1):length(indb1)]

indb2 <- which(y2 > max12/1000)
indb2a <- indb2[1:which.max(diff(indb2))]
indb2b <- indb2[(which.max(diff(indb2))+1):length(indb2)]

figName <- paste0(dir1, "fig2-beanplot-gamma-no-overlapPen.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(0, type = "n", xlim = c(-1,1) * max12, ylim = ylim)
abline(h = log(2), lty = 3, lwd = 1.5)

polygon(c(-y1[indb1a], y1[rev(indb1a)]), x1[c(indb1a, rev(indb1a))], col = col1)
polygon(c(-y1[indb1b], y1[rev(indb1b)]), x1[c(indb1b, rev(indb1b))], col = col1)

polygon(c(-y2[indb2a], y2[rev(indb2a)]), x2[c(indb2a, rev(indb2a))], col = col2)
polygon(c(-y2[indb2b], y2[rev(indb2b)]), x2[c(indb2b, rev(indb2b))], col = col2)
dev.off()



# 
# 
# figName <- paste0(dir1, "fig2-beanplot-alpha.png")
# png(figName, width = width, height = height, res = 300, units = "in")
# suppressWarnings(par(parCustom))
# BeanPlot(result2a, parameters = 1, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
#          xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
#          ylim = ylim)
# 
# abline(h = c(2 * pi, 4 * pi), lty = 3, lwd = 1.5)
# par(new = TRUE)
# BeanPlot(result2b, parameters = 1, col = rgb(0.75,0,0,0.33), overridePar = FALSE, xlab = expression(" "~gamma),
#          xaxt = "n", srt = 0, adj = 0.5, ylab = "relative sedimentation rate", maxIterations = 10000,
#          ylim = ylim)
# abline(h = log(2), lty = 3, lwd = 1.5)
# dev.off()
# 
# 
# ylim <- c(0.54, 0.76)
# 
# figName <- paste0(dir1, "fig2-beanplot-gamma-no-overlapPen.png")
# png(figName, width = width, height = height, res = 300, units = "in")
# suppressWarnings(par(parCustom))
# BeanPlot(result2a, parameters = 2, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(" "~gamma),
#          xaxt = "n", srt = 0, adj = 0.5, ylab = "relative sedimentation rate", maxIterations = 10000,
#          ylim = ylim)
# abline(h = log(2), lty = 3, lwd = 1.5)
# par(new = TRUE)
# BeanPlot(result2b, parameters = 2, col = rgb(0.75,0,0,0.33), overridePar = FALSE, xlab = expression(" "~gamma),
#          xaxt = "n", srt = 0, adj = 0.5, ylab = "relative sedimentation rate", maxIterations = 10000,
#          ylim = ylim)
# abline(h = log(2), lty = 3, lwd = 1.5)
# dev.off()


# figures
width <- 2.1
height <- 4
parOp <- par()
parCustom <- parOp
parCustom$mar <- c(4.1, 4.1, 1.6, 1.1)
parCustom$mgp <- c(2.25, 0.8, 0)
xlim <- range(data2a$signal$value)
ylim <- range(data2a$signal$height)
mainCex <- 1
ptCex <- 2/3
alpha <- 0.5

figName <- paste0(dir1, "fig2-aligned-section-al1.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(result2a, xlim = xlim, ylim = ylim, main = "alignment 1", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites = F, type = "n", alignment = 2)
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(result2a, xlim = xlim, ylim = ylim, main = "alignment 1", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites = F, alignment = 2)
dev.off()


####
### second signal
####

data2c <- result2c$data
#
# directory for saving images
dir1 <- "figures/validation-figure/"
# figures
width <- 2.1
height <- 4
parOp <- par()
parCustom <- parOp
parCustom$mar <- c(4.1, 4.1, 1.6, 1.1)
parCustom$mgp <- c(2.25, 0.8, 0)
xlim <- range(c(data2c$signal$signal1, data2c$signal$signal2))
ylim <- range(data2c$signal$height)
mainCex <- 1
ptCex <- 2/3
alpha <- 0.5


cols <- c(hcl.colors(2, "Dark 3", alpha = 0.5)[1],
  rgb(0.75, 0.5, 0.25, 0.5),
  hcl.colors(2, "Dark 3", alpha = 0.5)[2],
  rgb(0.05, 0.33, 0.75, 0.5))
pchs <- c(21, 23, 22, 24)

plot(1:4, col = cols, cex = 4, pch = 19)

figName <- paste0(dir1, "fig2signal-section1a.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(data2c, sites = 1, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, type = "n")
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(data2c, sites = 1, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, col = cols[1], pch = pchs[1])
dev.off()

figName <- paste0(dir1, "fig2signal-section1b.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(data2c, sites = 1, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, type = "n", signal = 2)
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(data2c, sites = 1, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, signal = 2, 
     col = cols[2], pch = pchs[2])
dev.off()

figName <- paste0(dir1, "fig2signal-section2a.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(data2c, sites = 2, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, type = "n", signal = 1)
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(data2c, sites = 2, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, signal = 1, 
     col = cols[3], pch = pchs[3])
dev.off()

figName <- paste0(dir1, "fig2signal-section2b.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(data2c, sites = 2, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, type = "n", signal = 2)
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(data2c, sites = 2, xlim = xlim, ylim = ylim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F, signal = 2, 
     col = cols[4], pch = pchs[4])
dev.off()



figName <- paste0(dir1, "fig2signal-aligned-section-al1.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(result2c, xlim = xlim, ylim = ylim, main = "signal 1", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites = F, type = "n")
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(result2c, xlim = xlim, ylim = ylim, main = "signal 1", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites = F)
dev.off()

figName <- paste0(dir1, "fig2signal-aligned-section-al2.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(result2c, xlim = xlim, ylim = ylim, main = "alignment 2", overridePar = F,
     cex = ptCex, alpha = alpha, separateSites= F, signal = 2, type = "n")
abline(h = c(0, 1, 2, 3) * 2 * pi, lty = 3, lwd = 1.5)
par(new = TRUE)
plot(result2c, xlim = xlim, ylim = ylim, main = "alignment 2", overridePar = F,
     cex = ptCex, alpha = alpha, signal = 2, separateSites = F,
     col = cols[c(2, 4)], pch = pchs[c(2,4)])
dev.off()



####
### Extra bean plot 2 signals
####

width <- 2.1
height <- 4
ylim <- c(5.75, 13.5)

col1 <- rgb(0, 0, 0, 0.25)
col2 <- rgb(0, 0.67, 0.2, 0.33)


b1 <- BeanPlot(result2a, parameters = 1, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)

b2 <- BeanPlot(result2c, parameters = 1, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)


y1 <- b1[[1]][[1]]$y
y2 <- b2[[1]][[1]]$y
x1 <- b1[[1]][[1]]$x
x2 <- b2[[1]][[1]]$x
max12 <- max(y1, y2)
indb1 <- which(y1 > max12/1000)
indb1a <- indb1[1:which.max(diff(indb1))]
indb1b <- indb1[(which.max(diff(indb1))+1):length(indb1)]

indb2 <- which(y2 > max12/1000)
indb2a <- indb2[1:which.max(diff(indb2))]
indb2b <- indb2[(which.max(diff(indb2))+1):length(indb2)]

figName <- paste0(dir1, "fig2signal-beanplot-alpha.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(0, type = "n", xlim = c(-1,1) * max12, ylim = ylim)

abline(h = c(2, 4) * pi, lty = 3, lwd = 1.5)

polygon(c(-y1[indb1a], y1[rev(indb1a)]), x1[c(indb1a, rev(indb1a))], col = col1)
polygon(c(-y1[indb1b], y1[rev(indb1b)]), x1[c(indb1b, rev(indb1b))], col = col1)

polygon(c(-y2[indb2a], y2[rev(indb2a)]), x2[c(indb2a, rev(indb2a))], col = col2)
polygon(c(-y2[indb2b], y2[rev(indb2b)]), x2[c(indb2b, rev(indb2b))], col = col2)
dev.off()

width <- 2.1
height <- 4
ylim <- c(0.5025, 0.715)

b1 <- BeanPlot(result2a, parameters = 2, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)

b2 <- BeanPlot(result2c, parameters = 2, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = expression(alpha),
               xaxt = "n", srt = 0, adj = 0.5, ylab = "reference height", maxIterations = 100,
               ylim = ylim)


y1 <- b1[[1]][[1]]$y
y2 <- b2[[1]][[1]]$y
x1 <- b1[[1]][[1]]$x
x2 <- b2[[1]][[1]]$x
max12 <- max(y1, y2)
indb1 <- which(y1 > max12/1000)
indb1a <- indb1[1:which.max(diff(indb1))]
indb1b <- indb1[(which.max(diff(indb1))+1):length(indb1)]

indb2 <- which(y2 > max12/1000)
indb2a <- indb2[1:which.max(diff(indb2))]
indb2b <- indb2[(which.max(diff(indb2))+1):length(indb2)]

figName <- paste0(dir1, "fig2signal-beanplot-gamma.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(0, type = "n", xlim = c(-1,1) * max12, ylim = ylim)
abline(h = log(2), lty = 3, lwd = 1.5)

polygon(c(-y1[indb1a], y1[rev(indb1a)]), x1[c(indb1a, rev(indb1a))], col = col1)
polygon(c(-y1[indb1b], y1[rev(indb1b)]), x1[c(indb1b, rev(indb1b))], col = col1)

polygon(c(-y2[indb2a], y2[rev(indb2a)]), x2[c(indb2a, rev(indb2a))], col = col2)
polygon(c(-y2[indb2b], y2[rev(indb2b)]), x2[c(indb2b, rev(indb2b))], col = col2)
dev.off()
