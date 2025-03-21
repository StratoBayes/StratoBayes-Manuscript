# Visualise results with Cambrian data from Morocco and Siberia
#
# load package
library(StratoBayes)
# 
# resultMorocco <- readRDS("D://OneDrive - Durham University/Strat/alignment-tests/resultMorocco3sites100k3.rds")

resC <- readRDS("D://Github/StratobayesManuscript/results/cambrian/cambrian-results.rds")
# resC <- StratPosterior("D://OneDrive - Durham University/Strat/alignment-tests/run-2024-12-11-2 - Copy")
modC <- resC$model # readRDS("results/cambrian/cambrian-model.rds")
datC <- resC$data # readRDS("results/cambrian/cambrian-data.rds")
  
clHdbScan <- readRDS("results/cambrian/cambrian-results-clusterHdbscan.rds")

burnIn <- resC$mcmcSummary$endAdapting[[1]]
fullSummaryResC2 <- summary(resC, burnIn = burnIn, alignment = "none", maxSamples = Inf)

siteCol <- 
c("#935f24",
  "#de7b9f",
  "#a6bc8a",
  "#0080b3")
DecreaseAlpha <- function(colours, alpha = 0.5) {
  rgbValues <- col2rgb(colours) / 255  # Convert to RGB values between 0 and 1
  apply(rgbValues, 2, function(col) rgb(col[1], col[2], col[3], alpha = alpha))
}
siteColTrans <- DecreaseAlpha(siteCol, alpha = 0.5)

siteColDarker <- c("#8C5B22",
                   "#B8537B",
                   "#5B9451",
                   "#0080b3"
                   
)
siteColDarkerTrans <- c(DecreaseAlpha(siteColDarker[1], alpha = 0.67),
                   DecreaseAlpha(siteColDarker[2], alpha = 0.67),
                   DecreaseAlpha(siteColDarker[3], alpha = 0.67),
                   DecreaseAlpha(siteColDarker[4], alpha = 0.67))

siteColDarkerTrans2 <- c(DecreaseAlpha(siteColDarker[1], alpha = 0.85),
                        DecreaseAlpha(siteColDarker[2], alpha = 0.85),
                        DecreaseAlpha(siteColDarker[3], alpha = 0.85),
                        DecreaseAlpha(siteColDarker[4], alpha = 0.85))

moroccoCol <- c(
  "#bf51ad", # Tifnout strom
  "#864063", # 
  "#da95a7", #
  "#c9415c",
  "#5e7e45",
  "#8a5d39", 
  "#ccb357",
  "#cf6436", 
  "#7dc967"
)
  
  
siberiaCol <- c("#2d6389",
                "#49afcd")
partsCol <- c(moroccoCol, siberiaCol)

morCol <- "#a27066"
sibCol <- "#5093ba"

continentCol <- c(morCol, sibCol)
  
continentColTrans <- DecreaseAlpha(continentCol, alpha = 0.5)

partsColtrans <- DecreaseAlpha(partsCol, alpha = 0.5)
partsColtrans30 <- DecreaseAlpha(partsCol, alpha = 0.30)
partsColtrans40 <- DecreaseAlpha(partsCol, alpha = 0.40)

siteOrder <- c(1,3,2,4)
partitionOrder <- c(9, 7, 6, 5, 4, 3, 2, 1, 8, 11, 10)
#
# directory for saving images
dir1 <- "figures/cambrian-figure/"
dir2 <- "figures/cambrian-figure/appendix/"

# figures
width <- 2.2
height <- 4.3
parOp <- par()
parCustom <- parOp
parCustom$mar <- c(4.1, 4.1, 1.6, 1.1)
parCustom$mgp <- c(2.25, 0.8, 0)
xlim <- range(datC$signal$d13C) + diff(range(datC$signal$d13C)) * c(-0.015, 0.03)
ylim <- range(datC$signal$height)
mainCex <- 1
ptCex <- 7/10
alpha <- 0.5
type = "p"
xlab = expression(delta^13*C)
ylab = "height (m)"
mainLab = NA
sitenames <- c("Tiout", "Oued Sdas",  "Talat n'Yssi", "Sukharikha")
sitesAbbr <- c("MS6", "MS7",  "MS3", "Sib")

lineCol <- "grey80"

# liths <- c("laminite / siltstone",  "grainstone / sandstone", "microbialaminite", "stromatolite",
#            "lower Sukharikha fm.", "upper Sukharika fm.", "Krasnoporog fm.")
liths <- c("Tifnout stromatolite", "Lie de Vin (l.)",
           "Lie de Vin (m.)", "Lie de Vin (u.)",
           "Igoudine (l.)", "Igoudine (u.)",
           "Amouslek", "Tifnout (l.)",
           "Isaafen",
           "Sukharikha", "Krasnoporog")


for (site in c(4,3,2,1)) {
figName <- paste0(dir1, paste0("figC-", sitesAbbr[site], ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))

plot(datC, sites = site, xlim = xlim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F,
     xlab = xlab, ylab = ylab, main = mainLab,
     col = partsCol, type = "n")
abline(v = c(-4,0,4), col = lineCol)

par(new = TRUE)
plot(datC, sites = site, xlim = xlim, cex.main = mainCex,
     cex = ptCex, alpha = alpha, overridePar = F,
     xlab = xlab, ylab = ylab, type = type, main = mainLab,
     col = partsCol)

if (site %in% c(1, 2, 3)) points(rep(xlim[2] + 0.04*diff(xlim), datC$tiesAttr$tiesN[site]), 
                            datC$tiesAttr$tiesH[1:datC$tiesAttr$tiesN[site],site],
                            pch = 1, xpd = TRUE, cex = 1.1)

# if (site == 4) points(rep(xlim[1] - 0.04*diff(xlim), datC$tiesAttr$tiesN[site]), 
#                              datC$tiesAttr$tiesH[1:datC$tiesAttr$tiesN[site],site],
#                              pch = 4, xpd = TRUE)
if (site == 4) points(c(xlim[1] - 0.04*diff(xlim), xlim[2] + 0.04*diff(xlim)),
                        rep(datC$partsAttr$Bound[2,site], 2),
                      type = "l", lty = 2, col = "grey33", lwd = 1.5)

points(xlim[2] + 0.075 * diff(range(xlim)), modC$alphaPosition[site], xpd = T, pch = "–",
       cex = 1, col = "grey50")
dev.off()
}

# legend
figName <- paste0(dir1, paste0("figC-", "legend", ".png"))
png(figName, width = width*1.5, height = 1.075 * height, res = 300, units = "in")
suppressWarnings(par(parCustom))
nCol <- 11
plot(rep(1,nCol), 14:4, pch = c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)[partitionOrder],
     bg = partsCol[partitionOrder],
     col = partsCol[partitionOrder],
     ylim = c(-1, 15), xlim = c(0.8, 1.6))
points(c(0.94, 1.06), c(3,3), type = "l", lty = 2, lwd = 1.25)
points(1, 2, type = "p", pch = 1, cex = 1.2)
points(1, 1, type = "p", pch = NA)

text(rep(1.08, 10), 14:0, c(liths[partitionOrder],
                         "sequence boundary",
                         "radiometric date",
                         "reference horizon",
                         "first trilobite fossils"),
     cex = 0.66, adj = 0)

dev.off()



# Priors and posteriors
width <- 3.5
height <- 3.2
xlab <- expression(alpha~"(Ma)")
figName <- paste0(dir1, paste0("figP-alpha", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(modC$priors, overridePar = F, parameters = c(1:3),
     xlim = c(534.2, 508), xaxt = "n", xlab = xlab,
     yaxt = "n", col = c(rgb(0,0,0,0.25), siteColTrans[2], siteColTrans[3]), alpha = 0.5)
axis(2, seq(0,0.3,0.1))
axis(1, seq(540,505,-5), labels = rep(c(1,NA),4) * seq(540,505,-5))
dev.off()

figName <- paste0(dir1, paste0("figP-gamma1", ".png"))
xlab <- expression(gamma~"(m/Myr)")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(modC$priors, overridePar = F, parameters = 6:10, log = T,
xaxt = "n", xlab = xlab, yaxt = "n", col = partsColtrans[c(5,1,2,3,4) + 1])#,
polygon(log(c(175, 630, 630, 175)),
        c(-0.23, -0.23, -0.06, -0.06),
        col = rgb(.8, .3, .3, .75), border = NA)
axis(2, seq(0,10,2))
axisTx <- c(200, 250, 350, 500)
axis(1, log(axisTx), labels = axisTx)
dev.off()


figName <- paste0(dir1, paste0("figP-gamma2", ".png"))
xlab <- expression(gamma~"(m/Myr)")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(modC$priors, overridePar = F, parameters = c(5, 14), log = T,
     xaxt = "n", xlab = xlab, yaxt = "n", col = continentColTrans)
polygon(log(c(175, 630, 630, 175)),
        c(-0.0165, -0.0165, -0.004, -0.004),
        col = rgb(.8, .3, .3, .75), border = NA, )
axis(2, seq(0,0.5,0.1))
axisTx <- c(1, 5, 25, 100, 500, 3000)
axis(1, log(axisTx), labels = axisTx)
dev.off()


figName <- paste0(dir1, paste0("figP-zeta", ".png"))
xlab <- expression(zeta)
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(modC$priors, overridePar = F, parameters = 16, log = T,
     xaxt = "n", xlab = xlab, yaxt = "n", col = "grey60")#,
axis(2, seq(0,2,0.5))
axisTx <- c(1/2, 3/4, 1, 1.5, 2)
axis(1, log(axisTx), labels = axisTx)
dev.off()

figName <- paste0(dir1, paste0("figP-gap", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(modC$priors, overridePar = F, parameters = c(19),
     col = "grey75", yaxt = "n",
      xlab = expression(delta~"(Myr)"))
axis(2, seq(0,1, 0.25), labels = c(0, NA, 0.5, NA, 1))
dev.off()


## alignment plot


xlab = expression(delta^13*C)
xlim <- range(datC$signal$d13C) + diff(range(datC$signal$d13C)) * c(-0.015, 0.03)

height <- 6.5
width <- 2.9
al <- 0
ptCex <- 2/3
figName <- paste0(dir1, paste0("figCambrianAlignment",al, ".png"))
png(figName, width = width, height = height, res = 300, units = "in")

ylim <- c(535, 516.75)

par(mfrow = c(1,1))

plot(resC, quantile = c(0.5), alignment = al, ylim = ylim, separateSites = F,
     overridePar = F, ylab = "age (Ma)",
     xlab = xlab, cex = ptCex, col = siteColTrans, type = "n", burnIn = burnIn,
     xlim = xlim, maxSamples = Inf,
     xaxt = "n", yaxs = "i", stratCluster = clHdbScan)

axis(1, seq(-8, 8, 2), c(NA, NA, -4, NA, 0, NA, 4, NA, 8))

abline(v = c(-4, 0, 4), col = "grey80")
abline(h = seq(532.5, 515, -2.5), lty = 1, col =  "grey80")

par(new = TRUE)

plot(resC, quantile = c(0.5), alignment = al, ylim = ylim, separateSites = F,
     overridePar = F, ylab = "age (Ma)",
     xlab = xlab, cex = ptCex, col = siteColTrans, burnIn = burnIn,
     xlim = xlim, maxSamples = Inf,
     xaxt = "n", yaxs = "i", stratCluster = clHdbScan)

dev.off()


## prior alignments
height <- 4.8
width <- 2.3
 set.seed(1)
 prior1 <- sapply(1:19, function(i) do.call(resC$model$priors$metro[[i]]$R, c(list(n = 1),  resC$model$priors$metro[[i]]$args)))
 set.seed(10)
 prior2 <- sapply(1:19, function(i) do.call(resC$model$priors$metro[[i]]$R, c(list(n = 1),  resC$model$priors$metro[[i]]$args)))
 
 
 figName <- paste0(dir1, paste0("figCambrianAlignmentprior_1", ".png"))
 png(figName, width = width, height = height, res = 300, units = "in")
 
 ylim <- NULL
 
 par(mfrow = c(1,1), mgp = c(2.25, 0.75, 0))
 
 plot(resC, parameters = prior1, separateSites = F,
      overridePar = F, ylab = "age (Ma)",
      xlab = xlab, cex = ptCex, col = siteColTrans, type = "n", burnIn = burnIn,
      xlim = xlim, maxSamples = Inf,
      xaxt = "n", stratCluster = clHdbScan)
 
 axis(1, seq(-8, 8, 2), c(NA, NA, -4, NA, 0, NA, 4, NA, 8))
 
 abline(v = c(-4, 0, 4), col = "grey80")
# abline(h = seq(532.5, 515, -2.5), lty = 1, col =  "grey80")
 
 par(new = TRUE)
 
 plot(resC, parameters = prior1, separateSites = F,
      overridePar = F, ylab = "age (Ma)",
      xlab = xlab, cex = ptCex, col = siteColTrans, burnIn = burnIn,
      xlim = xlim, maxSamples = Inf,
      xaxt = "n", stratCluster = clHdbScan)
 
 dev.off()
 
 
 
 
 # Posterior alignment 1
 ### Slow!
 # metroPar <- ExtractParameters(resC, burnIn = burnIn, stratCluster = clHdbScan, alignment = 1)
 # metroParAll <- ExtractParameters(resC, burnIn = burnIn, stratCluster = clHdbScan, alignment = "all")
 # 
 # convertedHeights <- t(
 #   vapply(seq_len(nrow(metroPar)), function(i) do.call(c,
 #                                                       lapply(modC[["siteIndices"]], function(site) {
 #                                                         StratoBayes:::.AgeConversionGeneral(
 #                                                           parameters = as.numeric(metroPar[i,]),
 #                                                           heights = modC[["heightsForClustering"]][[site]],
 #                                                           site = site,
 #                                                           stratData = datC,
 #                                                           stratModel = modC)
 #                                                         
 #                                                       })), numeric(sum(lengths(modC[["heightsForClustering"]][modC[["siteIndices"]]]))))
 # )
 # 
 # quantileIndex <-
 #   StratoBayes:::.SelectRowClosestToQuantile(convertedHeights, 0.5)
 # 
 # rowIndex <- which(apply(metroParAll, 1, function(x) all(x == metroPar[quantileIndex,])))
 # iter <- (rowIndex - 3 * 12000 ) * 50 + 150000 
 
 iter1 <- 522150
 plot(resC, separateSites = T, colourBy = "p", iterations = iter, runs = 4, sites = 1, col = partsCol,
      overridePar = F)
 
 plot(resC, separateSites = T, colourBy = "p", iterations = iter, runs = 4, sites = 1:4, col = partsCol)
 
 ### lithology bars
topH <- resC$data$partsAttr$Bound 
topH[1,4] <- topH[1,4] - 0.0001 
bottomH <- resC$data$partsAttr$BoundLow
bottomH[2,4] <- bottomH[2,4] - 0.0001 

nBound <- resC$data$partsAttr$nBound
topA <- list()
bottomA <- list()
al <- 1
for(al in 1:3) {
  topA[[al]] <- list()
  bottomA[[al]] <- list()
for (s in 1:4) {
  print(c(al, s))
  topA[[al]][[s]] <- StratMap(resC, heights = topH[1:nBound[s],s], site = s, alignment = al, stratCluster = clHdbScan,
  burnIn = burnIn)[,"50%"]
  bottomA[[al]][[s]] <- StratMap(resC, heights = bottomH[1:nBound[s],s], site = s, alignment = al, stratCluster = clHdbScan,
                              burnIn = burnIn)[,"50%"]
}
}

# saveRDS(list(bottomA = bottomA, topA = topA), "results/cambrian/cambrian-results-bound-ages.rds") 
boundAges <- readRDS("results/cambrian/cambrian-results-bound-ages.rds")

sedRate <- list()
for(al in 1:3) {
  sedRate[[al]] <- list()
  for(s in 1:4) {
    if (s == 1) sedRate[[al]][[s]] <- summaryFull[[al]]$summary[4 + 1:7, "50%"]
    if (s == 2) sedRate[[al]][[s]] <- summaryFull[[al]]$summary[4 + c(8,1:5), "50%"] +  summaryFull[[al]]$summary["zetaLog_MS7", "50%"]
    if (s == 3) sedRate[[al]][[s]] <- summaryFull[[al]]$summary[4 + c(3:7,9), "50%"] +  summaryFull[[al]]$summary["zetaLog_MS3", "50%"]
    if (s == 4) sedRate[[al]][[s]] <- summaryFull[[al]]$summary[c(14, NA, 15), "50%"]
    
  }
}

sedRange <- range(unlist(sedRate), na.rm = T)
sedDiff <- diff(range(unlist(sedRate)))
normalise <- function(v, x, y, a, b) {
  a + ((v - x) * (b - a)) / (y - x)
}
normaliseS <- function(v, x, y, a = 0.5, b = 1) {
  x = sedRange[1]
  y = sedRange[2]
  a + ((v - x) * (b - a)) / (y - x)
}
  

normaliseS(c(sedRange, 6))  

((sedRange - sedRange[1]) * 0.1)/(sedDiff *1.1)

ylim <- c(535, 516.75)


inds <- list(seq_len(nBound[1]),
             seq_len(nBound[1]),
             seq_len(nBound[1]),
             c(1,3)
)


height <- 6.5
width <- 1.6
al <- 0
ptCex <- 2/3
al = 2

figName <- paste0(dir1, paste0("figCambrianAlignment-boundaryAges",al, ".png"))
png(figName, width = width, height = height, res = 300, units = "in")

par(mfrow = c(1,1))

plot(0,ylab = "age (Ma)",
     xlab = xlab,type = "n", 
     xlim = c(0,5), 
     xaxt = "n", yaxs = "i",
     ylim = ylim)
abline(h = seq(532.5, 515, -2.5), lty = 1, col =  "grey80")


for(s in 1:4) {
xPos <- c(1,2,3,4)
xAd <- c(-0.5, 0.5)
numParts <- resC$data$partsAttr$sedRatesDF
sapply(inds[[s]], function(x) polygon(xPos[s] + c(xAd, rev(xAd)), #* normaliseS(sedRate[[al]][[s]][x]), 
                                               c(rep(boundAges[[1]][[al]][[s]][x], 2),
                                                 rep(boundAges[[2]][[al]][[s]][x], 2)),
                                               col = partsCol[numParts[x,s]], border = "black", lwd = 0.75))
}

dev.off()

##
## date plots
##
##

# dateColTrans <- c(rgb(0.8, 0, 0, 0.5), rgb(0.6, 0.3, 0, 0.5),
#                   rgb(0.4, 0.6, 0, 0.5), rgb(0.2, 0.9, 0, 0.5))

xPos <- list(
  c(-0.76, 0, 0.76),
  c(-0.76, -0.26, 0.26, 0.76),
  -0.76
)

height <- 6.5
width <- 2.5
figName <- paste0(dir1, paste0("datesAlignment", ".png"))

png(figName, width = width, height = height, res = 600, units = "in")

plot(0, type = "n", xlim = c(-1.5 + 0.1, 3 + 0.12 - 0.1), ylim = ylim, yaxs = "i", xaxt = "n")
abline(h = seq(532.5, 517.5, -2.5), col = "grey80")
indi <- list(1:3, 1:4, 1)

for(s in 1:3) {
for (d in indi[[s]]) {
densities <- NULL

 mu1 <- datC$tiesAttr$arg1[d,s]
 sd1 <- datC$tiesAttr$arg2[d,s]
# y1 <- seq(mu1 - 3 * sd1, mu1 + 3 * sd1, length.out = 1000)
# pd1 <- dnorm(y1, mu1, sd1)
# pd1 <- pd1 / max(pd1)
# polygon(c(pd1, - rev(pd1)), c(y1, rev(y1)), border = NA, col = partsColtrans30)

xWidth <- c(2/3, 1/5, 1/16)
 
col1 <- "grey77"
col2 <- "grey55"

#lines(c(-1, 1) * xWidth[1] + xPos[[s]][d], c(mu1, mu1), col = col2, lwd = 1.5, lty = 3)
y1 <- seq(mu1 - 3 * sd1, mu1 + 3 * sd1, length.out = 1000)
d1 <- dnorm(y1, mu1, sd1)
d1 <- d1 / max(d1) * 2/3

y1 <- y1
x1 <- c(-d1, d1) + xPos[[s]][d]
polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)

# y1 <- seq(mu1 - 1 * sd1, mu1 + 1 * sd1, length.out = 1000)
# x1 <- c(rep(-xWidth[1], 1000), rep(xWidth[1], 1000)) + xPos[[s]][d]
# polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)


# y1 <- seq(mu1 + 2 * sd1, mu1 + 1 * sd1, length.out = 1000)
# x1 <- c(rep(-xWidth[2], 1000), rep(xWidth[2], 1000)) + xPos[[s]][d]
# polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
# 
# y1 <- seq(mu1 - 3 * sd1, mu1 - 2 * sd1, length.out = 1000)
# x1 <- c(rep(-xWidth[3], 1000), rep(xWidth[3], 1000)) + xPos[[s]][d]
# polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
# 
# y1 <- seq(mu1 + 3 * sd1, mu1 + 2 * sd1, length.out = 1000)
# x1 <- c(rep(-xWidth[3], 1000), rep(xWidth[3], 1000)) + xPos[[s]][d]
# polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)

}
}

for(s in 1:3) {
  for (d in indi[[s]]) {
h1 <- datC$tiesAttr$tiesH[d,s]
a1 <- StratMap(resC, heights = h1, site = s, burnIn = burnIn, quantiles = "samples")

a1b <- as.numeric(a1[,2:2001])
med1 <- median(a1b)
d1 <- density(a1b)
d1$y <- d1$y / max(d1$y) * 2/3
densities[[1]] <- d1
i <- 1
# "bean plot" type density
# plot separate polygons for regions separated by nearly 0 density
densThreshold <- 1 / 1000
meetsThreshold <- densities[[i]][["y"]] >= densThreshold
starts <- which(diff(c(FALSE, meetsThreshold)) == 1)
ends <- which(diff(c(meetsThreshold, FALSE)) == -1) + 1

for (j in seq_along(starts)) {
  argsPoly <- NULL
  argsPoly$border <- NA
  argsPoly$col <- siteColDarkerTrans[[s]]
  indj <- starts[[j]]:ends[[j]]
  argsPoly[["y"]] <- c(densities[[i]][["x"]][indj],
                       rev(densities[[i]][["x"]][indj]))
  argsPoly[["x"]] <- c(densities[[i]][["y"]][indj] + xPos[[s]][d],
                       -rev(densities[[i]][["y"]][indj]) + xPos[[s]][d])
  do.call(polygon, argsPoly)
  
  #lines(c(-2/3, 2/3) + xPos[[s]][d], c(med1, med1), col = siteColDarker[[s]], lwd = 1.2)
}
}
}

xPos <- c(1.75, NA, NA, 2.25) + 0.12
for(s in c(1)) {
  densities <- NULL
  
  mu1 <- modC$priors$metro[[s]]$args$mean
  sd1 <- modC$priors$metro[[s]]$args$sd
  # y1 <- seq(mu1 - 3 * sd1, mu1 + 3 * sd1, length.out = 1000)
  # pd1 <- dnorm(y1, mu1, sd1)
  # pd1 <- pd1 / max(pd1)
  # polygon(c(pd1, - rev(pd1)), c(y1, rev(y1)), border = NA, col = partsColtrans30)
  
  xWidth <- c(2/3, 1/5, 1/16)
  
  #col1 <- rgb(0,0,0,0.25)
  #col2 <- "grey55"
  
  y1 <- seq(mu1 - 3 * sd1, mu1 + 3 * sd1, length.out = 1000)
  d1 <- dnorm(y1, mu1, sd1)
  d1 <- d1 / max(d1) * 2/3
  
  x1 <- c(-d1, d1) + xPos[[s]] + 0.25
  polygon(c(x1), c(y1, rev(y1)), border = NA, col = rgb(0.65, 0.65, 0, 0.25))
  # 
  # y1 <- seq(mu1 - 1 * sd1, mu1 + 1 * sd1, length.out = 1000)
  # x1 <- c(rep(-xWidth[1], 1000), rep(xWidth[1], 1000)) + xPos[[s]] + 0.25
  # polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
  # 
  # #lines(c(-1, 1) * xWidth[1] + xPos[[s]][d], c(mu1, mu1), col = col2, lwd = 1.5, lty = 3)
  # 
  # y1 <- seq(mu1 - 2 * sd1, mu1 - 1 * sd1, length.out = 1000)
  # x1 <- c(rep(-xWidth[2], 1000), rep(xWidth[2], 1000)) + xPos[[s]] + 0.25
  # polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
  # 
  # y1 <- seq(mu1 + 2 * sd1, mu1 + 1 * sd1, length.out = 1000)
  # x1 <- c(rep(-xWidth[2], 1000), rep(xWidth[2], 1000)) + xPos[[s]] + 0.25
  # polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
  # 
  # y1 <- seq(mu1 - 3 * sd1, mu1 - 2 * sd1, length.out = 1000)
  # x1 <- c(rep(-xWidth[3], 1000), rep(xWidth[3], 1000)) + xPos[[s]] + 0.25
  # polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
  # 
  # y1 <- seq(mu1 + 3 * sd1, mu1 + 2 * sd1, length.out = 1000)
  # x1 <- c(rep(-xWidth[3], 1000), rep(xWidth[3], 1000)) + xPos[[s]] + 0.25
  # polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
  # 
}
for (s in c(1,4)) {
  h1 <- modC$alphaPosition[s]
  a1 <- ExtractParameters(resC, parameters = s, burnIn = burnIn)
  
  a1b <- as.numeric(a1[,1])
  med1 <- median(a1b)
  d1 <- density(a1b)
  d1$y <- d1$y / max(d1$y) * 2/3
  densities[[1]] <- d1
  i <- 1
  # "bean plot" type density
  # plot separate polygons for regions separated by nearly 0 density
  densThreshold <- 1 / 1000
  meetsThreshold <- densities[[i]][["y"]] >= densThreshold
  starts <- which(diff(c(FALSE, meetsThreshold)) == 1)
  ends <- which(diff(c(meetsThreshold, FALSE)) == -1) + 1
  
  for (j in seq_along(starts)) {
    argsPoly <- NULL
    argsPoly$border <- NA
    argsPoly$col <- siteColDarkerTrans2[[s]]
    indj <- starts[[j]]:ends[[j]]
    argsPoly[["y"]] <- c(densities[[i]][["x"]][indj],
                         rev(densities[[i]][["x"]][indj]))
    argsPoly[["x"]] <- c(densities[[i]][["y"]][indj] + xPos[[s]],
                         -rev(densities[[i]][["y"]][indj]) + xPos[[s]])
    do.call(polygon, argsPoly)
    
    #lines(c(-2/3, 2/3) + xPos[[s]][d], c(med1, med1), col = siteColDarker[[s]], lwd = 1.2)
  }
}
axis(4)

dev.off()

### Alpha ages


height <- 6.5
width <- 2
figName <- paste0(dir1, paste0("alphaAlignment", ".png"))
png(figName, width = width, height = height, res = 600, units = "in")
plot(0, type = "n", xlim = c(-1.5,1.5), ylim = ylim, yaxs = "i")
abline(h = seq(535, 517.5, -2.5), col = "grey80")
xPos <- c(-0.76, -0.26, 0.26, 0.76)
for(s in c(1,4)) {
    densities <- NULL
    
    mu1 <- modC$priors$metro[[s]]$args$mean
    sd1 <- modC$priors$metro[[s]]$args$sd
    # y1 <- seq(mu1 - 3 * sd1, mu1 + 3 * sd1, length.out = 1000)
    # pd1 <- dnorm(y1, mu1, sd1)
    # pd1 <- pd1 / max(pd1)
    # polygon(c(pd1, - rev(pd1)), c(y1, rev(y1)), border = NA, col = partsColtrans30)
    
    xWidth <- c(2/3, 1/5, 1/16)
    
    col1 <- rgb(0,0,0,0.25)
    col2 <- "grey55"
    y1 <- seq(mu1 - 1 * sd1, mu1 + 1 * sd1, length.out = 1000)
    x1 <- c(rep(-xWidth[1], 1000), rep(xWidth[1], 1000)) + xPos[[s]]
    polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
    
    #lines(c(-1, 1) * xWidth[1] + xPos[[s]][d], c(mu1, mu1), col = col2, lwd = 1.5, lty = 3)
    
    y1 <- seq(mu1 - 2 * sd1, mu1 - 1 * sd1, length.out = 1000)
    x1 <- c(rep(-xWidth[2], 1000), rep(xWidth[2], 1000)) + xPos[[s]]
    polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
    
    y1 <- seq(mu1 + 2 * sd1, mu1 + 1 * sd1, length.out = 1000)
    x1 <- c(rep(-xWidth[2], 1000), rep(xWidth[2], 1000)) + xPos[[s]]
    polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
    
    y1 <- seq(mu1 - 3 * sd1, mu1 - 2 * sd1, length.out = 1000)
    x1 <- c(rep(-xWidth[3], 1000), rep(xWidth[3], 1000)) + xPos[[s]]
    polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
    
    y1 <- seq(mu1 + 3 * sd1, mu1 + 2 * sd1, length.out = 1000)
    x1 <- c(rep(-xWidth[3], 1000), rep(xWidth[3], 1000)) + xPos[[s]]
    polygon(c(x1), c(y1, rev(y1)), border = NA, col = col1)
    
}




dev.off()



## AgeDepthPlot
height <- 4.1
width <- 4.3
oriPar <- par()
ages <- list()
xlims <- list()
ylims <- list()
n <- 1000
# takes several minutes:
# for (s in 1:4) {
#   ages[[s]] <- list()
#   for (b in seq_len(datC$partsAttr$nBound[s])) {
#     rangeb <- range(c(datC$partsAttr$Bound[b,s], datC$partsAttr$BoundLow[b,s]), na.rm = T)
#     heights <- seq(rangeb[1], rangeb[2], length.out = n)
#     ages[[s]][[b]] <- StratMap(resC, heights = heights, site = s, burnIn = burnIn,
#                                quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975), maxSamples = Inf)
#   }
#   xlims[[s]] <- range(c(datC$partsAttr$Bound[,s], datC$partsAttr$BoundLow[,s]), na.rm = T)
#   ylims[[s]] <- c(max(ages[[s]][[1]]$`97.5%`), min(ages[[s]][[length(ages[[s]])]]$`2.5%`))
# }
ages <- readRDS("results/cambrian/cambrian-results-ages.rds")

for (s in 1:4) {
  xlims[[s]] <- range(c(datC$partsAttr$Bound[,s], datC$partsAttr$BoundLow[,s]), na.rm = T)
  ylims[[s]] <- c(max(ages[[s]][[1]]$`97.5%`), min(ages[[s]][[length(ages[[s]])]]$`2.5%`))
}



ablines <- list(525:519,
                seq(532.5, 520, -2.5),
                521:517,
                seq(532.5, 517.5, -2.5))

yaxes <- list(526:518,
             seq(535, 515, -2.5),
             523:516,
             seq(535, 515, -2.5))

yaxesLab <- list(c(NA, NA, 524, NA,  522, NA,  520, NA, NA),
             c(NA, NA, 530, NA, 525, NA, 520, NA, NA),
             c(NA, 522, NA, 520, NA, 518, NA, NA),
             c(NA, NA, 530, NA, 525, NA, 520, NA, NA))
             
for (s in 1:4) {
figName <- paste0(dir1, paste0("figCambrianAgeDepth", s, ".png"))

png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, 1, 1), mgp = c(2.5, 0.9, 0))

  
plot(0, xlim = xlims[[s]], ylim = ylims[[s]], type = "n",
     ylab = "age (Ma)", xlab = paste("height (m) at", sitenames[s]),
     yaxt = "n")
axis(2, yaxes[[s]], yaxesLab[[s]])
abline(h = ablines[[s]], col = "grey80")
sapply(seq_len(datC$partsAttr$nBound[s]), function(b) {
  if (nrow(unique(ages[[s]][[b]])) != 1) {
    if (s == 4 & b == 1) selRow <- 1:(n-1) else selRow <- 1:n
  polygon(c(ages[[s]][[b]]$height[selRow], rev(ages[[s]][[b]]$height[selRow])),
          c(ages[[s]][[b]]$`2.5%`[selRow], rev(ages[[s]][[b]]$`97.5%`[selRow])),
          border = NA, col = partsColtrans40[datC$partsAttr$sedRatesDF[b,s]])
  # polygon(c(ages[[s]][[b]]$height[selRow], rev(ages[[s]][[b]]$height[selRow])),
  #         c(ages[[s]][[b]]$`16%`[selRow], rev(ages[[s]][[b]]$`84%`[selRow])),
  #         border = NA, col = partsColtrans30[datC$partsAttr$sedRatesDF[b,s]])
  lines(c(ages[[s]][[b]]$height[selRow]),
        c(ages[[s]][[b]]$`50%`[selRow]),
        col = partsCol[datC$partsAttr$sedRatesDF[b,s]], lwd = 3)
}
})

if (s %in% 1:3) {
  heights <- datC$tiesAttr$tiesH[,s]
  heights <- heights[!is.na(heights)]
  sapply(seq_len(datC$tiesAttr$tiesN[s]), function(d) {
    dateBounds <- datC$tiesAttr$arg1[d,s] + 2 * c(-1, 1) * datC$tiesAttr$arg2[d,s]
    points(datC$tiesAttr$tiesH[d,s], datC$tiesAttr$arg1[d,s], lwd = 1.5)
    lines(rep(datC$tiesAttr$tiesH[d,s],2), dateBounds, lwd = 1.5)
  })
}

if (s %in% c(1, 4)) {
  heights <- modC$alphaPosition[s]
  triloAge <- StratMap(resC, heights = heights, site = s, burnIn = burnIn)
  points(heights, triloAge["50%"], lwd = 1.5, pch = 4, cex = .75)
  #lines(rep(heights,2), triloAge[c("2.5%", "97.5%")], lwd = 1.5)
}
dev.off()
}

par(oriPar)

StratMapOut <- StratMapPlot(resC, sites = 1)

xlab = expression(delta^13*C)

height <- 6.5
width <- 6
figName <- paste0(dir1, paste0("figCambrianAlignment", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")

ylim <- c(534.5, 517.5)

par(mfrow = c(2,2))

StratMapPlot(resC, sites = 1, prior = F, col = "black", overridePar = F)
StratMapPlot(resC, sites = 2, prior = F, col = "black", overridePar = F)
StratMapPlot(resC, sites = 3, prior = F, col = "black", overridePar = F)
StratMapPlot(resC, sites = 4, prior = F, col = "black", overridePar = F)

         
         quantile = c(0.5), alignment = 1, ylim = ylim, separateSites = F,
     overridePar = F, ylab = "age (Ma)",
     xlab = xlab)
abline(h = seq(535, 515, -2.5), lty = 3, col = "grey")


plot(resC, quantile = c(0.5), alignment = 2, ylim = ylim, separateSites = F,
     overridePar = F, ylab = "age (Ma)",
     xlab = xlab)
abline(h = seq(535, 515, -2.5), lty = 3, col = "grey")

dev.off()

#### Appendix
parNames <- c(expression(alpha["Tiout"]), expression(alpha["Oued Sdas"]),
              expression(alpha["Talat n'Yssi"]), expression(alpha["Shukharikha"]),
              expression(scriptstyle("ln")~gamma["Tifnout strom."]), 
              expression(scriptstyle("ln")~gamma["Lie de Vin (l.)"]),
              expression(scriptstyle("ln")~gamma["Lie de Vin (m.)"]),
              expression(scriptstyle("ln")~gamma["Lie de Vin (u.)"]),
              expression(scriptstyle("ln")~gamma["Igoudine (l.)"]),
              expression(scriptstyle("ln")~gamma["Igoudine (u.)"]),
              expression(scriptstyle("ln")~gamma["Amouslek"]),
              expression(scriptstyle("ln")~gamma["Tifnout (l.)"]),
              expression(scriptstyle("ln")~gamma["Isaafen"]),
              expression(scriptstyle("ln")~gamma["Sukharikha"]),
              expression(scriptstyle("ln")~gamma["Krasnoporog"]),
              expression(scriptstyle("ln")~zeta["Oued Sdas"]),
              expression(scriptstyle("ln")~zeta["Talat n'Yssi"]),
              NA,
              expression(delta["Sukharikha"]))
  
# TracePlots
height <- 2.2
width <- 5.4
indices <- c(1:17, 19)
for (i in indices) {
figName <- paste0(dir2, paste0("traceplot",i, ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .5, .5), mgp = c(2.25, 0.8, 0))
TracePlot(resC, burnIn = burnIn, maxSamples = 1000, parameters = i, type = "o", lwd = .5,
          stratCluster = clHdbScan, overridePar = F, alpha = .7, cex = .85, ylab = "")
mtext(parNames[i], 2, 3, cex = 1.33, padj = 0.75)
dev.off()
}

# Histograms
height <- 2.5
width <- 3.4
indices <- c(1:17, 19)
colourBy1 <- "r"
for (i in indices) {
  figName <- paste0(dir2, paste0("hist",i, ".png"))
  png(figName, width = width, height = height, res = 300, units = "in")
  par(mar = c(4.1, 4.1, .5, .5), mgp = c(2.25, 0.8, 0))
  hist(resC, burnIn = burnIn, maxSamples = Inf, parameters = i, type = "o", lwd = .5,
            stratCluster = clHdbScan, overridePar = F, alpha = .6, cex = .85, xlab = "",
       colourBy = colourBy1)
  mtext(parNames[i], 1, 3, cex = 1.33, padj = -0.75)
  dev.off()
}

height <- 2.5*1.01
width <- 3.4 *1.01
figName <- paste0(dir2, paste0("hist","_post", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .68, .68), mgp = c(2.25, 0.8, 0))
hist(resC, burnIn = burnIn, maxSamples = Inf, parameters = "posterior", type = "o", lwd = .5,
     stratCluster = clHdbScan, overridePar = F, alpha = .6, cex = .85, xlab = "",
     colourBy = colourBy1, freq = F)
mtext("ln posterior", 1, 3, cex = 1, padj = -0.75)
dev.off()

# ScatterPlots
height <- 3.6
width <- 3.6
cex1 <- 0.7
alpha1 <- 0.45

figName <- paste0(dir2, paste0("scatterplot_prior_p", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))

ScatterPlot(resC, burnIn = burnIn, maxSamples = Inf, parameters = c("prior_parameters", "posterior"), type = "p", lwd = .5, asp = 1,
          stratCluster = clHdbScan, overridePar = F, alpha = alpha1, cex = cex1,  colourBy = "al",
          xlab = "ln prior parameters",
          ylab = "ln posterior")
dev.off()

figName <- paste0(dir2, paste0("scatterplot_prior_o", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))

ScatterPlot(resC, burnIn = burnIn, maxSamples = Inf, parameters = c("prior_overlap", "posterior"), type = "p", lwd = .5, asp = 1,
            stratCluster = clHdbScan, overridePar = F, alpha = alpha1, cex = cex1,  colourBy = "al",
            xlab = "ln prior overlap",
            ylab = "ln posterior")
dev.off()

figName <- paste0(dir2, paste0("scatterplot_lik_s", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))

ScatterPlot(resC, burnIn = burnIn, maxSamples = Inf, parameters = c("likelihood_spline", "posterior"), type = "p", lwd = .5, asp = 1,
            stratCluster = clHdbScan, overridePar = F, alpha = alpha1, cex = cex1,  colourBy = "al",
            xlab = "ln likelihood spline",
            ylab = "ln posterior")
dev.off()

figName <- paste0(dir2, paste0("scatterplot_lik_t", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))

ScatterPlot(resC, burnIn = burnIn, maxSamples = Inf, parameters = c("likelihood_ties", "posterior"), type = "p", lwd = .5, asp = 1,
            stratCluster = clHdbScan, overridePar = F, alpha = alpha1, cex = cex1,  colourBy = "al",
            xlab = "ln likelihood ties",
            ylab = "ln posterior")
dev.off()

figName <- paste0(dir2, paste0("scatterplot_4_15", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))

ScatterPlot(resC, burnIn = burnIn, maxSamples = Inf, parameters = c(4, 15), type = "p", lwd = .5,
            stratCluster = clHdbScan, overridePar = F, alpha = alpha1, cex = cex1,  colourBy = "al",
            xlab = "",
            ylab = "")
mtext(parNames[4], 1, 3, cex = 1.33, padj = -0.5)
mtext(parNames[15], 2, 3, cex = 1.33, padj = 0.5)

dev.off()


figName <- paste0(dir2, paste0("scatterplot_1_11", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))

ScatterPlot(resC, burnIn = burnIn, maxSamples = Inf, parameters = c(1,11), type = "p", lwd = .5, 
            stratCluster = clHdbScan, overridePar = F, alpha = alpha1, cex = cex1,  colourBy = "al",
            xlab = "",
            ylab = "")
mtext(parNames[1], 1, 3, cex = 1.33, padj = -0.5)
mtext(parNames[11], 2, 3, cex = 1.33, padj = 0.5)
dev.off()

## Clustering plot
cols <- c(rgb(0,0,0,0.45), hcl.colors(3, "Dark 3", alpha = 0.45))

s <- 4
h4 <- modC$heightsForClustering[[s]][[2]]
a4 <- StratMap(resC, h4, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")

s <- 3
h3 <- modC$heightsForClustering[[s]][[6]]
a3 <- StratMap(resC, h3, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")

s <- 2
h2 <- modC$heightsForClustering[[s]][[6]]
a2 <- StratMap(resC, h2, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")

s <- 1
h1 <- modC$heightsForClustering[[s]][[7]]
a1 <- StratMap(resC, h1, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")


hist(as.numeric(a1[1,2:48001]))

figName <- paste0(dir2, paste0("scatterplot_tops_1", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))

dat1 <- cbind(as.numeric(a1[1,2:48001]), as.numeric(a2[1,2:48001]), clHdbScan$cluster)
dat1 <- dat1[sample(1:nrow(dat1), nrow(dat1)),]
plot(dat1[,1], dat1[,2], bg = cols[1 + dat1[,3]], col = NA,
     pch = c(24,21:23)[1+dat1[,3]],
     xlab = "age (Ma) at top of Tiout",
     ylab = "age (Ma) at top of Oued Sdas")
dev.off()

figName <- paste0(dir2, paste0("scatterplot_tops_2", ".png"))
png(figName, width = width, height = height, res = 300, units = "in")
par(mar = c(4.1, 4.1, .6, .6), mgp = c(2.25, 0.8, 0))
dat1 <- cbind(as.numeric(a3[1,2:48001]), as.numeric(a4[1,2:48001]), clHdbScan$cluster)
dat1 <- dat1[sample(1:nrow(dat1), nrow(dat1)),]
plot(dat1[,1], dat1[,2], bg = cols[1 + dat1[,3]], col = NA,
     pch = c(24,21:23)[1+dat1[,3]],
         xlab = "age (Ma) at top of Talat n'Yssi",
ylab = "age (Ma) at top of Shukharika")
dev.off()


### histograms
maxSamples <- Inf
colourBy <- "run"
parameter <- 1
cols <- "Viridis"
breaks = 100
for (parameter in 1:18) {
  
  parameter <- 19
hist(resC, burnIn = burnIn, overridePar = F, prior = T, freq = F,
     maxSamples = maxSamples, colourBy = colourBy, parameters = parameter,
     breaks = breaks, colourScale = "Viridis", alpha = 2/3)
  
  
}


hist(resC, burnIn = burnIn, overridePar = F, prior = T, freq = F,
     maxSamples = maxSamples, colourBy = colourBy, parameters = parameter)

par(new = TRUE)
hist(resC, burnIn = burnIn, overridePar = F, parameters = 2, xaxt = "n", yaxt = "n", prior = F, freq = F)




## alignment plot Appendix


xlab = expression(delta^13*C)
xlim <- range(datC$signal$d13C) + diff(range(datC$signal$d13C)) * c(-0.015, 0.03)

height <- 6.5
width <- 2.9
al <- 0
quant1 <- 0.5
ptCex <- 2/3
figName <- paste0(dir2, paste0("figCambrianAlignment_appendix",al,"_", quant1, ".png"))
png(figName, width = width, height = height, res = 300, units = "in")

ylim <- c(535, 516.75)

par(mfrow = c(1,1))

plot(resC, quantile = quant1, alignment = al, ylim = ylim, separateSites = F,
     overridePar = F, ylab = "age (Ma)",
     xlab = xlab, cex = ptCex, col = siteColTrans, type = "n", burnIn = burnIn,
     xlim = xlim, maxSamples = Inf,
     xaxt = "n", yaxs = "i", stratCluster = clHdbScan)

axis(1, seq(-8, 8, 2), c(NA, NA, -4, NA, 0, NA, 4, NA, 8))

abline(v = c(-4, 0, 4), col = "grey80")
abline(h = seq(532.5, 515, -2.5), lty = 1, col =  "grey80")

par(new = TRUE)

plot(resC, quantile = quant1, alignment = al, ylim = ylim, separateSites = F,
     overridePar = F, ylab = "age (Ma)",
     xlab = xlab, cex = ptCex, col = siteColTrans, burnIn = burnIn,
     xlim = xlim, maxSamples = Inf,
     xaxt = "n", yaxs = "i", stratCluster = clHdbScan)

dev.off()



### 2D density plots 

library(MASS)

lab1 <- "ln prior parameters"
lab2 <- "ln prior overlap"
par12 <- c("likelihood_spline", "likelihood_ties")

lab1 <- "ln prior overlap"
lab2 <- "ln likelihood"

par12 <- c("prior_overlap", "likelihood")
samples <- list()

samples <- list()
for(i in 1:4) {
  samples[[i]] <- ExtractParameters(
    resC, 
    parameters = par12, 
    burnIn      = burnIn, 
    stratCluster= clHdbScan, 
    alignment   = i - 1,
    maxSamples = Inf
  )
}

asp = 1

asp = 0


lab1 <- "age (Ma) top of Tiout"
lab2 <- "age (Ma) top of Oued Sdas"
s <- 1
h1 <- modC$heightsForClustering[[s]][[7]]
a1 <- StratMap(resC, h1, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")

s <- 2
h2 <- modC$heightsForClustering[[s]][[6]]
a2 <- StratMap(resC, h2, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")

for (i in 1:4) {
  samples[[i]] <- cbind(as.numeric(a1[1,(2:48001)[clHdbScan$cluster == i-1]]), as.numeric(a2[1,(2:48001)[clHdbScan$cluster == i-1]]))
}

lab1 <- "age (Ma) at top of Talat n'Yssi"
lab2 <- "age (Ma) at top of Sukharikha"

s <- 3
h1 <- modC$heightsForClustering[[s]][[6]]
a1 <- StratMap(resC, h1, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")

s <- 4
h2 <- modC$heightsForClustering[[s]][[2]]
a2 <- StratMap(resC, h2, s, burnIn = burnIn, maxSamples = Inf, quantiles = "samples")

for (i in 1:4) {
  samples[[i]] <- cbind(as.numeric(a1[1,(2:48001)[clHdbScan$cluster == i-1]]), as.numeric(a2[1,(2:48001)[clHdbScan$cluster == i-1]]))
}


samplesAll <- do.call(rbind, samples)
xlim <- range(samplesAll[,1])
ylim <- range(samplesAll[,2])

bandwidth <- c(
  bw.nrd0(samplesAll[,1]),
  bw.nrd0(samplesAll[,2])
) * 10

# Larger grid => smoother density
nGrid <- 300

#    Compute kde2d for each group (same bandwidth & grid), then
#    weight the density by the group's sample size.

kdeListWeighted <- list()
groupSizes      <- sapply(samples, nrow)

for (i in seq_along(samples)) {
  # Raw kernel density
  kdeObj <- MASS::kde2d(
    x = samples[[i]][,1],
    y = samples[[i]][,2],
    h = bandwidth,
    n = nGrid,
    lims = c(xlim, ylim)  # ensure the same bounding box for all 
    
  )
  # Multiply the z-matrix by the group’s sample size
  kdeObj$z <- kdeObj$z  * groupSizes[i]^(1/3)
  
  kdeListWeighted[[i]] = kdeObj
}

#    Global min & max from weighted densities

zAll <- unlist(lapply(kdeListWeighted, function(kde) kde$z))
zMin <- min(zAll)   # probably 0 or near 0
zMax <- max(zAll)

# Number of contour levels
nLevels <- 5
levels <- seq(zMin, zMax, length.out = nLevels + 1)
levels[1] <- levels[2]/10
# Midpoints for alpha scaling
levelMids <- (levels[-1] + levels[-length(levels)]) / 2

#    Global alpha scaling function
scaleAlphaGlobal <- function(zValue) {
  (zValue - zMin) / (zMax - zMin)  # [zMin..zMax] => [0..1]
}

#    Base colors
colors1 <- c(hcl.colors(9, "Viridis"))[c(9,1,3,7)]

# Plot 
width <- 3.25
height <- 3.25
figName <- paste0(dir1, paste0("figC-", lab1, lab2, "-b.png"))
png(figName, width = width, height = height, res = 400, units = "in")
par(mgp = c(2.1, .75, 0), mar = c(4.1, 4.1, .75, .75))
plot(NA, NA, xlim = xlim, ylim = ylim, asp = asp,
     xlab = lab1,
     ylab = lab2,
     xaxs = "i", yaxs = "i")

# Overlay each group's weighted density, masking low values

for (i in c(1,4,3,2)) {
  
  baseCol <- colors1[i]
  
  # Build alpha ramp for each contour band
  alphaVals <- 
    scaleAlphaGlobal(levelMids)^(if (i == 1) 5/12 else 1/2 )
  fillCols <- rgb(
    red   = col2rgb(baseCol)[1]/255,
    green = col2rgb(baseCol)[2]/255,
    blue  = col2rgb(baseCol)[3]/255,
    alpha = alphaVals
  )
  
  # Weighted z-matrix
  zMat <- kdeListWeighted[[i]]$z
  
  # Mask out near-zero densities so empty regions remain white
  zMaxGroup <- max(zMat)
  threshold <- 0.01 * zMaxGroup
  zMat[zMat < threshold] <- NA
  
  # Overplot filled contours
  .filled.contour(
    x      = kdeListWeighted[[i]]$x,
    y      = kdeListWeighted[[i]]$y,
    z      = zMat,
    levels = levels,
    col    = fillCols
  )
  # contour(
  #   x      = kdeListWeighted[[i]]$x,
  #   y      = kdeListWeighted[[i]]$y,
  #   z      = zMat,
  #   levels = levels,
  #   col    = fillCols,
  #   drawlabels = T,
  #   add = T,
  #   lwd = 2,
  # )
}

dev.off()
