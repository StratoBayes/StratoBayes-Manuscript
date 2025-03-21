# Conceptual figures with simulated data
#
# load package
library(StratoBayes)
# read data with StratData()
sData <- StratData(signal = "data/simulated/signalData01.csv")
# visualise sData
plot(sData)
#
# template for priors and model run
StratModelTemplate(stratData = sData, 
                   alignmentScale = "height",
                   sedModel = "site",
                   alphaPosition = "top")
#
# filling the priors template
## using min and max of site1 to delineate the range where alpha2 (middle of 
## site2) may be:
minAlpha <- min(SignalLookUp(sData, type = "height", site = "site1"))
maxAlpha <- max(SignalLookUp(sData, type = "height", site = "site1")) + 
  diff(range(SignalLookUp(sData, type = "height", site = "site2")))
priors <- structure(list(
  "alpha_site2" = UniformPrior(min = minAlpha, max = maxAlpha),
  "gammaLog_site2" = NormalPrior(mean = 0, sd = log(2))),
  class = c("StratPrior", "list"))
# create model object 
model <- StratModel(stratData = sData,
                    priors = priors,
                    alignmentScale = "height",
                    sedModel = "site",
                    alphaPosition = "top",
                    nKnots = 10)
# run model
result <- RunStratModel(stratObject = sData,
                        stratModel = model,
                        nRun = 1,
                        nIter = 1000)
#
# directory for saving images
dir1 <- "figures/conceptual-figure/"
# figures
width <- 2.1
height <- 4
parOp <- par()
parCustom <- parOp
parCustom$mgp <- c(2.25, 0.8, 0)
xlim <- range(sData$signal$value)
ylim <- range(sData$signal$height)
mainCex <- 1

figName <- paste0(dir1, "fig1-section1.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(sData, sites = 1, xlim = xlim, ylim = ylim, cex.main = mainCex)
dev.off()

figName <- paste0(dir1, "fig1-section2.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(sData, sites = 2, xlim = xlim, ylim = ylim, cex.main = mainCex)
dev.off()

figName <- paste0(dir1, "fig1-aligned-section1.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(result, sites = 1, xlim = xlim, ylim = ylim, main = "site1", overridePar = F)
dev.off()

figName <- paste0(dir1, "fig1-aligned-section2.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(result, sites = 2, xlim = xlim, ylim = ylim, main = "site2", overridePar = F)
dev.off()

width <- 2.33
height <- 2.33
figName <- paste0(dir1, "fig1-prior1.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(priors, parameters = 1, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = "",
     xlim = c(minAlpha - 2, maxAlpha + 2), ylim = c(0, 0.053))
mtext(expression(bold(alpha)), side = 1, line = 2, cex = 1.3)
#mtext("    site2", side = 1, line = 2)
dev.off()

figName <- paste0(dir1, "fig1-prior2.png")
png(figName, width = width, height = height, res = 300, units = "in")
suppressWarnings(par(parCustom))
plot(priors, parameters = 2, col = rgb(0,0,0,0.25), overridePar = FALSE, xlab = "",
     log = TRUE, xaxt = "n")
axis(1, at = c(-2, -1, 0, 1, 2), labels = c(NA, -1, 0, 1, NA))
mtext(expression(""~~~bold(gamma)), side = 1, line = 2, cex = 1.3)
#mtext("    site2", side = 1, line = 2)
dev.off()
