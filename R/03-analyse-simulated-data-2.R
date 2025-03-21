# Results with simulated data 02
#
# load package
library(StratoBayes)
# read data with StratData()
sData2a <- StratData(signal = "data/simulated/model-illustration.csv")
# visualise sData
plot(sData2a)
#
set.seed(0)
# template for priors and model run
StratModelTemplate(stratData = sData2a, 
                   alignmentScale = "height",
                   sedModel = "site",
                   alphaPosition = "top")
#
# filling the priors template
## using min and max of site1 to delineate the range where alpha2 (middle of 
## site2) may be:
minAlpha <- min(SignalLookUp(sData2a, type = "height", site = "site1"))
maxAlpha <- max(SignalLookUp(sData2a, type = "height", site = "site1"))
priors <- structure(list(
  "alpha_site2" = UniformPrior(min = 0 -  pi, max = 8 * pi),
  "gammaLog_site2" = NormalPrior(mean = 0, sd = 1)),
  class = c("StratPrior", "list"))
#
# create model object 
# using large individualSplineIter for sigmaFixed
model <- StratModel(stratData = sData2a,
                    priors = priors,
                    alignmentScale = "height",
                    sedModel = "site",
                    alphaPosition = "top",
                    nKnots = 20, 
                    sigmaFixed = TRUE,
                    individualSplineIter = 20000)

# run model
result2a <- RunStratModel(stratObject = sData2a,
                          stratModel = model,
                          nRun = 3,
                          runParallel = TRUE,
                          nIter = 60000,
                          endAdapting = 10000,
                          nThin = 25,
                          nChains = 16,
                          seed = 1:3)
# save results
#  saveRDS(result2a, "D://OneDrive - Durham University/Strat/alignment-tests/alignment-2-peaks.rds")
# result2a <- readRDS("D://OneDrive - Durham University/Strat/alignment-tests/alignment-2-peaks.rds")
saveRDS(result2a, "results/model-illustration/alignment-2-peaks.rds")

## additional tests (not in MS)
# no overlap penalty
# create model object 
set.seed(0)
model <- StratModel(stratData = sData2a,
                    priors = priors,
                    alignmentScale = "height",
                    sedModel = "site",
                    alphaPosition = "top",
                    nKnots = 20, 
                    sigmaFixed = TRUE, 
                    overlapPenalty = FALSE,
                    minSignalOverlap = sData2a$signalAttr$N[2])
# run model
result2a2 <- RunStratModel(stratData = sData2a,
                          stratModel = model,
                          nRun = 3,
                          runParallel = T,
                          nIter = 10000,
                          nThin = 5,
                          seed = 1:3)
# save results
saveRDS(result2a2, "results/model-illustration/alignment-2-peaks-noPenalty.rds")
# switch off overlapPenalty  / on 
TracePlot(result2a2, parameters = c(1,2,"posterior"), colourBy = "alignment", type = "p")


### second signal
#
# read data with StratData()
sData2b <- StratData("data/simulated/signalData02b.csv", signalColumn = c("signal1", "signal2"))
# visualise sData
plot(sData2b, signal = 1:2)
#
set.seed(0)
# template for priors and model run
StratModelTemplate(stratData = sData2b, 
                   alignmentScale = "height",
                   sedModel = "site",
                   alphaPosition = "top")
#
# filling the priors template
## using min and max of site1 to delineate the range where alpha2 (middle of 
## site2) may be:
minAlpha <- min(SignalLookUp(sData2b, type = "height", site = "site1"))
maxAlpha <- max(SignalLookUp(sData2b, type = "height", site = "site1"))
priors <- structure(list(
  "alpha_site2" = UniformPrior(min = minAlpha - 1/2 * pi, max = maxAlpha + 2 * pi),
  "gammaLog_site2" = UniformPrior(min = log(1/4), max = log(4))),
  class = c("StratPrior", "list"))
#
# create model object 
model <- StratModel(stratData = sData2b,
                    priors = priors,
                    alignmentScale = "height",
                    sedModel = "site",
                    alphaPosition = "top",
                    nKnots = 20, 
                    sigmaFixed = TRUE,
                    individualSplineIter = 20000)
# run model
result2b <- RunStratModel(stratObject = sData2b,
                          stratModel = model,
                          nRun = 3,
                          runParallel = TRUE,
                          nIter = 25000,
                          endAdapting = 5000,
                          nThin = 10,
                          nChains = 16,
                          seed = 1:3)

saveRDS(result2b, "results/model-illustration/alignment-second_signal.rds")
# switch off overlapPenalty  / on 
TracePlot(result2b, parameters = c(1,2,"posterior"), colourBy = "alignment", type = "p")
#
#
#
