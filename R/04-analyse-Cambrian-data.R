# Case study with Cambrian d13C data
#
# load package
library(StratoBayes)
#
## read data
# paths:
cambrianSignal <- "data/cambrian/cambrian_signal.csv"
cambrianTies <- "data/cambrian/cambrian_ties.csv"
cambrianParts <- "data/cambrian/cambrian_parts.csv"

# read raw data for modification
signalRaw <- read.csv(cambrianSignal)
partsRaw <- read.csv(cambrianParts)
tiesRaw <- read.csv(cambrianTies)

# remove lower parts in MS7 and Siberian section to make model estimation easier
partsRaw <- partsRaw[-which(partsRaw$formation == "Sukharikha (lower)"),]
partsRaw <- partsRaw[-1,]
signalRaw <- signalRaw[-which(signalRaw$site == "MS7" & signalRaw$height < 350),]
signalRaw <- signalRaw[-which(signalRaw$site == "Siberia" & signalRaw$height < 350),]

# read and format StratData from raw data
cData <- StratData(signal = signalRaw, ties = tiesRaw, 
                   parts = partsRaw, signalColumn = "d13C",
                   partitionColumn = "lithology_for_model_final",
                   referenceSite = "MS6", 
                   selectSites = c("MS7", "MS6", "MS3", "Siberia"))
plot(cData, alpha = 1, cex = 0.8)
saveRDS(cData, "results/cambrian/cambrian-data.rds")

## additional data: astronomical sedimentation rate estimates
# min (cycle more + 100.5 kyr SE duration)
minSedRate <- c(284.5771144, 
                222.2222222,
                194.4125526,
                253.968254,
                224.1424457)    
# max sed rate (cycle less + 92.5 kyr duration)
maxSedRate <- c(515.3153153,
                310.4247104,
                249.6314496,
                340.8585056,
                272.1780604)
# Average sed rate (m/Myr) - mean value for Normal distribution
meanSedRate <- c(370.4663212,
                 260.3626943,
                 219.343696,
                 292.3370603,
                 246.4018423)   # Desired mean 
# generate mean and sigma
zLower <- qnorm(0.025)
zUpper <- qnorm(0.975)
sigmaAstro <- (log(maxSedRate) - log(minSedRate)) / (zUpper - zLower)
muAstro <- log(minSedRate) - zLower * sigmaAstro


## generate model template
StratModelTemplate(cData, alignmentScale = "age", sedModel = "s x p")

## determine alpha positions
# use dates 
date1MS7height <- cData$ties[1,"height"]
date1MS3height <- cData$ties[5,"height"]
# priors on date ages
MS7tie1Mean <- cData$ties[1, "mean"]
MS7tie1SD <- cData$ties[1, "sd"]

MS3tieMean <- cData$ties[5, "mean"]
MS3tieSD <- cData$ties[5, "sd"]

# use trilobite appearances
tioutAlpha <- 1150 #cData$ties[8, "height"]
siberianAlpha <- 584 # height
# prior on trilobite appearance age
triloMean <- 520 # mean age (Ma)
triloSD <- 2 # sd (Myr)

# list alpha positions (heights)
alphaPositions <- c(tioutAlpha, date1MS7height, date1MS3height,
                    siberianAlpha)

## calculate plausible sedimentation rate priors
ageDifferencesMS7 <- abs(diff(cData$ties[1:4, "mean"]))  
heightDifferencesMS7 <- diff(cData$ties[1:4, "height"])
ageDifferencesMS6 <- abs(diff(cData$ties[6:8, "mean"]))  
heightDifferencesMS6 <- diff(cData$ties[6:8, "height"])

sedRatesMS7 <- heightDifferencesMS7 / ageDifferencesMS7
sedRatesMS6 <- heightDifferencesMS6 / ageDifferencesMS6
sedRateMor <- mean(log(c(sedRatesMS7, sedRatesMS6)))

# broad priors (for less constrained partitions)
Morocco_prior <- list(mu = sedRateMor, sigma = 3/4)
Siberia_prior <- list(mu = 3.3, sigma = 3/4)

## create priors object
priors <- structure(list(
  "alpha_MS6" = NormalPrior(mean = triloMean, sd = triloSD),
  "alpha_MS7" = NormalPrior(mean = MS7tie1Mean,  sd = triloSD),
  "alpha_MS3" = NormalPrior(mean = MS3tieMean, sd = triloSD),
  "alpha_Siberia" = NormalPrior(mean = triloMean, sd = triloSD),
  "gammaLog_Tifnout strom." = NormalPrior(mean = Morocco_prior$mu, sd = Morocco_prior$sigma),
  "gammaLog_Lie-de-vin lower" = NormalPrior(mean = muAstro[5], sd = sigmaAstro[5]),
  "gammaLog_Lie-de-vin middle" = NormalPrior(mean = muAstro[4], sd = sigmaAstro[4]),
  "gammaLog_Lie-de-vin upper" = NormalPrior(mean = muAstro[3], sd = sigmaAstro[3]),
  "gammaLog_Igoudine lower" = NormalPrior(mean = muAstro[2], sd = sigmaAstro[2]),
  "gammaLog_Igoudine upper" = NormalPrior(mean = muAstro[1], sd = sigmaAstro[1]),
  "gammaLog_Amouslek" = NormalPrior(mean = Morocco_prior$mu, sd = Morocco_prior$sigma),
  "gammaLog_Tifnout lower" = NormalPrior(mean = Morocco_prior$mu, sd = Morocco_prior$sigma),
  "gammaLog_Isaafen" = NormalPrior(mean = Morocco_prior$mu, sd = Morocco_prior$sigma),
  "gammaLog_Sukharikha" = NormalPrior(mean =Siberia_prior$mu, sd = Siberia_prior$sigma),
  "gammaLog_Krasnoporog" = NormalPrior(mean = Siberia_prior$mu, sd = Siberia_prior$sigma),
  "zetaLog_MS7" = NormalPrior(mean = 0, sd = 1/4),
  "zetaLog_MS3" = NormalPrior(mean = 0, sd = 1/4),
  "zetaLog_Siberia" = Fixed(0),
  "gap_Siberia_1" = ExponentialPrior(rate = 1)),
  class = c("StratPrior", "list"))

## create model
modelMorSib <- StratModel(stratData = cData,
                          priors = priors,
                          alignmentScale = "age",
                          sedModel = "s x p",
                          alphaPosition = alphaPositions,
                          nKnots = 40)
# saveRDS(modelMorSib, "results/cambrian/cambrian-model.rds")

## MCMC settings
# replace with own directory
modelDir <- tempdir() #"D://OneDrive - Durham University/Strat/alignment-tests/run-2024-20-12-sxp4"
# dir.create(modelDir)
modelName <- "model1"
nThin <- 50
nChains <- 24
nIterAdapt <- 150000
burnIn <- nIterAdapt
nIterAfterBurnIn <- 600000
totalIter <- burnIn + nIterAfterBurnIn
seed <- 1:4
nRun <- length(seed)
# collate settings
cambrianSettings <- list(
  nThin = nThin,
  nChains = nChains,
  nIterAdapt = as.character(format(nIterAdapt, big.mark = ",", scientific = FALSE)),
  burnIn = as.character(format(burnIn, big.mark = ",", scientific = FALSE)) ,
  nIterAfterBurnIn = as.character(format(nIterAfterBurnIn, big.mark = ",", scientific = FALSE)) ,
  totalIter = as.character(format(totalIter, big.mark = ",", scientific = FALSE)) ,
  nRun = nRun,
  nSamplesPerRun = as.character(format(nIterAfterBurnIn / nThin, big.mark = ",", scientific = FALSE)), 
  nSamples = as.character(format(nRun * nIterAfterBurnIn / nThin, big.mark = ",", scientific = FALSE)),
  
  MoroccoPriorMu = format(round(Morocco_prior[[1]], 2), nsmall = 2),
  MoroccoPriorSd = format(round(Morocco_prior[[2]], 2), nsmall = 2),
  SiberiaPriorMu = format(round(Siberia_prior[[1]], 2), nsmall = 2),
  SiberiaPriorSd = format(round(Siberia_prior[[2]], 2), nsmall = 2),
  
  MoroccoPrior025 = format(round(qlnorm(c(0.025), Morocco_prior$mu, Morocco_prior$sigma), 1), nsmall = 1),
  MoroccoPrior975 = format(round(qlnorm(c(0.975), Morocco_prior$mu, Morocco_prior$sigma), 0), nsmall = 0),
  SiberiaPrior025 = format(round(qlnorm(c(0.025), Siberia_prior$mu, Siberia_prior$sigma), 2), nsmall = 2),
  SiberiaPrior975 = format(round(qlnorm(c(0.975), Siberia_prior$mu, Siberia_prior$sigma), 1), nsmall = 1)
)

saveRDS(cambrianSettings,
        "results/text-results/cambrianSettings.rds")

## model run:
# using blocks of 50,000 iterations to avoid memory issues
# starting with three burn-in blocks (150,000 iterations total)
# careful, each block takes several hours to run!

res01 <- RunStratModel(stratObject = cData,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       endAdapting = 50000,
                       adaptProposalInterval = 200,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res01, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-1-cyclo.rds")

res02 <- RunStratModel(stratObject = res01,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       endAdapting = 50000,
                       adaptProposalInterval = 200,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res02, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-2-cyclo.rds")


res03 <- RunStratModel(stratObject = res02,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       endAdapting = 50000,
                       adaptProposalInterval = 200,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res03, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-20-cyclo.rds")

# 12 blocks after burn-in (600,000 iterations total)

res04 <- RunStratModel(stratObject = res03,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res04, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-21-cyclo.rds")

res05 <- RunStratModel(stratObject = res04,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res05, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-22-cyclo.rds")

res06 <- RunStratModel(stratObject = res05,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res06, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-23-cyclo.rds")


res07 <- RunStratModel(stratObject = res06,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res07, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-24-cyclo.rds")


res08 <- RunStratModel(stratObject = res07,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res08, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-25-cyclo.rds")

res09 <- RunStratModel(stratObject = res08,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res09, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-26-cyclo.rds")


res10 <- RunStratModel(stratObject = res09,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res10, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-27-cyclo.rds")

res11 <- RunStratModel(stratObject = res10,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res11, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-28-cyclo.rds")

res12 <- RunStratModel(stratObject = res11,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res12, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-29-cyclo.rds")



res13 <- RunStratModel(stratObject = res12,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res13, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-30-cyclo.rds")


res14 <- RunStratModel(stratObject = res13,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)
saveRDS(res14, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-31-cyclo.rds")

res15 <- RunStratModel(stratObject = res14,
                       stratModel = modelMorSib,
                       modelDir = modelDir,
                       modelName = modelName,
                       nChains = nChains,
                       runParallel = T,
                       nIter = 50000, # 
                       nThin = 50, 
                       visualise = list(500, 1:2, c(3, 5, 7)),
                       adapt = F,
                       seed = 1:4,
                       saveMCMC = TRUE,
                       clustWithParametersMCMC = TRUE,
                       parametersInitial = NULL
)

saveRDS(res15, "D://OneDrive - Durham University/Strat/alignment-tests/resultMorSibTioutLithShort1220-32-cyclo.rds")

