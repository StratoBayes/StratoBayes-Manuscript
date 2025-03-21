### results for the manuscript text
##
## Model illustration
library(StratoBayes)
result2a <- readRDS("results/model-illustration/alignment-2-peaks.rds")
burnIN <- 10000
clusterNew <- Cluster(result2a, burnIn = burnIN, clusterMethod = "hdbscan", 
                      maxSamples = 10000, minPts = 100)
ScatterPlot(result2a, colourBy = "a", stratCluster = clusterNew, burnIn = burnIN, maxSamples = 10000)
result2a$summary <- summary(result2a, burnIn = burnIN, stratCluster = clusterNew)
# results to display in text:
period1 <- which(vapply(1:2, function(x) abs(result2a$summary[[x]]$summary["alpha_site2", "mean"] - 2 * pi) < 1, logical(1)))
period2 <- which(c(1,2) != period1)
probValidation1AlignPeriod1 <- format(round(100*result2a$summary[[period1]]$proportionAmongAllSamples, 1),nsmall = 1)
probValidation1AlignPeriod2 <- format(round(100*result2a$summary[[period2]]$proportionAmongAllSamples, 1),nsmall = 1)

probValidation1AlignPeriod1Runs <- vapply(1:3, function(x) {
  clusterNewX <- Cluster(result2a, burnIn = burnIN, clusterMethod = "hdbscan", 
                        maxSamples = 10000, minPts = 100, runs = x)
  sumR <- summary(result2a, stratCluster = clusterNewX, runs = x)
  period1 <- which(vapply(1:2, function(x) abs(sumR[[x]]$summary["alpha_site2", "mean"] - 2 * pi) < 1, logical(1)))
  sumR[[period1]]$proportionAmongAllSamples
}, numeric(1))

probValidation1AlignPeriod1SEM <- format(round(100*sd(probValidation1AlignPeriod1Runs)/sqrt(3), 0),nsmall = 0)
probValidation1AlignPeriod2SEM <- format(round(100*sd(1 - probValidation1AlignPeriod1Runs)/sqrt(3), 0),nsmall = 0)

# summarise gamma across both alignments
s2 <- summary(result2a, alignment = "none", burnIn = burnIN)

gammaValidation1_median <- format(round(exp(s2[[1]]$summary$`50%`[2]), 2), nsmall = 2)
gammaValidation1_025 <- format(round(exp(s2[[1]]$summary$`2.5%`[2]), 2), nsmall = 2)
gammaValidation1_975 <- format(round(exp(s2[[1]]$summary$`97.5%`[2]), 2), nsmall = 2)
gammaValidation1 <- c(gammaValidation1_median, gammaValidation1_025, gammaValidation1_975)

  
validation1 <- list(probAlignPeriod1 = probValidation1AlignPeriod1, 
                    probAlignPeriod2 = probValidation1AlignPeriod2,
                    probAlignSEM = probValidation1AlignPeriod1SEM,
                    nIter = formatC(result2a$mcmc$nIter, format = "f", big.mark = ",", digits = 0),
                    nChains = formatC(result2a$mcmc$nChains, format = "f", big.mark = ",", digits = 0),
                    burnIn = formatC(burnIN, format = "f", big.mark = ",", digits = 0),
                    nSamples = 6000,
                    nThin = 25,
                    gamma = gammaValidation1)
saveRDS(validation1, "results/text-results/validation1.rds")



##
## Cambrian study
nIter <- 750000
burnIn <- 150000
nThin <- 50
nRun <- 4
nIterAfterBurnIn <- nIter - burnIn
nSamplesAfterBurnIn <- nIterAfterBurnIn / nThin * nRun
resC <- readRDS("results/cambrian/cambrian-results.rds")
#clHdbScan <- Cluster(resC, clusterMethod = "hdbscan", burnIn = burnIn, 
#                     maxSamples = Inf, minPts = nSamplesAfterBurnIn / 100)
#saveRDS(clHdbScan, "results/cambrian/cambrian-results-clusterHdbscan.rds")
clHdbScan <- readRDS("results/cambrian/cambrian-results-clusterHdbscan.rds")

#summaryFull <- summary(resC, burnIn = burnIn, maxSamples = Inf,
#                       stratCluster = clHdbScan)
#saveRDS(summaryFull, "results/cambrian/cambrian-results-summaryFull.rds")
summaryFull <- readRDS("results/cambrian/cambrian-results-summaryFull.rds")

summaryTotal <- summary(resC, alignment = "none", burnIn = burnIn, maxSamples = Inf)
ESS <- format(round(summaryTotal$alignmentAll$multiESS, 0), nsmall = 0)
nSamples <- format(round(summaryTotal$general$nSamplesUsed, 0), nsmall = 0) 
indices <- c(1:17, 19)
psrf <- format(round(range(summaryTotal$general$psrf[indices]), 2), nsmall = 0) 
al1 <- format(round(100*summaryFull$alignment1$proportionAmongAllSamples, 0),nsmall = 0)
al2 <- format(round(100*summaryFull$alignment2$proportionAmongAllSamples, 1),nsmall = 0)
al3 <- format(round(100*summaryFull$alignment3$proportionAmongAllSamples, 1),nsmall = 0)
unaligned <- format(round(100*summaryFull$alignment0$proportionAmongAllSamples, 1),nsmall = 0)

trilo_tiout <- format(round(summaryTotal$alignmentAll$summary[1,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)
trilo_sib <- format(round(summaryTotal$alignmentAll$summary[4,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)

trilo_tiout_sib_diff <- format(round(quantile(ExtractParameters(resC, 4, burnIn = burnIn)[,1] - ExtractParameters(resC, 1, burnIn = burnIn)[,1], probs = c(0.5, 0.025, 0.975)), 2), nsmall = 2)
  
trilo_tiout_1 <- format(round(summaryFull$alignment1$summary[1,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)
trilo_sib_1 <- format(round(summaryFull$alignment1$summary[4,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)

trilo_tiout_sib_diff1 <- format(round(quantile(ExtractParameters(resC, 4, burnIn = burnIn, stratCluster = clHdbScan, alignment = 1)[,1] - ExtractParameters(resC, 1, burnIn = burnIn, stratCluster = clHdbScan, alignment = 1)[,1], probs = c(0.5, 0.025, 0.975)), 2), nsmall = 2)

trilo_tiout_2 <- format(round(summaryFull$alignment2$summary[1,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)
trilo_sib_2 <- format(round(summaryFull$alignment2$summary[4,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)

trilo_tiout_sib_diff2 <- format(round(quantile(ExtractParameters(resC, 4, burnIn = burnIn, stratCluster = clHdbScan, alignment = 2)[,1] - ExtractParameters(resC, 1, burnIn = burnIn, stratCluster = clHdbScan, alignment = 2)[,1], probs = c(0.5, 0.025, 0.975)), 2), nsmall = 2)

trilo_tiout_3 <- format(round(summaryFull$alignment3$summary[1,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)
trilo_sib_3 <- format(round(summaryFull$alignment3$summary[4,c("50%", "2.5%", "97.5%")], 2),nsmall = 2)

trilo_tiout_sib_diff3 <- format(round(quantile(ExtractParameters(resC, 4, burnIn = burnIn, stratCluster = clHdbScan, alignment = 3)[,1] - ExtractParameters(resC, 1, burnIn = burnIn, stratCluster = clHdbScan, alignment = 3)[,1], probs = c(0.5, 0.025, 0.975)), 2), nsmall = 2)

trilo_tiout_0_range <- format(round(range(ExtractParameters(resC, 1, burnIn = burnIn, stratCluster = clHdbScan, alignment = 0)), 2),nsmall = 2)
trilo_sib_0_range <- format(round(range(ExtractParameters(resC, 4, burnIn = burnIn, stratCluster = clHdbScan, alignment = 0)), 2),nsmall = 2)

trilo_tiout_sib_diff3 <- format(round(quantile(ExtractParameters(resC, 4, burnIn = burnIn, stratCluster = clHdbScan, alignment = 3)[,1] - ExtractParameters(resC, 1, burnIn = burnIn, stratCluster = clHdbScan, alignment = 3)[,1], probs = c(0.5, 0.025, 0.975)), 2), nsmall = 2)

t1 <- ExtractIterations(resC, burnIn = burnIn, stratCluster = clHdbScan, alignment = "all")
sapply(t1, length)
t2 <- ExtractIterations(resC, burnIn = burnIn, stratCluster = clHdbScan, alignment = 1)
sapply(t2, length)
t3 <- ExtractIterations(resC, burnIn = burnIn, stratCluster = clHdbScan, alignment = 3)
sapply(t3, length)

siteMultiplierTalat1 <- format(round(exp(summaryFull$alignment1$summary[17,c("50%", "2.5%", "97.5%")]), 2), nsmall = 2)
siteMultiplierTalat2 <- format(round(exp(summaryFull$alignment2$summary[17,c("50%", "2.5%", "97.5%")]), 2), nsmall = 2)
siteMultiplierTalat3 <- format(round(exp(summaryFull$alignment3$summary[17,c("50%", "2.5%", "97.5%")]), 2), nsmall = 2)

siteMultiplierOued1 <- format(round(exp(summaryFull$alignment1$summary[16,c("50%", "2.5%", "97.5%")]), 2), nsmall = 2)
siteMultiplierOued2 <- format(round(exp(summaryFull$alignment2$summary[16,c("50%", "2.5%", "97.5%")]), 2), nsmall = 2)
siteMultiplierOued3 <- format(round(exp(summaryFull$alignment3$summary[16,c("50%", "2.5%", "97.5%")]), 2), nsmall = 2)

sedAmouslek1 <- sapply(
  exp(summaryFull$alignment1$summary["gammaLog_Amouslek", c("50%", "2.5%", "97.5%")]), 
  function(x) {
    # Round to three significant digits
    roundedValue <- signif(x, 3)
    
    # Format numbers based on their size
    if (roundedValue >= 10000) {
      format(roundedValue, nsmall = 0, big.mark = ",", scientific = FALSE)
    } else {
      as.character(round(roundedValue, 0)) # Convert smaller values to integers
    }
  }
)

AmouslekBottomAgeSamples <- StratMap(resC, site = 3, heights = resC$data$partsAttr$BoundLow[5,3],
                                     quantiles = "samples", stratCluster = clHdbScan, alignment = 1,
                                     maxSamples = 4000)
AmouslekTopAgeSamples <- StratMap(resC, site = 3, heights = resC$data$partsAttr$Bound[5,3],
                                     quantiles = "samples", stratCluster = clHdbScan, alignment = 1,
                                  maxSamples = 4000)
AmouslekHeight <- format(round(resC$data$partsAttr$Bound[5,3] - resC$data$partsAttr$BoundLow[5,3],0), nsmall = 0)
  
AmouslekDuration1 <- c(format(round(quantile(as.numeric(AmouslekBottomAgeSamples[1,2:4001]) - as.numeric(AmouslekTopAgeSamples[1,2:4001]), probs = c(0.5, 0.975)) * 1000,0),nsmall = 0),
                       format(round(quantile(as.numeric(AmouslekBottomAgeSamples[1,2:4001]) - as.numeric(AmouslekTopAgeSamples[1,2:4001]), probs = c(0.025)) * 1000,1),nsmall = 1))[c(1,3,2)]
  
  
talatDateAge1 <- format(round(StratMap(resC, site = 3, heights = resC$data$tiesAttr$tiesH[1,3], stratCluster = clHdbScan, alignment = 1, burnIn = burnIn)[c("50%", "2.5%", "97.5%")], 1), nsmall = 1)

TioutDate23diffHeight <- diff(resC$data$tiesAttr$tiesH[2:3,1])
n <- 10^7
TioutDate23diffAge <- rnorm(n,resC$data$tiesAttr$arg1[2,1], resC$data$tiesAttr$arg2[2,1]) -
      rnorm(n,resC$data$tiesAttr$arg1[3,1], resC$data$tiesAttr$arg2[3,1])

dateAmouslekSedRate <- quantile(TioutDate23diffHeight / TioutDate23diffAge, probs = c(0.5, 0.025, 0.975))
dateAmouslekSedRate <- sapply(1:3, function(x) format(round(dateAmouslekSedRate[x], c(0,1,0)[x])))

AgeIIpeak <- StratMap(resC, resC$data$signal$height[resC$data$signalAttr$siteInd[[2]]][542], site = 2, stratCluster = clHdbScan, burnIn = burnIn, alignment = 1)

sigmaFixed <- format(round(resC$model$sigmaFixed, 2), nsmall = 2)

cambrian <- list(ESS = ESS,
                 nSamples = nSamples,
                 psrf = psrf,
                 alignmentProp = c(al1, al2, al3, unaligned),
                 sigmaFixed = sigmaFixed,
                 trilo_tiout = trilo_tiout,
                 trilo_tiout_1 = trilo_tiout_1,
                 trilo_sib = trilo_sib, 
                 trilo_tiout_sib_diff = trilo_tiout_sib_diff,
                 trilo_sib1 = trilo_sib_1,
                 trilo_sib2 = trilo_sib_2,
                 trilo_sib3 = trilo_sib_3,
                 trilo_tiout_0_range = trilo_tiout_0_range,
                 trilo_tiout_sib_diff1 = trilo_tiout_sib_diff1,
                 trilo_tiout_sib_diff2 = trilo_tiout_sib_diff2,
                 trilo_tiout_sib_diff3 = trilo_tiout_sib_diff3,
                 siteMultiplierOued1 = siteMultiplierOued1,
                 siteMultiplierOued2 = siteMultiplierOued2,
                 siteMultiplierOued3 = siteMultiplierOued3,
                 siteMultiplierTalat1 = siteMultiplierTalat1,
                 siteMultiplierTalat2 = siteMultiplierTalat2,
                 siteMultiplierTalat3 = siteMultiplierTalat3,
                 talatDateAge1 = talatDateAge1,
                 sedAmouslek1 = sedAmouslek1,
                 AmouslekHeight = AmouslekHeight,
                 AmouslekDuration1 = AmouslekDuration1,
                 dateAmouslekSedRate = dateAmouslekSedRate
                 
                 )
saveRDS(cambrian, "results/text-results/cambrian.rds")
