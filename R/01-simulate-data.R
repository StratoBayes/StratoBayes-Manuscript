################################################################################
# Create simulated data
################################################################################
#
# load package
library(StratoBayes)
#
#
################################################################################
#
# 1 - create simple data set for demonstration
#
################################################################################
#
# sine function to generate underlying signal y for any heights x
GenerateSine <- function(x) {
sin(x) - 1/5*(x - 2 * 2 * pi)
}
#
# set random seed to ensure consistency of results
set.seed(0)
# reference section data
# number of periods
periods <- 3
# number of data points for reference section
n1 <- 70
x1 <- runif(n1, 0, 2 * pi * periods)
y1 <- GenerateSine(x1) + rnorm(n1, 0, 1/4)
# second section data
n2 <- 30
x2 <- runif(n2, 2 * pi, 4 * pi)
y2 <- GenerateSine(x2) + rnorm(n2, 0, 1/3)
#
# visualise signal
plot(x1, y1, bg = rgb(.8, 0, 0), pch = 21)
points(x2, y2, bg = rgb(0,.6,.8), pch = 22)
#
# scale x2 (change sedimentation rate)
sedRate2 <- 3
x2 <- x2 * sedRate2
# shift x2 to start at 0
x2 <- x2 - min(x2)
#
# format data for reading it with StratData()
signalData <- data.frame(site = rep(c("site1", "site2"), c(n1, n2)),
                         height = c(x1, x2),
                         value = c(y1, y2))
# save data
write.csv(signalData, "data/simulated/signalData01.csv")
#
#
#
################################################################################
#
# 2 - create data set to demonstrate multiple alignments
#
################################################################################
#
# sine function to generate underlying signal y for any heights x

GenerateSine <- function(x, shift, sd, dampening) {
  y <- rep(NA, length(x))
  ind1 <- which(x <= (2 + shift)*pi)
  ind2 <- which(x <= (4 + shift)*pi & x > (2 + shift)*pi)
  ind3 <- which(x > (4 + shift)*pi )

  y[ind1] <- dampening[1] * sin(x[ind1])  + rnorm(length(ind1), 0, sd) - (1 - dampening[1])
  y[ind2] <- dampening[2] * sin(x[ind2])  + rnorm(length(ind2), 0, sd) - (1 - dampening[2])
  y[ind3] <- dampening[3] * sin(x[ind3])  + rnorm(length(ind3), 0, sd) - (1 - dampening[3])
 # return
  y
}

GenerateSineConsistentNoise <- function(x, shift, sd, dampening, nPerCycle) {
  y <- rep(NA, length(x))
  period <- 2 * pi  # The fundamental period of sine
  
  # generate a single period's worth of noise at fixed points
  refX <- seq(0, period, length.out = nPerCycle)  # Reference phase points
  refNoise <- rnorm(nPerCycle, 0, sd)  # Generate stable noise per phase point
  
  # consistent noise based on phase
  getNoise <- function(vals) {
    phase <- (vals + shift * pi) %% period  # phase position
    index <- findInterval(phase, refX, all.inside = TRUE)  # closest reference noise
    refNoise[index]  # stable periodic noise
  }
  
  # apply to different dampening regions
  ind1 <- which(x <= (2 + shift) * pi)
  ind2 <- which(x <= (4 + shift) * pi & x > (2 + shift) * pi)
  ind3a <- which(x <= (5 + shift) * pi & x > (4 + shift) * pi)
  ind3b <- which(x > (5 + shift) * pi)
  ind3 <- c(ind3a, ind3b)
  
  y[ind1] <- dampening[1] * sin(x[ind1]) + getNoise(x[ind1]) + (1 - dampening[1])
  y[ind2] <- dampening[2] * sin(x[ind2]) + getNoise(x[ind2]) + (1 - dampening[2])
  y[ind3a] <- dampening[2] * sin(x[ind3a]) + getNoise(x[ind3a]) + (1 - dampening[2])
  y[ind3b] <- dampening[3] * sin(x[ind3b]) + getNoise(x[ind3b]) + (1 - dampening[3])
  
  return(y)
}


set.seed(0)
# sd
sd <- 1/5
# number of data points for reference section
periods <- 3
nPerCycle <- 250
cycles <- 1/2 + periods
n1 <- cycles * nPerCycle
shift <- -0.5
x1 <- seq((-0.5 + shift)*pi, pi * (2 * periods + shift + 0.5), length.out = n1)
y1 <- GenerateSineConsistentNoise(x1, shift, sd, c(1, 1, 0.75), nPerCycle)

set.seed(0)

n2a <- 1* nPerCycle
x2a <- x1[x1 >= shift * pi & x1 <= (shift + 2) * pi]
y2a <- GenerateSineConsistentNoise(x2a, shift, sd, c(1, NA, NA), nPerCycle)
# shift and stretch
x2a <- 2*x2a - min(2*x2a) 

set.seed(0)

n2b <- 0.5 * nPerCycle
x2b <- x1[x1 >= shift * pi & x1 <= (shift + 2) * pi]
y2b <- GenerateSineConsistentNoise(x2b, shift, sd, c(0.75, NA, NA), nPerCycle) - 1
# shift and stretch
x2b <- 2*x2b - min(2*x2b) 

n2b <- 1* nPerCycle
x2b <- x1[x1 >= shift * pi & x1 <= (shift + 2) * pi]
y2b <- rnorm(n2b, 0, 1/2 * sd)
# shift and stretch
x2b <- 2*x2b - min(2*x2b) 

# to move the start of the period to 0
x1 <- x1 - shift*pi

zpart1 <- x1 <  pi
z1 <- rnorm(n1, 0, sd) + (zpart1 * 1/5 * (x1 - pi))
z2 <- rnorm(n2a, 0, sd)

nz2 <- ceiling(0.5 * n1)


x1c <- x1[seq(1, n1, 2)]

plot(y1, x1)
points(y2a, x2a, bg = "dodgerblue", col = NA, pch = 21)
points(z1, x1, bg = "orange", col = NA, pch = 21)
points(z2, x2a, bg = "red", col = NA, pch = 21)

signalData_2a <- data.frame(site = rep(c("site1", "site2"), c(n1, n2a)),
                            height = c(x1, x2a),
                            value = c(y1, y2a))

# save signalData_2a
write.csv(signalData_2a, "data/simulated/model-illustration.csv")

plot(StratData(signalData_2a))



# format signalData_2b
signalData_2b <- data.frame(site = rep(c("site1", "site2"), c(n1, n2b)),
                            height = c(x1, x2a),
                            signal1 = c(y1, y2a),
                            signal2 = c(z1, z2))

# save signalData_2b
write.csv(signalData_2b, "data/simulated/signalData02b.csv")

plot(StratData(signalData_2b, signalColumn = c("signal1", "signal2")), signal = 1:2)
