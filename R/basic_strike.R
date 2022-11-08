#=========================================================================
#
#   Program to model the duration of strikes per annum in the U.S.
#   Data are from Kennan (1985)
#
#=========================================================================

rm(list = ls(all = TRUE))
graphics.off()

# load any compatibility requirements
source("EMTSUtil.R")

# Read the data: US strike data from 1968 to 1976
load("strike.RData")
yt <- as.matrix(strike.number)               

# Maximum likelihood estimate of strike duration 
theta <- mean(yt)

cat( '\nMLE of mean strike duration (in days) = ', theta, '\n')


# Plot the estimated distribution of strike durations  
y <- seq(0,300,10)
f <- (1/theta)*exp(-y/theta)

figure()
plot(y,f, type="l",
     col="darkblue",
     main = "Distribution of duration of strikes", 
     ylab = expression(paste("f(y;", theta, ")")),
     xlab = "Duration of Strikes",
     bty="l")

# Plot the histogram 
figure()
par(xaxs="i", yaxs="i")
hist(yt,breaks=21,
     col="darkblue",
     main = "Histogram of the Duration of strikes",    
     xlab = "Duration of Strikes")
