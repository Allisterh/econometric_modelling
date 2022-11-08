#=========================================================================
#
#   Program to model the number of strikes per annum in the U.S.
#   Data are from Kennan (1985)
#
#=========================================================================

# clear all
rm(list = ls(all = TRUE))
graphics.off()

# load any compatibility requirements
source("EMTSUtil.R")

# Number of strike per annum
yt <- c(8, 6, 11, 3, 3, 2, 18, 2,  9) 


# Maximum likelihood estimate of the number of strikes 
theta <- mean(yt)
cat('\nMLE of mean number of strikes (in years) =', theta, '\n')

# Plot the estimated distribution of strike numbers  
y <- seq(0, max(yt), 1)                          
f <- dpois(y,theta)

figure()
plot(y, f, type="l",
     col="darkblue",
     main = "Distribution of Number of strikes",
     ylab = expression(paste("f(y;", theta, ")")),
     xlab = "Number of Strikes",
     bty="l")

figure()
par(xaxs="i", yaxs="i")
hist(yt, breaks=11, 
     col="darkblue",
     main = "Histogram of Number of strikes",     
     xlab = "Number of Strikes")

