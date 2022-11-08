#=========================================================================
#
#   Maximum likelihood estimation of the stationary distribution of the 
#   Vasciek model of interest rates using Ait Sahalia's (1996) data.
#
#=========================================================================

rm(list = ls(all = TRUE))
graphics.off()

# load any compatibility requirements
source("EMTSUtil.R")

# Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
load("eurodata.RData")
rt <- eurodata[,4] * 100


# Maximum likelihood estimates of the stationary distribution 
mu_r   <- mean( rt )
sig2_r <- mean( (rt - mu_r)^2 )

cat( '\nMLE of mean of stationary distribution:      ', mu_r, '\n')
cat('MLE of variance of stationary distribution: ', sig2_r, '\n')


# Compute stationary density from -5% to 25%
r     <- seq(-5, 25, 0.1)
fnorm <- dnorm( ( r - mu_r )/sqrt(sig2_r) )/sqrt(sig2_r)
prob  <- pnorm( ( 0 - mu_r )/sqrt(sig2_r) )

cat( '\nProbability of a negative interest rate:     ', prob, '\n')

#***     Generate graph    ***#
figure()

plot(r,fnorm, type="l",     
     col="darkblue",   
     main="Stationary Distribution for Vasicek Model",
     ylab = "f(r)",
     xlab = "Interest Rate", 
     bty="l")
