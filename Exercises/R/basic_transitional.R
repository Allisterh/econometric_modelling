#=========================================================================
#
#   Maximum likelihood estimation of the transitional distribution of the 
#   Vasciek model of interest rates using Ait Sahalia's (1996) data.
#
#=========================================================================
rm(list = ls(all=TRUE))
graphics.off()

# load any compatibility requirements
source("EMTSUtil.R")

# Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
load("eurodata.RData")
rt <- eurodata[,4] * 100

# Regress r(t) on r(t-1)
y  <- rt[-1]
x  <- rt[-length(rt)]

theta <- as.matrix( lm(y~x)$coef )
ones <- array(1, c(length(y), 1))
x <- cbind(ones, x)

e     <- y - x %*% theta
sig2  <- mean( e^2 )

cat( '\nMLE of alpha:      ', theta[1], '\n')
cat( '\nMLE of rho:        ', theta[2], '\n')
cat( '\nMLE of beta:       ', theta[2] - 1, '\n')
cat( '\nMLE of sigma^2:    ', sig2, '\n')


cat( '\n' )

minrt <- min(rt)
maxrt <- max(rt)
medianrt <- median(rt)

cat( '\nMinimum interest rate:      ', minrt, '\n' )
cat( '\nMedian interest rate:       ', medianrt, '\n')
cat( '\nMaximum interest rate:      ', maxrt, '\n')

cat( '\n' )

cat( '\nMLE of mean based on the transitional distribution:   ', -theta[1]/(theta[2]-1), '\n') 
cat( '\nMLE of variance based on the transitional distribution:   ', -sig2/((theta[2]-1)*(2+theta[2]-1)), '\n')

# Compute transitional densities
r <- seq(0, 30, 0.1)

mu_min <- theta[1] + theta[2] * minrt
fmin   <- dnorm( ( r - mu_min )/sqrt(sig2) )/sqrt(sig2)
 
mu_med <- theta[1] + theta[2] * medianrt
fmed   <- dnorm( ( r - mu_med )/sqrt(sig2) )/sqrt(sig2)

mu_max <- theta[1] + theta[2] * maxrt
fmax   <- dnorm( ( r - mu_max )/sqrt(sig2) )/sqrt(sig2)

#***     Generate graph ***#

figure()
plot(r, fmin, type="l", lty=2,
     xlab="r",
     ylab="f(r)", bty="l")
lines(r, fmed, lty=1)
lines(r, fmax, lty=3)

# Add a legend to the plot  
legend("topright",                       
       legend=c("Minimum", "Median", "Maximum"),             
       lty=c(2,1,3),                    
       lwd=c(1,1,1))
