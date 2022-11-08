#============================================================================
#
#     Demonstration of bias and variance of nonparametric estimators 
#     using different bandwidths. 
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")

#
#-----------------------  Sampling Properties -------------------------------
#


# Simulate the data from a normal distribution 
t   <- 500        
mu  <- 0
sig <- 1
yt  <- mu + sig*rnorm(t)

# Generate the population (normal) distribution 
y      <- seq(-5, 5, 0.1)
f_norm <- dnorm( ( y - mu )/sig )/sig

# Estimate the nonparametric density for alternative bandwidths
f1 <- rep(0, length(y))
h  <- 1

for (i in seq(y) ) {
  z     <- ((y[i] - yt)/h)    
  f1[i] <- mean( dnorm(z)/h ) 
}

f2 <- rep(0, length(y) )
h <- 0.1
for (i in seq(y)) {
  z    <- ((y[i] - yt)/h)    
  f2[i] <- mean( dnorm(z)/h )
  
}       
          

#****************************************************************************
#
#     Generate graphs
#
#****************************************************************************


figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))


#--------------------------------------------------------#
# Panel (a)
plot(y,f_norm,type="l",
     main = "(a) Larger Bandwidth",
     xlab = expression(y),
     ylab = expression(f(y[t])),
     ylim = c(0, 0.5),
     bty="l")

lines(y,f1, lty=2, col="blue")

#--------------------------------------------------------#
# Panel (b)
plot(y,f_norm,type="l",
     main = "(b) Smaller Bandwidth",
     xlab = expression(y),
     ylab = expression(f(y[t])),
     ylim = c(0, 0.5),
     bty="l")

lines(y,f2, lty=2, col="blue")

