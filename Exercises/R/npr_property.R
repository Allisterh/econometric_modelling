#============================================================================
#
#    Sampling properties of the nonparametric kernel regression 
#    estimator
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")

#
#---------------  Properties of the Nadaraya-Watson Estimator ---------------
# 

# Simulate the model
t  <- 500
ut <- 0.1*rnorm(t)                # N(0,0.1^2) 
xt <- -2 + 4*runif(t)              # U(-2,2) 
mx <- 0.3*exp( -4*(xt + 1)^2 ) + 0.7*exp( -16*(xt - 1)^2 )  
yt  <- mx + ut

# Construct true conditional mean
x   <- seq(-2,2,0.01)
mxt <- 0.3*exp( -4*(x + 1)^2 ) + 0.7*exp( -16*(x - 1)^2 )  


#  Compute kernel regression for h = 0.5  
h  <- 0.2
fx <- rep(0, length(x))
fxy <- rep(0, length(x))
for (i in seq(x)) {
  z      <- ((x[i] - xt)/h)    
  fx[i]  <- mean( dnorm(z)/h )
  fxy[i] <- mean( dnorm(z)*yt/h )  
}
mx1 <- fxy / fx

# Compute kernel regression for h = 0.1 
h  <- 0.02
for (i in seq(x) ) {
  z      <- ((x[i] - xt)/h)    
  fx[i]  <- mean( dnorm(z)/h )
  fxy[i] <- mean( dnorm(z)*yt/h ) 
}
mx2 <- fxy/ fx

#**********************************************************************
#***
#***     Generate graphs
#***
#**********************************************************************

figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

#--------------------------------------------------------#
# Panel (a)
plot(x,mx1, type="l", lty=2,
     main = "(a) Larger Bandwidth",
     xlab = expression(x[t]),
     ylab = expression(paste("y"[t]*", ",   m(x[t]))),
     ylim = c(-0.1, 0.8),
     bty = "l")
lines(x,mxt, lty=1)


#--------------------------------------------------------#
# Panel (b)
plot(x,mx2, type="l",
     main = "(b) Smaller Bandwidth",
     xlab = expression(x[t]),
     ylab = expression(paste("y"[t]*", ",   m(x[t]))),
     ylim = c(-0.1, 0.8),
     bty = "l")
lines(x,mxt, lty=2)
