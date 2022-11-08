#============================================================================
#
#   Semi-parametric Look Ahead Estimator (LAE): 
#   nonlinear threshold autoregressive model
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

# load required library - ones, zeros, meshgrid, reshape
library("matlab")

#
#-----------------  Semi-Parametric Density Estimation ----------------------
#

# Parameters
theta <- 0.6
t     <- 5
# Simulate tar(1) model   
yt <- theta*rep(1,t)

for (i in 2:t) {
  yt[i] <- theta*abs(yt[i-1]) + sqrt(1 - theta^2)*rnorm(1,1)  
}

# or set values to those in text

yt <- c(0.600, -0.279, -1.334, 0.847, 0.894)

# Grid
y <- seq(-4,4,0.1)
n <- length(y)

# Generate the stationary distribution using the LAE 
s <- sqrt(1 - theta^2)
    
f <- array(0, c(length(yt),length(y)) )
for (j in seq(yt)) {
  for (i in seq(y) ) {
    z    <- (y[i] - theta*yt[j])/s
    f[j,i]<- (1/sqrt(2*pi*s^2))*exp(-0.5*(z^2))   
  }
}
f_lae <- colMeans(f)



# Generate the true stationary distribution using the analytical expression      
delta  <- theta/sqrt(1 - theta^2)
f_true <- 2*dnorm(y)*pnorm(delta*y)    

# # Kernel density estimate of the stationary distribution    **/  
h     <- 1.06*sd(yt)*t^(-1/5)  
    
fx <- rep(0, length(y))
for (i in seq(y)) {
  z      <- ((y[i] - yt)/h)    
  fx[i]  <- mean( dnorm(z)/h )  
}

#**********************************************************************
#***
#***     Generate graphs
#***
#**********************************************************************
figure()
par(mfrow=c(1,2))

#--------------------------------------------------------#
# Panel (a)
plot(y,f[1,],type="l", lty=2, col="blue",
     main = "(a) Component Transitional Densities",
     xlab = expression(y),
     ylab = expression(f(y)),
      ylim = c(0, 0.5),
     bty="l")
lines(y,f[3,], lty=1)
lines(y,f[5,], lty=6, col="darkblue")


#--------------------------------------------------------#
# Panel (b)
plot(y,f_lae,type="l", lty=2, col="blue",
     main = "(b) True, LAE and Kernel",
     xlab = expression(y),
     ylab = expression(f(y)),
     ylim = c(0, 0.5),
     bty="l")
lines(y,f_true,lty=1)
lines(y,fx, lty=6, col="darkblue")

