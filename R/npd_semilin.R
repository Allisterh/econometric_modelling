#=========================================================================
#
#   Semi-parametric Look Ahead Estimator (LAE): 
#   linear autoregressive model
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(124, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")


#
#-----------------  Semi-Parametric Density Estimation ----------------------
#
# Parameters
mu    <- 1.0
phi   <- 0.5
sig2  <- 1.0
t     <- 5

# Grid
y <- seq(-4, 6, 0.1)
n <- length(y)

#  Simulate AR(1) model   
yt <- mu*rep(1,t)

for (i in 2:t) {
   yt[i] <- mu + phi*(yt[i-1]) + sqrt(sig2)*rnorm(1)  
}

# Generate the stationary distribution using the LAE 
m <- mu + phi*yt
s <- sqrt(sig2)
    
f <- array(0, c(length(yt),length(y)) )
for (j in seq(yt) ) {   
   for (i in seq(y)) {
     z    <- (y[i] - m[j])/s
     f[j,i]<- (1/sqrt(2*pi*s^2))*exp(-0.5*(z^2))     
   }     
}  
f_lae <- colMeans(f)

# Generate the true stationary distribution using the analytical expression      
m <- mu/(1 - phi)
s <- sqrt( sig2/(1 - phi^2))
z <- ( y - m )/ s 

f_true <- (1/sqrt(2*pi*s^2))*exp(-0.5*(z^2))  


# Kernel density estimate of the stationary distribution    
h     <- 1.06*sd(yt)*t^(-1/5)  
    
fx <- rep(0, length(y) )
for (i in seq(y)){
  z      <- ((y[i] - yt)/h)    
  fx[i]  <- mean( dnorm(z)/h )
}

#**********************************************************************
#***
#***     Generate graph
#***
#**********************************************************************

figure()
plot(y,f_lae,type="l", lty=2, col="blue",
     main = "True, LAE and Kernel Densities",
     ylab = expression(f(y)),
     xlab = expression(y),
     ylim = c(0, 0.8),
     bty="l")
lines(y,f_true, lty=1, col="black")
lines(y,fx, lty=4, col="green") 
legend("topright",                       
      legend=c("LAE Density",
                "True Density",
                "Kernel Density"),             
      lty=c(2, 1, 4),                    
      lwd=c(1,1,1,1),
      col = c("blue", "black", "green"))




