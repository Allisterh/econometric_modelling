#==========================================================================
#
#   Simulating a regression model with autocorrelation 
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
# load required functions - figure
source("EMTSUtil.R")

#
#--------------------------- Simulations  -----------------------------------
#

auto_simulate <- function(){
  t <- 200  
  # Model parameters     
  beta0  <- 2
  beta1  <- 1
  rho1   <- 0.95
  delta1 <- 0.95
  sigma  <- 3
  
  # Generate the exogenous variable
  x <- 0.5*seq(t+100) + rnorm(t+100)             
  
  # Simulate a regression model with an AR(1) disturbance term      
  v <- sigma*rnorm(t+100)
  zeros <- rep(0, t+100)
  u <- zeros                                                                           
  y <- zeros                                                                          
  
  for (i in 2:t+100) {
    u[i] <- rho1*u[i-1] + v[i]                                                                    
    y[i] <- beta0 + beta1*x[i] + u[i]  
  }
  y_ar1 <- y
  
  # Simulate a regression model with a MA(1) disturbance term      
  v <- sigma*rnorm(t+100)                                                                  
  u <- zeros                                                                             
  y <- zeros                                                                             
  
  for (i in 2:t+100) {
      u[i] <- v[i] + delta1*v[i-1]                                                                 
      y[i] <- beta0 + beta1*x[i] + u[i]                                                             
  }
  y_ma1 <- y
  
  # Trim data to overcome startup problems      
  y_ar1 <- y_ar1[-100]
  y_ma1 <- y_ma1[-100]
  x     <- x[-100]                      
  mu    <- beta0 + beta1*x                                               
  
  #**********************************************************************
  #***
  #***     Plot the series
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
  
  # Panel (a)
  plot(x,y_ar1, pch=19,
       main = "(a) AR(1) Regression Model",
       xlab = expression(x[t]),
       ylab = expression(y[t]),
       xlim = c(40, 160),
       ylim = c(20, 160),
       bty = "l")
  lines(x,mu)
  
  # Panel (b)
  plot(x,y_ma1, pch=19,
       main = "(b) MA(1) Regression Model",
       xlab = expression(x[t]),
       ylab = expression(y[t]),
       xlim = c(40, 160),
       ylim = c(20, 160),
       bty = "l")
  lines(x,mu)
  
  
  
  # Simulate a regression model with an AR(2) disturbance term      
  t <- 200                                                                                          
  beta0 <- 2
  beta1 <- 1
  rho1  <- 0.1
  rho2  <- -0.9
  sigma <- 3
  
  x <- 0.5*seq(t+100) + rnorm(t+100)
  v <- sigma*rnorm(t+100)                                                               
  
  u0 <- c(0,0)
  rho <- c(rho1,rho2)
  for (i in 3:length(v+1)) {
     u0[i] <- v[(i-1)] + rho[1]*u0[(i-1)] + rho[2]*u0[(i-2)]  
  }
  u     <- u0
  y     <- beta0 + beta1*x + u
  y_ar2 <- y[-100]
  
  # Simulate a regression model with an ARMA(2,2) disturbance term      
  t <- 200                                                                                           
  beta0  <- 2
  beta1  <- 1
  rho1   <- 0.1
  rho2   <- -0.9
  delta1 <- 0.3
  delta2 <- 0.2
  sigma  <- 3
  
  x <- 0.5*seq(t+100) + rnorm(t+100)        	# xt is generated from a trend with normal additive errors    
  v <- sigma*rnorm(t+100)                       # vt is N(0,sigma^2)                                          
  
  vend <- length(v)
  tmpvd1d2 <- v[3:vend] + delta1*v[2:(vend-1)] + delta2*v[1:(vend-2)]
  tmp <- c(0, 0, tmpvd1d2)
  u0 <- c(0,0)
  rho <- c(rho1,rho2)
  for (i in 3:length(tmp+1)) {
    u0[i] <- tmp[(i-1)] + rho[1]*u0[(i-1)] + rho[2]*u0[(i-2)]  
  }
  u      <- u0
  y      <- beta0 + beta1*x + u
  y_arma <- y[-100]
  
}

