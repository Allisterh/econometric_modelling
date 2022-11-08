#============================================================================
#
#      Program to simulate a garch model 
#
#============================================================================
  
rm(list = ls(all=T))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
# ------------------------ Helper Function ----------------------------------
#
# Load required functions - trimr, figure
source("EMTSUtil.R")

#
# ------------------------ Garch Model Properties ---------------------------
#
garch_simulate <- function( )
{
  nobs <- 100000
  mlag <- 20
  
  z  <- rnorm( nobs+1000)
  u  <- z                         
  y  <- z
  h  <- rep(1, nobs+1000)
  
  # Choose parameter values  
  #a0 <- 0.1 a1 <- 0.7 b1 <- 0.2
  #a0 <- 0.05 a1 <- 0.15 b1 <- 0.8
  a0 <- 0.05 
  a1 <- 0.05 
  b1 <- 0.9
  
  h[1] <- a0/(1-a1-b1)
  u[1] <- rnorm(1, 0,sqrt(h[1]))
  
  # Generate data 
  pb <- txtProgressBar(min=0, max=nobs+1000, style=3)
  for (t in 2:nobs+1000) {  
    h[t] <- a0 + a1*u[t-1]^2 + b1*h[t-1]    # Conditional variance  
    u[t] <- z[t]*sqrt( h[t] )               # Disturbance term     
    y[t] <- u[t]                            # Returns  
    setTxtProgressBar(pb, t)
  }
  close(pb)
  y <- trimr( y,1000,0 )  
  y <- y - mean(y)
  
  
  #load test
  # Compute the autocorrelation of returns squared
  acfy1 <- acf(y^2, mlag, plot=F)$acf
  
  
  # Parameter values second model
  a0 <- 0.05
  a1 <- 0.15
  b1 <- 0.80
  
  
  # Generate data 
  pb <- txtProgressBar(min=0, max=nobs+1000, style=3)
  for (t in 2:nobs+1000) {  
    h[t] <- a0 + a1*u[t-1]^2 + b1*h[t-1]    # Conditional variance  
    u[t] <- z[t]*sqrt( h[t] )               # Disturbance term     
    y[t] <- u[t]                            # Returns  
    setTxtProgressBar(pb, t)
  }
  close(pb)
  y <- trimr( y,1000,0 )  
  y <- y - mean(y)
  
  
  #load test
  # Compute the autocorrelation of returns squared
  acfy2 <- acf(y^2, mlag, plot=F)$acf
  
  #*********************************************************************
  #**     Generate graphs
  #*********************************************************************
  figure()
  par(mfrow=c(2, 1), xaxs="i", yaxs="i", mar=c(5,5,5,5))  
  lags <- seq(mlag+1)
  
  matplot(lags,cbind(acfy1, rep(0, mlag+1)), type="l",
          main = "(a)",
          xlab = "p",
          xlim = c(0, 10),
          ylim = c(-0.5, 1),
          ylab = expression(ACF(y1[t]^2)),
          bty= "l")

  matplot(lags,cbind(acfy2, rep(0, mlag+1)), type="l",
          main = "(b)",
          xlab = "p",
          xlim = c(0, 10),
          ylim = c(-0.5, 1),
          ylab = expression(ACF(y2[t]^2)),
          bty= "l")
}


