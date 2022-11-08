#============================================================================
#
#   Program to simulate bilinear time series models
#
#============================================================================
  
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - figure
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Function to simulate a bilinear model
#----------------------------------------------------------------------------
bilinear <- function(t,phi,theta,gam)
{
  y <- rep(0, t+100)
  u <- rnorm(t+100)
  
  for (i in 2:(t+100)) {
    y[i] <- phi*y[i-1] + u[i] + theta*u[i-1] + gam*y[i-1]*u[i-1]
  }
  
  y <- y[101:length(y)]
  return(y)
}

#
#------------------------------- Bilinearity  -------------------------------
#


nlm_blinear <- function()
{
  # Parameters
  t <- 200
  
  # Simulate the models
  # Gamma = 0.0
  y1 <- bilinear(t,0.4,0.2,0.0)       
  
  # Gamma = 0.4
  y2 <- bilinear(t,0.4,0.2,0.4)       
  
  # Gamma = 0.8
  y3 <- bilinear(t,0.4,0.2,0.8)          
  
  # Gamma = 1.2
  y4 = bilinear(t,0.4,0.2,1.2) 
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(mfrow=c(2,2))
  
  s = seq(t)
  
  #--------------------------------------------------------#
  # Panel (a)
  plot(s, y1, type="l", 
       main=expression(paste("(a) ", gamma, " = 0.0")),
       xlab = "t",
       ylab = "y",
       bty = "l")
  
  #--------------------------------------------------------#
  # Panel (b)
  plot(s, y2, type="l", 
       main=expression(paste("(b) ", gamma, " = 0.4")),
       xlab = "t",
       ylab = "y",
       bty = "l")
  
  #--------------------------------------------------------#
  # Panel (c)
  plot(s, y3, type="l", 
       main=expression(paste("(c) ", gamma, " = 0.8")),
       xlab = "t",
       ylab = "y",
       bty = "l")
  
  #--------------------------------------------------------#
  # Panel (d)
  plot(s, y4, type="l", 
       main=expression(paste("(d) ", gamma, " = 1.2")),
       xlab = "t",
       ylab = "y",
       bty = "l")
}
