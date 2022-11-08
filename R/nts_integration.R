#==========================================================================
#
#     Simulate stochastic trends with different orders of integration
#
#==========================================================================

rm(list = ls(all=T))
graphics.off()
set.seed(145, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, recserar, seqa
source("EMTSUtil.R")

#
#------------------------- Integration --------------------------------------
#

nts_integration <- function() {
  # Generate simulate variables
  t     <- 200
  delta <- 0.1
  
  y_noise <- delta + rnorm(t)                                #      White noise
  y_i0    <- recserar(cbind(delta + rnorm(t)),cbind(0.0),cbind(0.5))              #      Stationary process (I(0))
  y_i1    <- recserar(cbind(delta + rnorm(t)),cbind(0.0),cbind(1.0))              #      Nonstationary process (I(1))
  y_i2    <- recserar(cbind(delta + rnorm(t)),rbind(0,0),rbind(2.0,-1.0))#      Nonstationary process (I(2))
  
  
  #**********************************************************************
  #***
  #***     Generate graphs
  #***
  #**********************************************************************
  
  figure()
  par(mfrow=c(2,2))
  
  #--------------------------------------------------------#  
  plot(seqa(1,1,t),y_noise, type="l",
       main = "White noise",
       xlab = "t",
       ylab = expression(y[t]),
       bty = "l")
  

  plot(seqa(1,1,t),y_i0, type="l",
       main="Stationary: I(0)",
       xlab = "t",
       ylab = expression(y[t]),
       bty = "l")
  
  plot(seqa(1,1,t),y_i1, type="l", 
       main = "Nonstationary: I(1)",
       xlab = "t",
       ylab = expression(y[t]),
       bty = "l")
  
  plot(seqa(1,1,t),y_i2, type="l",
       main = "Nonstationary: I(2)",
       xlab = "t",
       ylab = expression(y[t]),
       bty = "l")
  
}
