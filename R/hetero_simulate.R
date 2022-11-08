#==============================================================================
#
#     Simulating a model of heteroskedasticity 
#
#==============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#------------------------- Helper Functions------------------------------------
#
#load required functions - figure
source("EMTSUtil.R")

#
#------------------------- Model Simulation -----------------------------------
#

hetero_simulate <- function() {
  t <- 5000
  alpha <- 1
  beta  <- 2
  delta <- 1
  gam   <- 0.5
  x <- rnorm(t)                         #     xt is a time trend                      
  w <- seq(0.0, 0.1*(t-1), 0.1)         #     wt is a time trend   
  u <- sqrt(delta + gam*w)*rnorm(t) 	  #     ut is N(0,sig2)                         
  y <- alpha + beta*x + u             	#     yt                                      
             
  
  
  #**********************************************************************
  #***
  #***     Generate Graphs
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(2,2), mar=c(5, 5, 4, 1))
  
  #--------------------------------------------------------#
  # Panel (a) 
  plot(x,y, type="l",
       main = "(a)",
       xlab = expression(x[t]),
       ylab = expression(y[t]),
       xlim = c(-4, 4),
       bty = "l")

  
  #--------------------------------------------------------#
  # Panel (b)
  plot(w,y, type="l",
       main = "(b)",
       xlab = expression(w[t]),
       ylab = expression(y[t]),
       xlim = c(0, 500),
       bty = "l")
  
  #--------------------------------------------------------#
  # Panel (c)
  plot(w,y^2, type="l",
       main = "(c)",
       xlab = expression(w[t]),
       ylab = expression(y[t]^2),
       xlim = c(0, 500),      
       bty = "l")
 
  
  #--------------------------------------------------------#
  # Panel (d)
  x <- cbind(array(1, c(t,1)), x)
  u <- lm(y ~ x- 1)$res
  
  plot(w,u^2, type="l",
       main = "(d)",
       xlab = expression(w[t]),
       ylab = expression(hat(u)[t]^2),
       xlim = c(0, 500),      
       bty = "l")
}
