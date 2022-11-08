#=========================================================================
#
#     Simulation of nonlinear and linear exponential models.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
source("EMTSUtil.R")


#
#--------------------------- Simulations  -----------------------------------
#

nls_simulate <- function() {
    
  # Simulate data 
  t   <- 50                 
  b0  <- 1.0
  b1  <- 0.05
  sig <- 0.5
  u   <- sig*rnorm(t)
  x   <- 1:t
  
  y1  <- b0*exp( b1*x + u)
  y2  <- b0*exp( b1*x ) + u
  
  #***************************************************************************
  #***
  #***     Generate graphs
  #***
  #***************************************************************************
  figure()
  par(mfrow=c(1,2), xaxs="i", yaxs="i" )
  
  # Panel (a)  
  plot(x, y1, type="l", col = "red",
       main = "(a) Levels",
       xlab = expression(x[t]),
       xlim = c(0,50),       
       ylab = expression(y[t]),
       ylim = c(0,20),
       bty="l")  
  lines (x, y2, lty=5, col = "blue")
  legend("topright", 
         c(expression(y[1]), expression(y[2])), 
         lty = c(1,5),
         col = c("red", "blue"), lwd = 1)
  
  # Panel (b)  
  plot(x, log(y1), type="l", col = "red",
       main = "(b) Logs",
       xlab = expression(x[t]),
       xlim = c(0,50),       
       ylab = expression(paste("ln y",""[t]*"")),
       ylim = c(-1,4),
       bty="l")  
  lines (x, log(y2), lty=5, col = "blue")
  legend("topright", 
         c(expression(y[1]), expression(y[2])), 
         lty = c(1,5),
         col = c("red", "blue"), lwd = 1)  
}



