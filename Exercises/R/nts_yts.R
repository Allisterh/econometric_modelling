#============================================================================
#
#   Program to construct the step function of Y[ts]
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure
source("EMTSUtil.R")


#
#------------------------- Central Limit Theorem-----------------------------
#

nts_yts <- function() {
  # Generate simulated values for t = 5
  t  <- 5                              
  y0 <- 2
  vt5 <- rnorm(t)
  yt5 <- rep(0, t)
  
  for (j in seq(t)) {
    yt5[j] <- y0 + sum( vt5[1:j] )  
  }
  
  
  yts5 <- rep(0, 100)
  
  k <- 1
  n <- seq(0.01, 1, 0.01)
  for (s in n) {
    p <- floor(s*t)
    if (p == 0)
      yts5[k] <- y0
    else
      yts5[k] <- y0 + sum( vt5[1:p] )
    k <-k+1
  }
  
  # Generate simulated values for t=40
  t  <- 40                              
  y0 <- 2
  vt <- rnorm(t)
  yt <- rep(0, t)
  
  for (j in seq(t)) {
    yt[j] <- y0 + sum( vt[1:j] )
  }
  
  yts <- rep(0, 100)
  
  k <- 1
  for (s in n) {
    p <- floor(s*t)
    if (p == 0)
      yts[k] <- y0
    else
      yts[k] <- y0 + sum( vt[1:p] )
    k <-k+1  
  }
      
     
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(mfrow=c(2,2), xaxs="i", yaxs="i")
  
  #--------------------------------------------------------#
  # Panel (a)
  plot(seq(5),yt5, pch=19,
       main = '(a) Discrete Representation T=5',
       xlab = 't',
       ylab = expression(y[t]),
       xlim = c(0,5),
       ylim = c(1,3),
       bty = "l")
       
  #--------------------------------------------------------#
  # Panel (b)
  plot(seq(40),yt, pch=19,
       main = '(b) Discrete Representation T=40',
       ylab = expression(y[t]),
       xlab = 't',
       xlim = c(0,40),
       bty = "l")
    
  #--------------------------------------------------------#
  # Panel (c)
  plot(n,yts5,type="l",
       main = '(c) Continuous Representation T=5',
       ylab = expression(y["Ts"]),
       xlab = 't',
       xlim = c(0,1),
       ylim = c(1,3),
       bty = "l")
  
  
  #--------------------------------------------------------#
  # Panel (d)
  plot(n,yts,type="l",
       main = '(d) Continuous Representation T=40',
       ylab = expression(y["Ts"]),
       xlab = 't',
       xlim = c(0,1),
       bty = "l")

  
}


