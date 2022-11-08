# ========================================================================
#
#    Simulate Brownian motion and compare continuous and discrete data
#
# ========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure
source("EMTSUtil.R")

#
#------------------------- Brownian Motion ----------------------------------
#

sim_brown <- function() {
  dt <- 1/60                      # Minute data
  n  <- 240                       # 10 days = 240hrs 
  t <- n/dt                       # Sample size
  
  # Parameters
  mu   <- 0.0                     # Mean
  sig2 <- 1.0                     # Variance
  
  # Wiener process
  dw <- sqrt(dt)*rnorm(t)
  
  # Generate continous data (minutes)
  y <- rep(0, t)
  y[1] <- 0.0
  
  for (i in seq(t-1)) {
    y[i+1]<-y[i]+mu*dt + sqrt(sig2)*dw[i]  
  }
      
  y10 <- y[seq(10,t,10)]               # Generate 10 minute data by choosing every 10th observation
  yhr <- y[seq(10,t,60)]               # Generate hourly data by choosing every 60th observation 
  ydy <- y[seq(1440,t,1440)]           # Generate daily data by choosing every 1440th observation        
  
  
  #**************************************************************************
  #**
  #**     Generate graph
  #**
  #**************************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(2,2))
  #--------------------------------------------------------#
  # Panel (a)
  ind <- seq(t)/1440
  
  plot(ind,y,type="l",
       main = '(a) Minute Data',
       xlab = 't days',
       ylab = expression(y[t]),
       xlim = c(0, 10),
       ylim = c(-10, 30),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (b)
  
  ind <- seq(10, t, by=10)/1440
  plot(ind,y10,type="l",
       main = '(b) Ten Minute Data',
       ylab = expression(y[t]),
       xlab = 't days',
       xlim = c(0, 10),
       ylim = c(-10, 30),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (c)
  ind <-  seq(10, t, by=60)/1440
  plot(ind,yhr,type="l",
       main = '(c) Hourly Data',     
       ylab = expression(y[t]),
       xlab = 't days',
       xlim = c(0, 10),
       ylim = c(-10, 30),
       bty="l")
  #--------------------------------------------------------#
  # Panel (d)
  ind <-  seq(10, t, by=1440)/1440
  plot(ind,ydy,type="l",
       main = '(c) Daily Data',
       ylab = expression(y[t]),
       xlab = 't days',
       xlim = c(0, 10),
       ylim = c(-10, 30),
       bty="l")
}
