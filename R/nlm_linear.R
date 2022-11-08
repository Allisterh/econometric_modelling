#============================================================================
#
#   Nonparametric estimator of a linear autoregressive model
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(3, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - figure
source("EMTSUtil.R")

#
#--------------------------- Nonparameteric Autoregression ------------------
#

nlm_linear <- function( ) {
  
  theta <- 0.8 
  k <- 1                  
  t <- 5000
  y <- seq(-3, 3, 0.1)
  
  # Simulate AR(1) model
  mu  <- 0.0
  phi <- 0.8
  yt <- mu + phi*rep(0, t)
  
  for (i in 2:t) {
    yt[i] <- mu + phi*yt[i-1] + rnorm(1)
  }
  
  # Kernel regression of y_t on y_t-k 
  
  fx  <- rep(0, length(y))
  fxy <- rep(0, length(y))
  
  h <- 1.06*sd(yt)*t^(-1/5)
  
  
  for (i in seq(y)) {
    z      <- ((y[i] - trimr(yt,0,k))/h)    
    fx[i]  <- mean( dnorm(z)/h )
    fxy[i] <- mean( dnorm(z)*trimr(yt,k,0)/h )    
  }
  m <- fxy / fx
  
  # Compute true conditional mean
  m_true <- theta*y
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  figure()
  
  #--------------------------------------------------------#
  matplot(y,cbind(m,m_true), type="l",
          main = "",
          ylab = expression(m(y[t-1])),
          xlab = expression(y[t-1]),
          bty = "l")
  
}