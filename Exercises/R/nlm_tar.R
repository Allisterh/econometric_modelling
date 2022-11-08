#============================================================================
#
#   Nonparametric density estimator of a threshold autoregressive model
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
nlm_tar <- function( ) {
  
  theta <- 0.8 
  k <- 1                  
  t <- 5000
  y <- seq(-3,3,0.1)
  
  
  # Simulate TAR(1) model
  yt <- theta*rep(0, t)
  
  for (i in 2:t) {
    yt[i]<- theta*abs(yt[i-1]) + sqrt(1 - theta^2)*rnorm(1)
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
  
  m <- fxy/ fx
  
  # Compute true conditional mean
  m_true <- theta*abs(y)
  
  
  # Compute linear conditional mean 
  m_linear <- y*theta^k
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  
  #--------------------------------------------------------#
  matplot(y,cbind(m,m_true,m_linear), type="l",
          main = "",
          ylab = expression(m(y[t-1])),
          xlab = expression(y[t-1]),
          bty = "l")
                             
}
          