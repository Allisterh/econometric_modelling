#=========================================================================
#
#      Estimate an artificial neural network using the neural algorithm  
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr
source("EMTSUtil.R")

#
#----------------------- Artificial Neural Networks -------------------------
#

nlm_neural <- function() 
{
  # Parameters
  nobs <- 2000                         
  phi    <- 0.2                    
  gam    <- 2.0
  delta0 <- -2.0
  delta1 <- 2.0
  sig2   <- 0.1
  
  
  # Simulate  data  					     
  y <- rep(0, nobs+100)
  w <- rep(0, nobs+100)
  u <- sqrt(sig2)*rnorm(nobs+100)
  
  for (t in 2:(nobs+100)) {
    w[t-1] <- 1/( 1 + exp( - ( delta0 + delta1*y[t-1] ) ) )
    y[t]   <- phi*y[t-1] + gam*w[t-1] + u[t]
  }
  
  y <- trimr(y,100,0)
  
  # Estimate the model using neural algorithm
  n    <- 10000                                     
  pmat <- array(0, c(n,5))
  ylag <- trimr(y,0,1)         
  y    <- trimr(y,1,0)         
  
  pb <- txtProgressBar(min=0, max=n, style=3)
  
  for (iter in seq(n)) {
    d0 <- -5 + 10*runif(1)
    d1 <- -5 + 10*runif(1)
    w  <- 1/( 1 + exp( - ( d0 + d1*ylag ) ) )   
    x  <- cbind(ylag,   w)
    b  <- lm(y ~ x - 1)$coef
    v  <- y - x %*% b
    s2 <- t(v) %*% v/length(y)
    setTxtProgressBar(pb, iter)
    pmat[iter,] <- c(b, d0, d1, s2) 
    
  }
  close(pb)
  
  # Sort the parameter estimates based on s2    
  pmat <- t(apply(pmat, 1, sort))
  
  cat('\nParameter estimates')
  cat('\nphi         = ',pmat[1,1])
  cat('\ngam         = ',pmat[1,2])
  cat('\ndelta0      = ',pmat[1,3])
  cat('\ndelta1      = ',pmat[1,4])
  cat('\nsigma^2     = ',pmat[1,5])
}




