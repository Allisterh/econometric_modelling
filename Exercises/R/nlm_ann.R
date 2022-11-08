#============================================================================
#
#      Estimate an artificial neural network by maximum likelihood  
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Log-likelihood function
#----------------------------------------------------------------------------
neglog <- function(b,y)
{
  nobs <- length(y)
  v    <- rep(0, nobs)
  w    <- rep(0, nobs)
  
  for (t in 2:nobs) {
    w[t-1] <- 1/( 1 + exp( - ( b[3] + b[4]*y[t-1] ) ) )
    v[t]   <- y[t] - b[1]*y[t-1] - b[2]*w[t-1]
  }
  v  <- trimr(v,1,0)   
  s2 <- mean(v^2) 
  z  <- v/sqrt(s2)
  lf <- -mean( -0.5*log(2*pi)-0.5*log(s2)-0.5*z^2 )
  return(lf)  
}


#
#----------------------- Artificial Neural Networks -------------------------
#
nlm_ann <- function()
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
  
  for (t in 2:nobs+100) {
    w[t-1] <- 1/( 1 + exp( - ( delta0 + delta1*y[t-1] ) ) )
    y[t]   <- phi*y[t-1] + gam*w[t-1] + u[t]
  }
  
  y <- trimr(y,100,0)
  
  # Estimate the model using BFGS
  theta_0 <- c(phi,gam, delta0, delta1)
  estResults <- optim(theta_0, neglog, y=y, method="BFGS")
  theta <- estResults$par
  lnl <- estResults$val
  
  cat('\n      True  Estimated\n')
  print(unname(cbind(theta_0, theta)))
  cat('\nlnl = ', -lnl)
}





