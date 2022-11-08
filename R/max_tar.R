#============================================================================
#
#   Program to estimate a threshold autoregressive model based on 
#   the stationary and transitional distributions
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12357, kind="Mersenne-Twister")

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr
source("EMTSUtil.R")  

#----------------------------------------------------------------------------
#  Log-likelihood of a tar model: transitional distribution
#----------------------------------------------------------------------------
lnltrans <- function(b,y) {
  m <- b*abs(trimr(y,0,1))
  s <- 1
  z <- ( trimr(y,1,0) - m ) / s
  f <- log( dnorm(z)/s )                                  
  
  lf <- -mean(f)
  return (lf)
}

#-------------------------------------------------------------------------
  #  Log-likelihood of a tar model: stationary distribution
#-------------------------------------------------------------------------
lnlstat<- function (b,y){
  m <- 0.0
  s <- 1/sqrt(1 - b^2)
  z <- ( y - m )/ s
  f <- 2*dnorm(z)*(1/s)*pnorm(b*y)    
  lf <- -mean( log(f) )  
  return(lf)
}

max_tar <- function ( ) {
  t     <- 250 
  theta <- 0.5  
  
  # Simulate data 
  y <- rep(0, t)
  
  for (i in 2:t) {
    
    y[i] <- theta*abs(y[i-1]) + rnorm(1)
  }
  
  # Estimate by least squares       
  xvar <- trimr(abs(y),0,1)
  yvar <- trimr(y,1,0) 
  bols <- lm(yvar ~ xvar - 1)$coef
  
  # Estimate by Ml applied to the stationary distribution       
  start <- bols
  estResults <- optim(start, lnlstat, y=y, method="BFGS", hessian=T)
  bhat <- estResults$par
  hess <- estResults$hessian
  
  
  
  
  cat('\nTrue value of theta', theta)
  
  cat('\nOLS Results', bols)
  
  cat('\nStationary distribution results', bhat)
  
  # Estimate by Ml applied to the transtional distribution       
  estResults <- optim(start, lnltrans, y=y, method="BFGS", hessian=T)
  bhat <- estResults$par
  hess <- estResults$hessian
  
  
  cat('\nTransitional distribution results', bhat)
  
}


