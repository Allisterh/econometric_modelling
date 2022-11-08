#============================================================================
#
#   Program to estimate a nonlinear regression model
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

set.seed(1, kind="Mersenne-Twister")

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv
source("EMTSUtil.R")

#---------------------------------------------------------------------------
# Log-likelihood function
#---------------------------------------------------------------------------
neglog <- function(b,y,x) {
  m   <- 1/(x - b)
  s2  <- 1
  lt <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m)^2)/s2
  lf <- -mean(lt)  
  return(lf)
}


# 
#--------------------------------- Regression Model ---------------------
#

nls_regression1 <- function( ) {
  t <- 100
  
  # Generate the data                   
  theta <- 2
  x <- 2 + runif(t)
  u <- rnorm(t)
  y <- 1/( x - theta) + u
  
  # Method of Scoring   
  start <- 1.5              
  g <- mean( (y - 1/(x - start)) * (1/(x - start)^2) )
  i <- mean( (1/(x - start)^4) )
  
  theta1 <- start + inv(i) %*% g
  
  cat('\nMethod of Scoring')
  cat('\nStarting value of theta = ',start)
  cat('\nUpdated value of theta  = ',theta1)
  
  # Estimate the model using BGS and compute Hessian se  
  estResults <- optim(start, neglog, y=y, x=x, method="BFGS", hessian=T)
  bhat <- estResults$par
  hess <- estResults$val
  
  vc <- (1/t)*inv(hess)
  cat('\n ')
  cat('\nBFGS estimate of theta  = ',bhat)
  cat('\nStd. error  of theta    = ',sqrt(vc))
  
}



