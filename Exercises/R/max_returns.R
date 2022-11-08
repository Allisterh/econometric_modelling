#============================================================================
#
#   Program to Program to estimate the distribution of asset returns 
#   using a standardised Student t distribution
#
#   Asset price data from 6 August 2010 to 2 January 2001 (note that the
                                                           #   data are in reverse order ie from recent to past)
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - figure, inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Log-likelihood of a Student t distribution
#----------------------------------------------------------------------------
lnl <- function(b,z) {
  m   <- b[1]
  s   <- abs(b[2])              
  v   <- abs(b[3])              
  s2  <- s^2                   
  tmp <- log( gamma((v+1)/2) ) - 0.5*log(pi) - log(gamma(v/2)) - 0.5*log(v-2) - 0.5*log(s2) - ((v+1)/2)*log( 1 + ((z - m)^2)/(s2*(v-2)) )
  
  lf  <- -mean(tmp)
  return (lf)
}

#----------------------------------------------------------------------------
#  Log-likelihood of a normal distribution
#----------------------------------------------------------------------------
lnl1 <- function (b,z) {
  m   <- b[1]                  
  s   <- abs(b[2])                                           
  tmp <- log(dnorm(z,m,s))
  
  lf  <- -mean(tmp)
  return (lf)
}

max_returns <- function ( ) {
  # Load data
  load("diversify.RData")
  
  t <- 2413
  
  # Select appropriate sample 
  pt_apple     <- as.matrix(pt_apple)[1:t]
  pt_ford      <- as.matrix(pt_ford)[1:t]
  
  # Compute percentage returns  
  r_apple <- 100*diff(log(pt_apple)) 
  r_ford  <- 100*diff(log(pt_ford))
  
  # Select a return 
  y <- r_apple
  
  # Compute desrciptive statistics   
  m <- mean(y)
  s <- sd(y)
  z <- (y - m)/s
  
  skew <- mean( z^3)
  kurt <- mean( 4^3)
  
  cat('\nMean     =   ',m)
  cat('\nStd dev. =   ',s)
  cat('\nSkewness =   ',skew)
  cat('\nKurtosis =   ',kurt) 
  cat('\nMaximum  =   ',max(z))
  cat('\nMinimum  =   ',max(z))
  
  figure()  
  hist(z, breaks=21, col="blue", prob=T,
       main = '',
       xlab = '',
       ylab = '')
  x <- seq(min(z), max(z), by=0.01)
  curve(dnorm(x, mean=mean(z),sd= sd(z)), add=T)
  
  # Estimate the model by MLE with Student t distribution 
  start <- c(m, s, 5)
  estResults <- optim(start, lnl, z=z, method="BFGS", hessian=T)
  bhat <- estResults$par
  lf <- estResults$val
  hess <- estResults$hess
  
  lf <- -lf
  vc <- 1/(t-1)*inv(hess)
  
  cat('\nUnrestricted model results')
  cat('\nLog-likelihood function = ',lf, '\n')
  print(cbind("Param"=bhat, "Std Err"= sqrt(diag(vc))))
  
  # Estimate the model by MLE with normal distribution 
  estResults <- optim(start[1:2], lnl1, z=z, method="BFGS", hessian=T)
  bhat <- estResults$par
  lf <- estResults$val
  hess <- estResults$hess
  
  lf <- -lf
  vc <- 1/(t-1)*inv(hess)
  
  
  cat('\nRestricted model results')
  cat('\nLog-likelihood function = ',lf, '\n')
  print(cbind("Param"=bhat, "Std Err"= sqrt(diag(vc))))
}



