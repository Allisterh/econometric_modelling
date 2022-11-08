#============================================================================
#
#   Program to perform tests of the stationary distribution of the 
#   interest rate based on the gamma distribution
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Likelihood function for stationary distribution of CIR model
#-------------------------------------------------------------------------
neglog <- function(b,y) {
  nu <- abs(b[1])
  om <- abs(b[2])
  f  <- nu*log( om ) - log(gamma( nu )) + (nu-1)*log(y) - om*y  
  lf <- -mean(f)
  return(lf)
}

test_interest <- function(  ) {
  # Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
  #   1. year
  #   2. day
  #   3. date stamp
  #   4. interest rates
  load("eurodollar.RData")
  
  r <- eurodata[,4]
  t <- length(r)
  
  # Estimate the model
  start <- c(1,1)
  estResults <- optim(start, neglog, y=r, hessian=T)
  bhat <- estResults$par
  hess <- estResults$hess
  
  vc    <- (1/t)*inv(hess)
  nu    <- bhat[1]
  omega <- bhat[2]
  
  cat('\nParameter estimates')
  cat('\nnu           = ',nu)
  cat('\nomega        = ',omega)
  
  cat('\n ')
  cat('\nHessian matrix\n')
  print( hess )
  
  # --------------------------------------------------------------
  # Estimate the mean and its standard error    
  mu <- nu/omega
  d  <- cbind((1/omega), (-nu/omega^2))
  vr <- d %*% vc %*% t(d)
  se <- sqrt(vr)
  
  cat('\n ')
  cat('\nMean              = ',mu)    
  cat('\nStd Error of Mean = ',se) 
  
  # Wald test of the mean 
  c <- mu
  q <- 0.1
  wd <- t*t(c - q) %*% inv(d %*% (inv(hess)) %*% t(d)) %*% (c - q)
  
  cat('\n ')
  cat('\nWald statistic = ',wd) 
  cat('\np-value        = ',1-pchisq(wd,1)) 
  
  # --------------------------------------------------------------
  # Estimate the variance and its standard error   
  s2 <- nu/omega^2
  d  <- cbind((1/omega^2),   (-2*nu/omega^3))
  vr <- d %*% vc %*% t(d)
  se <- sqrt(vr)
  
  cat('\n ')
  cat('\nVariance              = ',s2)   
  cat('\nStd Error of Variance = ',se) 
  
  # Wald test of the variance   
  c <- s2
  q <- 0.001
  wd <- t* t(c - q) %*% inv(d %*% (inv(hess)) %*% t(d)) %*% (c - q)
  
  cat('\n ')
  cat('\nWald statistic = ',wd)
  cat('\np-value        = ',1-pchisq(wd,1))
  
  # --------------------------------------------------------------
  # Estimate skewness and its standard error    
  k3   <- 2/sqrt(nu)
  d    <- cbind((-1/nu^(3/2)),  0)
  vr <- d %*% vc %*% t(d)
  se <- sqrt(vr)
  
  cat('\n ')
  cat('\nSkewness              = ',k3)    
  cat('\nStd Error of skewness = ',se)
  
  # Wald test of the skewness   
  c <- k3
  q <- 0.0
  wd <- t*t(c - q) %*% inv(d %*% (inv(hess)) %*% t(d)) %*% (c - q)
  
  cat('\n ')
  cat('\nWald statistic = ',wd)
  cat('\np-value        = ',1-pchisq(wd,1))
  
  # --------------------------------------------------------------
  # Estimate kurtosis and its standard error    
  k4   <- 3 + 6/nu
  d    <- cbind(-6/nu^2, 0)
  vr <- d %*% vc %*% t(d)
  se <- sqrt(vr)
  
  cat('\n ')
  cat('\nKurtosis              = ',k4)    
  cat('\nStd Error of kurtosis = ',se) 
  
  # Wald test of the kurtosis   
  c <- k4
  q <- 3
  wd <- t*t(c - q) %*% inv(d %*% (inv(hess)) %*% t(d)) %*% (c - q)
  
  cat('\n ')
  cat('\nWald statistic = ',wd)
  cat('\np-value        = ',1-pchisq(wd,1))
  
}

