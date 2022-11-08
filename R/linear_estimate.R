#=========================================================================
#
#    Program to estimate linear regression by maximum likelihood.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
# 
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Log-likelihood function 
#-------------------------------------------------------------------------

neglog <- function(theta,y,x1,x2) {
  lf <- -mean( lnlt(theta,y,x1,x2) )
  return(lf)
}

#-------------------------------------------------------------------------
# Log-likelihood function concentrated
#-------------------------------------------------------------------------
neglogc <- function(theta,y,x1,x2) {
  lf <- -mean( lnltc(theta,y,x1,x2) )
  return(lf)
}

#-------------------------------------------------------------------------
# Log-likelihood function at each observation
#-------------------------------------------------------------------------
lnlt <- function(theta,y,x1,x2) {
  m  <- theta[1] + theta[2]*x1 + theta[3]*x2                             
  s2 <- theta[4]                              	                   
  z  <- (y - m)/sqrt(s2)     
	lf <- -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z^2
  return(lf)
}

#-------------------------------------------------------------------------
# Concentrated Log-likelihood function at each observation
#-------------------------------------------------------------------------
lnltc <- function(theta,y,x1,x2) {
  m  <- theta[1] + theta[2]*x1 + theta[3]*x2   	
  u  <- y - m                                  	                  
  s2 <- t(u) %*% u/length(y)                                 		
	z  <- (y - m)/sqrt(s2)    
	lf <- -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z^2
  return(lf)
}


#
#------------------------ ML Estimation of Regression Model ----------------------------
#

linear_estimate <- function() {
  t <- 200
  
  # Parameter values
  beta0  <- 1.0                                                         
  beta1  <- 0.7
  beta2  <- 0.3
  sig    <- sqrt(4.0)
   
  # Simulate the data
  x1 <- rnorm(t)                                                                 
  x2 <- rnorm(t)
  u  <- sig*rnorm(t)

  y  <- beta0 + beta1*x1 + beta2*x2 + u 
  ones <- rep(1, t)    
  x  <- matrix(c(ones, x1, x2), ncol = 3)

  # Initial guess
  #theta_0 <- runif(4)
  theta_0 <- c(1.0136,
               1.0745,
               0.5363,
               3.7022)
   
  # Estimate the model
  results <- optim(theta_0, neglog, y=y, x1=x1, x2=x2, method="BFGS")
  theta <- results$par    

  cat('\nEstimated Parameters')
  cat('\n', theta )
    
  ht <- numhess(neglog,theta,y,x1,x2)
  cat('\n\nEstimated covariance matrix\n')
  print((1/t) * inv(ht))

  # Estimate the concentraned model using Newton-Raphson    
  theta_0 <- runif(3)   # Initial guess
  
  # Estimate the model
  results <- optim(theta_0, neglogc, y=y, x1=x1, x2=x2, method="BFGS")
  theta <- results$par

  cat('\nEstimated Parameters (Concentrated)')
  cat('\n', theta)
  
  htc <- numhess(neglogc,theta,y,x1,x2)
  cat('\n\nEstimated covariance matrix (concentrated)\n')
  print((1/t)*inv(htc))    

  # Compute OLS estimates
  theta <- lm(y~x[,-1])$coef
  cat('\nOLS Parameter estimates')
  cat('\n',theta)
  
  # Compute the covariance matrices    
  e       <- y - x %*% theta           
  sig2hat <- t(e) %*% e/t
  vcov    <- as.numeric(sig2hat) * inv( t(x) %*% x)
    
  cat('\n\nCovariance matrix (OLS)\n\n')
  print(vcov)

  cat('\nVariance estimate (OLS)          = ', sig2hat)
  cat('\nVariance estimate standard error = ', 2*sig2hat^2/t)  
}

