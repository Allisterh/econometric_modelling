#============================================================================
#
#   Program to generate the finite sampling distribution of the 
#   maximum likelihood estimator of the Poisson autoregressive model 
#   using the binomial thinning operator.
#
#   Results also given for the conditional least squares estimator.
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions ------------------------------
# 
# Load required functions -  trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Log-likelihood function of Poisson autoregression
#----------------------------------------------------------------------------
neglog <- function(b,y) {
 
  # Restrict domain of parametes
  rho <- pnorm(b[1])   
  lam <- exp(b[2])           
  
  t <- length(y)
  f <- rep(0, t)
  
  for (i in 2:t) {
    sum <- 0.0      
    tmp <- c(y[i],y[i-1])
    
    for (k in 0:min(tmp)) {
      sum1 <- factorial(y[i-1])/(factorial(k)*factorial(y[i-1]-k))
      sum2 <- rho^k*(1-rho)^(y[i-1]-k)
      sum3 <- lam^(y[i]-k)*exp(-lam)/factorial(y[i]-k)
      sum  <- sum + sum1*sum2*sum3      
    }
    f[i] <- sum
  } 
 
  # Exclude first observation from the likelihood
  f  <- trimr(f,1,0)               
  lf <- -mean( log(f) )
}


#
#--------------- Finite Sample Properties of Binomial thinning model  -------
#
discrete_poissonauto <- function() {
  # Parameters
  t      <- 100                                                                      
  ndraws <- 500                          
  rho    <- 0.3                              
  lam    <- 3.5
  theta0 <- c(rho, lam)
  
  # Loop to generate sampling distributions   
  theta_mle <- array(0, c(ndraws,2))
  theta_cls <- array(0, c(ndraws,2))
  
  # Starting values (transformed to satisfy domain restrictions)
  theta_0 <- c(qnorm(rho,0,1), log(lam))
  
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (j in seq(ndraws)) {
    
    # Generate data
    u <- rpois(t, lam)
    
    # Initialize y by choosing the median of draws from the 
    # unconditional distribution to ensure that the starting value is positive
    y <- rep(1,t)*median(rpois(t, lam/(1-rho))) 
    
    i <- 2
    while (i <= t) {
      e <- runif(y[i-1])   
      b_thin <- sum( e < rho )                          
      y[i] <- b_thin + u[i]   
      
      # Ensure that y(i) is a positive 
      if (y[i] > 0)
        i <- i + 1      
    }                       
    # Estimate by MLE
    estResults <- optim(theta_0, neglog, y=y, method="BFGS")
    theta <- estResults$par
    
    theta_mle[j,] <- c(pnorm(theta[1]), exp(theta[2]))                     
    
    # Estimate by CLS     
    xx <- cbind(trimr(y,0,1),  rep(1, t-1))
    yy <- trimr(y,1,0)
    
    
    b_cls <- lm(yy ~ xx - 1)$coef
    theta_cls[j,] <- b_cls
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  # Compute statistics of sampling distributions  
  
  mse_mle    <- colMeans( t( apply(theta_mle, 1, '-', theta0))^2)
  rmse_mle   <- sqrt( mse_mle )
  mse_cls    <- colMeans( t( apply(theta_cls, 1, '-', theta0))^2)
  rmse_cls   <- sqrt( mse_cls )
  cat('\n')
  cat('\nMSE  (MLE)                 = ',mse_mle)
  cat('\nRMSE (MLE)                 = ',rmse_mle)
  cat('\n')
  cat('\nMSE  (CLS)                 = ',mse_cls)
  cat('\nRMSE (CLS)                 = ',rmse_cls)
  cat('\n')
  cat('\nEfficiency (mle/cls)       = ',(mse_mle/mse_cls))  
}
