#=========================================================================
#
#   Program to demostrate Anderson's CLT
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr
source("EMTSUtil.R")

#
#------------------------ Finite Sample Distribution--------------------------
#

stsm_finite <- function() {
  # Generate the data
  t      <- 50
  sigv   <- 1.0
  phi1   <- 0.9
  ndraws <- 10000
  
  ybar   <- rep(0, ndraws)
  theta_mle <- rep(0, ndraws)
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (i in seq(ndraws)) {
    vt <- sqrt(sigv)*rnorm(t+101)
    yt <- rep(length(vt))
    
    # Simulate the AR(1) model
    for (j in 2:length(vt)) {
       yt[j] <- phi1*yt[j-1] + vt[j]      
    }
    # Get rid of first 100 observations
    yt <- trimr( yt,100,0)  
    
     # Conditional mle
     theta_mle[i] <- lm(trimr( yt,1,0 ) ~ trimr( yt,0,1 ) - 1)$coef    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # Compute statistics of sampling distribution
  mse_mle  <- mean( (theta_mle - phi1)^2 )     
  rmse_mle <- sqrt( mse_mle )
        
  cat('\nPopulation parameter    =',  phi1 )
  cat('\nMean (cond. mle)        =',  mean(theta_mle)  )
  cat('\nBias (cond. mle)        =',  mean( theta_mle )- phi1 )
  cat('\nBias (Shenton-Johnson)  =',  -2*phi1/t )
  cat('\nRMSE (cond. mle)        =',  rmse_mle   )  
}
