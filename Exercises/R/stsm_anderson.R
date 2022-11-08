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
#------------------------ Anderson's Central Limit Theorem ------------------
#

stsm_anderson <- function( ) {
      
  # Generate the data
  t      <- 200
  lags   <- 10
  sig2   <- 0.1
  mu     <- 0.1
  phi1   <- 0.9
  ndraws <- 10000
  ybar   <- rep(0, ndraws)
    
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (i in seq(ndraws)) {          
    vt <- sqrt(sig2)*rnorm( t+101)
    yt <- rep(0, length( vt ))
    
    # Simulate the AR(1) model
    for (j in 2:length( vt )) {
      yt[j] <- mu + phi1*yt[j-1] + vt[j]      
    }
    # Get rid of first 100 observations
    yt <- trimr( yt,100,0)        
    ybar[i] <- mean(yt)
    setTxtProgressBar(pb, i)                  
  }
  close(pb)
  tvar <- (sig2/(1-phi1^2))*(2/(1-phi1)-1)/t

  cat('\nSimulation mean      =', mean(ybar) )
  cat('\nTheoretical mean     =', mu/(1-phi1) )
  cat('\nSimulation variance  =', mean( (ybar - mu/(1 - phi1))^2 ))
  cat('\nTheoretical variance =', tvar)
  
}

