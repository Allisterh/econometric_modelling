#=========================================================================
#
#   Sample moments of stochastic and deterministic trend models
#
#=========================================================================
rm(list = ls(all=T))
graphics.off()
set.seed(145, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, recserar, seqa, trimr
source("EMTSUtil.R")

#
#----------------- Properties of Moments for Alternative Models -------------
#

nts_moment <- function ()
{
  #  Choose model
  # itrend = 0 for stochastic trend model 
  # itrend = 1 for deterministic trend model
  itrend <- 0
  
  # Parameters
  t      <- 800
  ndraws <- 50000
  sig2   <- 1.0
  y0     <- 0.0
  
  # Parameters of the stochastic trend model
  delta <- 0.0                
  phi   <- 1.0
  
  #  Parameters of the deterministic trend model
  beta1 <- 0.2
  beta0 <- 0.1                
  
  # Initialise arrays
  m1 <- rep(0, ndraws)
  m2 <- rep(0, ndraws)        
  m3 <- rep(0, ndraws) 
  m4 <- rep(0, ndraws)       
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (i in seq(ndraws))
  {
    if (itrend == 0) {
      # AR(1) model with first observation discarded
      y    <- trimr( recserar(cbind(delta + sqrt(sig2)*rnorm(t+1)),cbind(y0),cbind(phi)) , 1 , 0)           
      ylag <- trimr(y,0,1)
      
      if (phi < 1.0) {
        # Standardization based on t for stationary case
        m1[i] <- sum(ylag^1)/t^1                                                            
        m2[i] <- sum(ylag^2)/t^1                    
        m3[i] <- sum(ylag^3)/t^1                
        m4[i] <- sum(ylag^4)/t^1        
      }else {
        # Standardization for nonstationary case
        m1[i] <- sum(ylag^1)/t^1.5                                                              
        m2[i] <- sum(ylag^2)/t^2                       
        m3[i] <- sum(ylag^3)/t^2.5                  
        m4[i] <- sum(ylag^4)/t^3
      } 
    }else if (itrend == 1){
      trend <- seqa(1,1,t)
      y     <- beta0 + beta1*trend + sqrt(sig2)*rnorm(t)
      
      # Standardization based on t for stationary case
      m1[i] <- sum(trend^1)/t^2                                                               
      m2[i] <- sum(trend^2)/t^3
      m3[i] <- sum(trend^3)/t^4
      m4[i] <- sum(trend^4)/t^5 
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  
  
  cat('\nSample size    = ',t ,'\n')
  cat('\nMean of m1     = ',mean(m1))  
  cat('\nVariance of m1 = ',var(m1) ,'\n')  
  cat('\nMean of m2     = ',mean(m2))
  cat('\nVariance of m2 = ',var(m2) ,'\n')  
  cat('\nMean of m3     = ',mean(m3))
  cat('\nVariance of m3 = ',var(m3) ,'\n')
  cat('\nMean of m4     = ',mean(m4))
  cat('\nVariance of m4 = ',var(m4) ,'\n')
  
}
