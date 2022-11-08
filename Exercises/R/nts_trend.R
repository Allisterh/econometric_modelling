#=========================================================================
#
#   Investigate the effects of misspecifying the trend 
#
#=========================================================================
rm(list = ls(all=T))
graphics.off()
set.seed(145, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  recserar, trimr, figure, inv
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Regression estimates of trend models
#-------------------------------------------------------------------------
ar <- function(y,ntrend) {
  t <- length(y)
  # ntrend <- 0 (constant and lag)
  if (ntrend == 0.0) {
    x <- cbind(rep(1, t-1), trimr(y,0,1))
    y <- trimr(y,1,0)
  }else if (ntrend == 1.0) {       # ntrend = 1 (constant, time trend and lag)    
    x <- cbind(rep(1, t-1), seqa(0,1,t-1)/100,  trimr(y,0,1))
    y <- trimr(y,1,0)   
  }else {     # ntrend = 2 (constant and time trend)
    x <- cbind(rep(1, t),  seqa(0,1,t)/1)                        
    y <- trimr(y,0,0)
  }  
  
  b <- lm(y ~ x - 1)$coef
  e <- y - x %*% b                                 
  s2 <- t(e) %*% e/length(e) 
  s2 <- as.numeric(s2)
  se <- sqrt( diag( s2 * inv(t(x) %*% x) ) )   
  return (list (b=b, se=se))
}

#
#----------------------------- Trend Estimation -----------------------------
#

nts_trend <- function()
{
  # DGP: Simulate I(1) variable
  t     <- 100
  delta <- 0.5
  sig2  <- 1.0
  y0    <- 0.0
  y     <- recserar(cbind(delta + sqrt(sig2)*rnorm(t)),cbind(y0),cbind(1.0))          
  
  # Estimate the deterministic trend model    
  model <- ar(y,2)
  b <- model$b
  se <- model$se
  
  cat('\n')
  print(cbind(b=b, "Std Err"=se, "t-stat"= b/se))
  
  # Sampling distribution of the slope parameter estimate 
  # of the deterministic trend model   
  ndraws <- 50000
  t      <- 5000
  
  b1 <- rep(0, ndraws)
  
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (i in seq(ndraws)) {
    # DGP is I(1)
    y      <- recserar(cbind(delta + sqrt(sig2)*rnorm(t)),cbind(y0),cbind(1.0))                                           
    b <- ar(y,2)$b                               
    b1[i] <- b[2]
    setTxtProgressBar(pb, i)
  }
  
  # Compute the standardized statistic   
  stat <- sqrt(t)*(b1 - delta)                        
  
  cat('\n')
  cat('\nMean of scaled statistic     = ', mean(stat)) 
  cat('\nMean (theoretical)           = ', 0.0)
  cat('\nVariance of scaled statistic = ', sd(stat)^2) 
  cat('\nVariance (theoretical)       = ', sig2*6/5)
  
  
  hist(stat,31, col="darkblue", 
       main="", xlab="")  
}

