#============================================================================
#
#   Program to demonstrate properties of the LM test
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, recserar
source("EMTSUtil.R")

#
#-------------- Properties of LM Tests of ARMA Models -----------------------
#

stsm_lm <- function() {
  # Sample sizes
  t <- 200
  
  # Simulate the data (discard first 100 simulated observations)
  mu   <- 0.0            
  phi1 <- 0.8           
  si1  <- 0.3
  vt   <- rnorm(t+101)
  yt   <- recserar( cbind(mu + trimr(vt,1,0) + si1*trimr(vt,0,1)) , cbind(0.0) , cbind(phi1))  
  yt   <- trimr(yt,100,0)   
  
  # Demonstrate the equivalence result  
  # LM test of AR(1)
  # First regression
  y <- yt   
  x <- rep(1, length(y))
  v <- lm(y ~ x - 1)$residuals
  
  # Second regression
  y   <- trimr(v,1,0)    
  x   <- cbind(rep(1, length(y)), trimr(v,0,1))
  e   <- lm(y ~ x - 1)$residuals
  y   <- y - mean(y)
  rsq <- 1 - (t(e) %*% e)/(t(y) %*% y)     
  lm  <- (t-1)*rsq
  
  cat('\n ')
  cat('\nLM statistic (phi1 = 0.0)      = ', lm) 
  cat('\np-value                        = ', 1-pchisq(lm,1))
  cat('\n ')
  
  # LM test of MA(1)
  # First regression
  y <- yt 
  x <- rep(1, length(y))
  e <- lm(y ~ x - 1)$residuals
  
  # Second regression
  y   <- trimr(e,1,0)    
  x   <- cbind(rep(1, length(y)),   trimr(e,0,1))
  e   <- lm(y ~ x - 1)$residuals
  y   <- y - mean(y)
  rsq <- 1 - (t(e) %*% e)/(t(y) %*% y)     
  lm  <- (t-1)*rsq
  
  cat('\n')
  cat('\nLM statistic (si1 = 0.0)       = ', lm) 
  cat('\np-value                        = ', 1-pchisq(lm,1))
  cat('\n')
  
  
  # Demonstrate singularity result
  # LM test (joint AR and MA)   
  # First regression
  y <- yt   
  x <- rep(1, length(y))
  v <- lm(y ~ x - 1)$residuals
  
  # Second regression
  y <- trimr(v,1,0)                   
  x <- cbind(rep(1, length(y)), trimr(yt,0,1), trimr(v,0,1))
  
  cat('\nCorrelation of explanatory variables at second stage\n')
  print( cor( x[,c(2,3)]) )
}

