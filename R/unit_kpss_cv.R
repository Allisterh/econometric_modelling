#============================================================================
#
#   KPSS asymptotic critical values
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#
# Load require functions - seqa
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  KPSS test
#----------------------------------------------------------------------------
kpsscv <- function(x,t,reps) {
  kpss <- rep(0, reps)
  for (j in seq(reps)) {
    y <- rnorm(t)
    z <- residuals(lm(y ~ x - 1))
    S <- cumsum(z)    
    kpss[j] <- (t(S)%*% S)/(t(z) %*% z %*% t)    
  } 
  cv <- quantile(kpss,c(0.90, 0.95, 0.99))  
  return (cv)
}


#
#------------- Testing the Null Hypothesis of Stationarity ------------------
#

unit_kpss_cv <- function() {  
  reps <- 100000
  t    <- 1000
  
  # Constant
  cv <- kpsscv(rep(1, t),t,reps)
  cat('\n Constant ' )
  cat('\n    0.10      0.05      0.01 ')
  cat('\n    -------------------------')
  cat('\n', cv, '\n')
  
  # Constant and trend
  cv <- kpsscv(cbind(rep(1, t), seqa(1,1,t)),t,reps) 
  cat('\n Constant and trend' )
  cat('\n    0.10      0.05      0.01 ')
  cat('\n    -------------------------')
  cat('\n', cv, '\n')
  
  tau=0.25 
  cv = kpsscv(cbind(rep(1, t),  (seqa(1,1,t)>floor(tau*t))),t,reps)
  cat('\n tau = 0.25 : Constant and trend' )
  cat('\n    0.10      0.05      0.01 ')
  cat('\n    -------------------------')
  cat('\n', cv, '\n')
  
  tau<-0.5 
  cv <- kpsscv(cbind(rep(1,t), (seqa(1,1,t)>floor(tau*t))),t,reps)
  cat('\n tau = 0.5 : Constant and trend' )
  cat('\n    0.10      0.05      0.01 ')
  cat('\n    -------------------------')
  cat('\n', cv, '\n')
                
  tau<-0.75 
  cv <- kpsscv(cbind(rep(1, t), (seqa(1,1,t)>floor(tau*t))),t,reps)
  cat('\n tau = 0.75 : Constant and trend' )
  cat('\n    0.10      0.05      0.01 ')
  cat('\n    -------------------------')
  cat('\n', cv, '\n')

   
}
    