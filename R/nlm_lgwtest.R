#=========================================================================
#
#      Lee, Granger, White neural network test J Ects (1993) 
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr
source("EMTSUtil.R")
#----------------------------------------------------------------------------
#  Estimate and test a neural network
#----------------------------------------------------------------------------
neural <- function(y,x,q) {
  
  activate <- array(0, c(length(y),q))
  
  for (i in seq(q)) {   
    gam    <- 4*runif(ncol(x))-2.0  #  Uniform random numbers [-2,2]  )
    lambda <- x %*% gam
    
    activate[,i] <- (1+exp(-lambda))^(-1)      
  }
    
  xx    <- cbind(x, activate)
  uhat  <- lm(y ~ x - 1)$residuals
  ehat  <- lm(uhat ~ xx - 1)$residuals
  ntest <- length(y)*(1 - sum(ehat^2)/sum(uhat^2))
  pred  <- xx %*% (lm(y ~ xx -1)$coef)         
  pval  <- 1 - pchisq(ntest,q)
  return(list(pred=pred, ntest=ntest, pval=pval))
}

#
#--------------------------- Neural newtwork test ---------------------------
#

nlm_lgwtest <- function( ) {
  
  # Parameters
  reps <- 2000       # Number of draws                 
  q    <- 2          # Number of activation functions    
  nobs <- 250 
  
  
  #  Null hypothesis: AR Model   
  ntest <- rep(0, reps)
  
  pb <- txtProgressBar(min=0, max=reps, style=3)
  for (r in seq(reps)) {
    y <- rep(0, nobs+100)
    v <- rnorm(nobs+100)
    for (t in 2:(nobs+100)) {
      y[t] <- 0.6*y[t-1] + v[t]
    }  
    y <- trimr(y,100,0)   
    x <- trimr(y,0,1)    
    y <- trimr(y,1,0)  
    
    results <- neural(y,cbind(rep(1, length(x)),x), q)
    ypred <- results$pred
    ntest[r] <- results$ntest
    pv <- results$pval
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  
  dist <- quantile(ntest,probs=c(.010, .025, .050, .100, .500, .900, .950, .975, .990))
  
  cat('\n')
  cat('\nEmpirical distribution (under the null of linearity)')
  
  cat('\n1.0 per cent      =  ', dist[1])
  cat('\n2.5 per cent      =  ', dist[2])
  cat('\n5.0 per cent      =  ', dist[3])
  cat('\n10.0 per cent     =  ', dist[4])
  cat('\n50.0 per cent     =  ', dist[5])
  cat('\n90.0 per cent     =  ', dist[6])
  cat('\n95.0 per cent     =  ', dist[7])
  cat('\n97.5 per cent     =  ', dist[8])
  cat('\n99.0 per cent     =  ', dist[9])
  
  cv <- dist[7]
  
  cat('\n')
  
  
  #  Alternative hypothesis: Bilinear Model   
  ntest <- rep(0, reps)
  
  pb <- txtProgressBar(min=0, max=reps, style=3)
  for (r in seq(reps)) {
    y <- rep(0, nobs+100)
    v <- rnorm(nobs+100)
    
    for (t in 3:(nobs+100)) {
      y[t] <- 0.7*y[t-1]*v[t-2] + v[t]
    }  
    y <- trimr(y,100,0)     
    x <- trimr(y,0,1)    
    y <- trimr(y,1,0)   
    
    results <- neural(y,cbind(rep(1, length(x)),x), q)
    ypred <- results$pred
    ntest[r] <- results$ntest
    pv <- results$pval
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  
  
  dist <- quantile(ntest,c(.010, .025, .050, .100, .500, .900, .950, .975, .990))
  
  cat('\n')
  cat('\nEmpirical distribution (under the alternative)')
  
  cat('\n1.0 per cent      =  ', dist[1])
  cat('\n2.5 per cent      =  ', dist[2])
  cat('\n5.0 per cent      =  ', dist[3])
  cat('\n10.0 per cent     =  ', dist[4])
  cat('\n50.0 per cent     =  ', dist[5])
  cat('\n90.0 per cent     =  ', dist[6])
  cat('\n95.0 per cent     =  ', dist[7])
  cat('\n97.5 per cent     =  ', dist[8])
  cat('\n99.0 per cent     =  ', dist[9])
  
  cv <- dist[7]
}
