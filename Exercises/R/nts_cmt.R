#==========================================================================
#
#   Properties of the Functional Central Limit Theorem for various moments
#
#==========================================================================
rm(list = ls(all=T))
graphics.off()
set.seed(145, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  recserar, seqa, trimr
source("EMTSUtil.R")
# load ks library - kde
library(ks)
#
#-------------------- Continuous Mapping Theorem ----------------------------
#

nts_cmt <- function()
{
  ndraws  <- 50000
  imoment <- 1             
  t       <- 500
  
  # Parameters
  delta <- 0.0
  phi   <- 1.0
  sig2  <- 1.0
  sig   <- sqrt(sig2)
  y0    <- 0.0
  
  m <- rep(0, ndraws)
  pb <- txtProgressBar(min=, max=ndraws, style=3)
  for (i in seq(ndraws)) {
    # Random walk with first observation discarded
    y  <- trimr( recserar(cbind(delta + sqrt(sig2)*rnorm(t+1)),cbind(y0),cbind(phi)) , 1 , 0)   
    
    if (imoment == 1) {
      # Sample mean with standardization based on t^(-0.5)
      m[i] <- (1/sig)*sum(y)*t^(-3/2)
    } else if (imoment == 2) {
      m[i] <- (1/sig)*sum( seqa(1,1,t)*y)*t^(-5/2)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat('\nSample size =              = ', t)
  cat('\n\nSample mean of m           = ', mean(m))
  cat('\nTheoretical mean of m      = ', 0.0)
  cat('\nSample variance of m       = ', sd(m)^2)
  cat('\nTheoretical variance of m  = ', 1/3)
  
  x <- seqa(-5,0.1,101)  
  fhat <- kde(m, h=0.07, eval.points=x)
  fnorm <- dnorm(x)
  plot(fhat, col="blue", xlab="", ylab="")
  lines(x,fnorm, col='red')
}

