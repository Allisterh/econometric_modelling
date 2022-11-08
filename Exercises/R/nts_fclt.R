#============================================================================
#
#   Program to demonstrate the properties of the 
#   Functional Central Limit Theorem: 
#   Standardization of a random walk
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, trimr, recserar
source("EMTSUtil.R")


#
#------------------------- Central Limit Theorem-----------------------------
#
       
nts_fclt <- function()
{
  ndraws <- 50000
  t      <- 500
  
  alpha <- 0.0
  phi   <- 1.0
  sig2  <- 1.0
  y0    <- 0.0
  
  yts1 <- rep(0, ndraws)    # s = 1/4
  yts2 <- rep(0, ndraws)    # s = 1/2
  yts3 <- rep(0, ndraws)    # s = 1
  
  pb <- winProgressBar(min=0, max=ndraws, title="Simulating model", "%0")
  for (i in seq(ndraws)) {
    # Random walk with first observation discarded
    y  <- trimr( recserar(cbind(alpha + sqrt(sig2)*rnorm(t+1)),cbind(y0),cbind(phi)),1,0 )   
    
    yts1[i] <- y[round(t*0.25)]*t^(-0.5)/sqrt(sig2)                                                
    yts2[i] <- y[round(t*0.5)]*t^(-0.5)/sqrt(sig2)                                              
    yts3[i] <- y[length(y)]*t^(-0.5)/sqrt(sig2)  
    setWinProgressBar(pb, i, paste(round(i/ndraws * 100), "%"), title="Simulating model")
  }
  close(pb)
  
  cat('\nSample size = ',t)
  cat('\n')
  
  cat('\nSample mean of yts1 (T/4)               = ',  mean(yts1) )
  cat('\nTheoretical mean of yts1 (T/4)          = ',  0.0)
  cat('\nSample variance of yts1                 = ',  sd(yts1)^2)
  cat('\nTheoretical variance of yts1            = ',  1/4)  
  cat('\n')
  
  cat('\nSample mean of yts2 (T/2)               = ',  mean(yts2))
  cat('\nTheoretical mean of yts2 (T/2)          = ',  0.0)
  cat('\nSample variance of yts2                 = ',  sd(yts2)^2)
  cat('\nTheoretical variance of yts2            = ', 1/2)  
  cat('\n')
  
  cat('\nSample mean of yts3 (T)                 = ', mean(yts3))
  cat('\nTheoretical mean of yts3 (T)            = ',  0.0)
  cat('\nSample variance of yts3                 = ', sd(yts3)^2)
  cat('\nTheoretical variance of yts3            = ',  1/1)  
  cat('\n')
  
  
  xi <- seq(-5,5, 0.1)
  
  fhat1 <- density(yts1,from=xi[1], to=xi[length(xi)], n=length(xi))$y  
  fhat2 <-density(yts2,from=xi[1], to=xi[length(xi)], n=length(xi))$y  
  fhat3 <- density(yts2,from=xi[1], to=xi[length(xi)], n=length(xi))$y  
  fnorm <- dnorm(xi)
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i")
  
  #--------------------------------------------------------#
  matplot(xi,cbind(fnorm,fhat1,fhat2,fhat3), type="l",
          xlab = expression(paste("Y", ""[T]*"", "(s)")),
          ylab = expression(paste("f(Y", ""[T]*"", "(s))")),
          xlim = c(-4, 4),  
          bty="l")
}
