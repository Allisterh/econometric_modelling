#============================================================================
#
#   Program to simulate the OLS estimator of the AR(1) model
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(5, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, recserar, figure
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# AR1 Procedure  
#----------------------------------------------------------------------------

ar1sim <- function( phi,t,nreps) {
   phihat <- rep(0,nreps)
   y <- rep(0, t)
   #  Loop over reps simulating and estimating AR(1) model
   for (j in seq(nreps)) {
     v         <- rnorm(t+1)
     y         <- trimr(recserar(cbind(v),cbind(0),cbind(phi)),1,0)   
     phihat[j] <- lm(y[2:t] ~ y[1:(t-1)] - 1 )$coef
   }
   if (phi < 1)
     z <- sqrt(t)*(phihat - phi)/sqrt(1 - phi^2) # Stationary case
    else
      z <- t*(phihat - phi)                       # Nonstationary case    
   return(z)  
}
 
#
#------------------ Sampling Distribution AR(1) LSE--------------------------
#

nts_distribution <- function () {
    
  nreps <- 50000           
  t     <- c(50,100,200,400,800,1600)
   
  # Compute statistics on the sampling distribution (stationary)
  cat('\nStatistics of the sampling distribution for phi = ', 0.8)
    
  xi <- seq(-5,5,0.05)
  fs  <- array(0, c(length(xi),length(t) ))  
  
  pb <- winProgressBar(min=0, max=length(t), title=paste("Sampling distribution (phi = 0.8, t = 0)"), "%0")
  for (j in seq(t) ) {    
    zs   <- ar1sim(0.8,t[j],nreps)
    bias <- mean(zs ) 
    stdv <- sd(zs) 
    skew <- mean( zs^3 )/stdv^3 
    kurt <- mean( zs^4 )/stdv^4   
    
    cat('\nT                         = ', t[j] )
    cat('\nBias                      = ', bias )
    cat('\nStandard deviation        = ', stdv )
    cat('\nSkewness                  = ', skew )
    cat('\nKurtosis                  = ', kurt )    
    cat('\n')    
    fs[,j] <- density(zs,from=xi[1], to=xi[length(xi)], n=length(xi))$y   
    setWinProgressBar(pb, j, paste(round(j/length(t) * 100), "%"), title=paste("Sampling distribution (phi = 0.8, t = ", t[j], ")"))
    
  }
  close(pb)
  
  tmps <- dnorm(xi)

  # Compute statistics on the sampling distribution (stationary)
  cat('\nStatistics of the sampling distribution for phi = ',1.0)
  
  xin <- seq(-15,5,0.01)
  fn  <- array(0, c(length(xin),length(t) ))
  
  pb <- winProgressBar(min=0, max=length(t), title=paste("Sampling distribution (phi = 1.0, t = 0)"), "%0")
  for (j in seq(t)) {
    zn <- ar1sim(1.0,t[j],nreps)     
    bias <- mean( zn ) 
    stdv <- sd( zn ) 
    skew <- mean( zn^3 )/stdv^3 
    kurt <- mean( zn^4 )/stdv^4
    cat('\nT                         = ', t[j] )
    cat('\nBias                      = ', bias )
    cat('\nStandard deviation        = ', stdv )
    cat('\nSkewness                  = ', skew )
    cat('\nKurtosis                  = ', kurt )    
    cat('\n')
    fn[,j] <- density(zn,from=xin[1], to=xin[length(xin)], n=length(xin))$y    
    setWinProgressBar(pb, j, paste(round(j/length(t) * 100), "%"), title=paste("Sampling distribution (phi = 1.0, t = ", t[j], ")"))
  }
  close(pb)

  tmpn <- dnorm(xin)
    
  #**********************************************************************
  #***
  #***     Generate graphs
  #***
  #**********************************************************************

  figure()
  par(mfrow=c(1,2))

  #--------------------------------------------------------#
  # Panel (a)
  plot(xi,tmps,type="l",
       main = "(a) phi = 0.8",
       xlab = expression(z),
       ylab = expression(f(z)),
       bty="l")
  lines(xi,fs[,1],lty=2)
  lines(xi,fs[,3],lty=3)

  #--------------------------------------------------------#
  # Panel (b)  
  plot(xin,tmpn,type="l",
       main = "(a) phi = 1.0",
       xlab = expression(z),
       ylab = expression(f(z)),
       bty = "l")
  lines(xin,fn[,1],lty=2)
  lines(xin,fn[,3],lty=3)  
}
