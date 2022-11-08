#=========================================================================
#
#   Maximum likelihood estimation of the transitional distribution of the 
#   CIR model of interest rates using Ait Sahalia's (1996) data.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
# 

#-------------------------------------------------------------------------
# Wrapper function to calculate inverse of a given matrix
#-------------------------------------------------------------------------
inv <- function (M) {
  return(solve(M))
}


#-------------------------------------------------------------------------
# Likelihood function for transitional distribution of CIR model (Bessel)
#-------------------------------------------------------------------------
cir1 <- function( p,data ) {
  alpha <- abs(p[1])
  mu    <- abs(p[2])
  sigma <- abs(p[3])
  
    
  rnow <- data[-1]
  rlag <- data[-length(data)]
  dt   <- 1/250
    
  c <- 2*alpha/(sigma^2 * (1-exp(-alpha*dt)))
  q <- 2*alpha*mu/sigma^2-1
  
  u <- c*exp(-alpha*dt)*rlag
  v <- c*rnow    
    
  lf <- log(c)-u-v+0.5*q*log(v/u)+log(besselI(2*sqrt(u*v), q, expon.scaled=TRUE))+2*sqrt(u*v)
  f  <- -sum( lf )
  cat('\np = [', p, ']','\tf(p) =', f)
  return(f)
}

    
#-------------------------------------------------------------------------
# Likelihood function for transitional distribution of CIR model
# (Chi-square)
#-------------------------------------------------------------------------
cir2 <- function ( p,data ) {
  alpha <- abs(p[1])
  mu    <- abs(p[2])
  sigma <- abs(p[3])
    
  rnow <- data[-1]
  rlag <- data[-length(data)]
  dt   <- 1/250
 
  c  <- 2*alpha/(sigma^2*(1-exp(-alpha*dt)))
  q  <- 2*alpha*mu/sigma^2-1
  u  <- c*exp(-alpha*dt)*rlag
  v  <- c*rnow 
  nc <- 2*u
  df <- 2*q+2 
  s  <- 2*v
    
  gpdf <- dchisq( s,df,nc )
  ppdf <- 2*c*gpdf

  f <- -sum(log( ppdf ) )
  cat('\np = [', p, ']','\tf(p) =', f)
  return(f)
}

# Load any compatibility requirements
source("EMTSUtil.R")

#
#--------------------------- Algorithms -----------------------------------
#

max_transitional <- function() {
   
  # Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
  #   1. year
  #   2. day
  #   3. date stamp
  #   4. interest rates
  load ("eurodollar.RData")
  
  dt <- 1/250
  rt <- eurodata[,4]
  
  # Starting values
  x  <- rt[-length(rt)]
  dx <- diff(rt)           
  dx <- dx/x^0.5  

  regressors <- cbind(c(dt/x^0.5), c(dt*x^0.5))
  reg <- lm(dx ~ regressors[,1] + regressors[,2] - 1)
  
  drift <- reg$coef                                   # Get the coefficients from the model  
  res <- reg$residuals                                # Get the residuals from the model

  alpha <- -drift[2]
  mu <- -drift[1]/drift[2]
  sigma <- sqrt(var(res)/dt)

  p0 <- c(abs(alpha), abs(mu), abs(sigma))
  
  # Estimation based on scaled Bessel function  
  results <- optim(p0, cir1, data = rt, method="BFGS", hessian=TRUE)  
  
  phat1 <- results$par
  hessian <- results$hessian
  
  hessinv <- inv(hessian)
  
  cat('\nParameter estimates\n')
  cat('\t', phat1)
  cat('\nStandard errors based on inverse Hessian\n')
  cat('\t', c(sqrt( hessinv[1,1] ), sqrt( hessinv[2,2] ), sqrt( hessinv[3,3]) ) )
  
  # Estimation based on non-central Chi-square function
  results <- optim(phat1, cir2, data = rt, method="BFGS", hessian=TRUE)  
  
  phat <- results$par
  hessian <- results$hessian
  
  hessinv <- inv(hessian) 
  
  cat('\nParameter estimates\n')
  cat('\t', phat)
  cat('\nStandard errors based on inverse Hessian\n')  
  cat('\t', c(sqrt( hessinv[1,1] ), sqrt( hessinv[2,2] ), sqrt( hessinv[3,3] )) )
  
  #********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************
  
  figure()
  par(mar=c(5, 5, 4, 1))
  rlag <- rt[-length(rt)]
  rnow <- rt[-1]
  
  tmp0 <- seq(0.03, 0.25, 0.01)
  tmp1 <- (phat[3]^2)*tmp0 
  plot(rlag, rnow^2, yaxt='n', pch=19,
       xlab = expression(r[t-1]),
       ylab = expression(r[t]^2))
  lines(tmp0, tmp1)
}   


