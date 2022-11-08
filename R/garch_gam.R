#============================================================================
#
#  Efficiency of QMLE for a garch model with conditional gamma
#
#============================================================================

rm(list = ls(all=T))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, figure, inv
source("EMTSUtil.R")



#----------------------------------------------------------------------------
# Likelihood wrapper function
# This function calls the lnlt function and returns the average log-likelihood.
#----------------------------------------------------------------------------
lnl <- function( b,y ) {
  logl <- -mean( lnlt( b,y ) )
  return(logl)
}

#----------------------------------------------------------------------------
# Likelihood function for a GARCH-gam(1,1) model
#----------------------------------------------------------------------------
lnlt <- function( b,y ) {
  v  <- y                                                                     
  s2 <- recserar(cbind(1-b[1]-b[2] + b[1]*trimr(c(0.0, v^2),0,1)),cbind(sd(y)^2),cbind(b[2]))     
  c  <- abs(b[3])                                             
  z  <- v/sqrt(s2)                                                            
  
  loglt <- 0.5*log(c) - log(gamma(c)) + (c-1)*log(sqrt(c)*z + c)- 0.5*log(s2) - sqrt(c)*z - c
  return(loglt)
}
#----------------------------------------------------------------------------
# Likelihood wrapper function
# This function calls the lnlt function and returns the average log-likelihood.
#----------------------------------------------------------------------------
lnl1 <- function( b,y ) {
  logl <- -mean( lnlt1( b,y ) )
  return(logl)
}

#----------------------------------------------------------------------------
# Likelihood function for a GARCH-N(1,1) model
#----------------------------------------------------------------------------
lnlt1 <- function( b,y ) {
  b <- abs(b)
  t <- length( y )
  u <- y  
  h <- sd( y )^2*rep(1, t)
  
  for (i in 2:t) {
    h[i] <- 1-b[1]-b[2] + b[1]*u[i-1]^2 + b[2]*h[i-1]
  }
  loglt1 <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u/sqrt( h ))^2 
  return(loglt1)
}

#
# ------------------------ GARCH GAM model ----------------------------------
#
garch_gam <- function( ) {
  t <- 1000000
  z    <- rnorm(t,1)
  v    <- z                                                
  y    <- z
  s2   <- rep(1, t)
  
  # Parameters
  a1 <- 0.1
  b1 <- 0.8
  c  <- 50
  
  #  Simulate the dgp of a garch model with conditional gamm disturbance     
  for (i in 2:t) {
    z[i] <- (rgamma(1, shape=c,scale=1) - c)/sqrt(c)                             
    s2[i] <- 1-a1-b1 + a1*v[i-1]^2 + b1*s2[i-1]         
    v[i] <- z[i]*sqrt(s2[i])                            
    y[i] <- v[i]
  }
  
  
  # Analaytical expressions (based on Engle and Gonzalez-Rivera (JBES, 1991)     
  ds2a <- recserar( cbind(c(0,trimr(y^2,0,1)) - 1), cbind(0.0) , cbind(b1) )
  ds2b <- recserar( cbind(c(0,trimr(s2,0,1)) - 1), cbind(0.0) , cbind(b1) )
  ds2  <- cbind(ds2a,   ds2b)
  
  tmp <- ds2/s2
  cov_h    <- inv( 0.5*t((tmp)) %*% (tmp)/t )
  beta2    <- 3 + 6/c
  cov_j    <- inv( (beta2-1)*0.25*t((tmp)) %*% (tmp)/t )
  cov_qmle <- cov_h %*% inv(cov_j) %*% cov_h  
  
  tmp <- ds2 * (c - c*z^2)/( c + sqrt(c)*z )
  dl  <- -0.5*tmp / s2
  j0    <- t(dl) %*% dl/t
  cov_0 <- inv(j0)                                        
  
  # Compute relative efficiency of qmle (based on analytical derivatives)     
  
  cat('\nRelative efficiency of qmle using analytical derivatives')
  cat('\n    alpha1    beta1\n')
  print((diag(cov_0)/diag(cov_qmle)))
  cat('\n\n alpha1 = ',a1)
  cat('\n beta1  = ',b1)
  cat('\n shape  = ',c)
  cat('\n' )
  
  # GARCH(1,1) gamma distribution:  
  start <- c(a1, b1, c)
  estResults <- optim(start, lnl, y=y, method="BFGS", hessian=T)
  theta <- estResults$par
  hess <- estResults$hess
  
  vch = inv(hess)
  tt  = diag(vch)
  
  # GARCH(1,1) normal distribution: qmle covariance
  start  = c(a1, b1)
  estResults <- optim(start, lnl1, y=y, method="BFGS", hessian=T)
  theta0 <- estResults$par
  h0 <- estResults$hess
  
  ih   <- inv(h0)
  g0 <- numgrad(lnlt1,theta0,y)
  j0   = t(g0) %*% g0/t
  vc0  = ih %*% j0 %*% ih
  tt0  = diag(vc0)
  
  # Efficiency ratio based on numerical derivatives
  cat('\nNumerical results')
  cat('\n', (tt[1:2]/tt0) )
}



    
                     