#============================================================================
#
#   Simulating GARCH - Student t Model
#  
#============================================================================
rm(list = ls(all=T))
graphics.off()
set.seed(1235, kind="Mersenne-Twister")

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
# Likelihood function for a GARCH-t(1,1) model
#----------------------------------------------------------------------------
lnlt <- function( b,y ){
  v  <- y                                                                     
  s2 <- recserar(cbind(1-b[1]-b[2] + b[1]*trimr(c(0.0, v^2),0,1)),cbind(sd(y)^2),cbind(b[2]))   
  
  nu <- abs(b[3])                                             
  z  <- v/sqrt(s2)                                                            
  
  const <- gamma( (nu+1)/2 ) / ( sqrt(pi*(nu-2)) * gamma( nu/2 ) )           
  loglt <- log(const) - 0.5*log(s2) - 0.5*(nu+1)*log( 1 + (z^2)/(nu-2) )
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
# ------------------------ GARCH Studt model --------------------------------
#
garch_studt <- function( ) {
  t  <- 10000 
  z  <- rnorm(t)
  u  <- z                         
  y  <- z
  h  <- rep(1, t)
  
  a1 <- 0.1
  b1 <- 0.8
  nu <- c(5,8,12)
  
  # Make only one call to tinv and rand
  w  <- qt(array(runif(t*3), c(t,3)), df=array(1, c(t,3))%*% diag(nu))

  # Parameter values first model
  for (k in seq(nu)) {
    cat('\n')
    cat('\nDegrees of freedom: nu  = ',nu[k])
    # Generate data 
    for (i in 2:t) {      
      z[i] <- sqrt((nu[k]-2)/nu[k])*w[i,k]                # Student t disturbance
      h[i] <- 1-a1-b1 + a1*u[i-1]^2 + b1*h[i-1]         # Conditional variance  
      u[i] <- z[i]*sqrt( h[i] )                         # Disturbance term     
      y[i] <- u[i]                                      # Returns
    }
    
    # Efficiency based on analytical derivatives Engle and Gonzalez-Rivera (JBES, 1991)
    ds2a <- recserar( cbind(c(0,trimr(y^2,0,1)) - 1), cbind(0.0) , cbind(b1) )
    ds2b <- recserar( cbind(c(0, trimr(h,0,1)) - 1), cbind(0.0) , cbind(b1) )
    ds2  <- cbind(ds2a,   ds2b)
    
    tmp   <- ds2/h
    tmp1  <- t(tmp) %*% tmp/t
    cov_h <- inv( 0.5*tmp1 )
    beta2 <- 3 + 6/(nu[k]-4)
    cov_j <- inv( (beta2-1)*0.25*tmp1 )
    
    # Analytical qmle covariance matrix of the misspecified model 
    cov_qmle <- cov_h %*% inv(cov_j) %*% cov_h                      
    
    # Analytical covariance matrix of the true model
    tmp2  <- ds2 * (1 - ((nu[k]+1)*y^2)/( y^2 + h*(nu[k]-2) ) )
    dl    <- -0.5*tmp2/h
    j0    <- t(dl) %*% dl/t
    cov_0 <- inv(j0)                                                      
    
    cat('\n\nAnalytical results\n')
    cat( diag(cov_0)/ diag(cov_qmle) )
    
    
    # GARCH(1,1) t-distribution: OPG covariance    
    start <- c(a1, b1, nu[k])
    estResults <- optim(start, lnl, y=y, method="BFGS", hessian=T)
    theta <- estResults$par
    hess <- estResults$hess
    
    vch <- inv(hess)
    tt  <- diag(vch)
    
                   
    # GARCH(1,1) normal distribution: qmle covariance
    start  <- c(a1, b1)
    estResults <- optim(start, lnl1, y=y, method="BFGS", hessian=T)
    theta0 <- estResults$par
    h0 <- estResults$hess
    
    ih   <- inv(h0)
    g0   <- numgrad( lnlt1,theta0,y)
    j0   <- t(g0) %*% g0/t
    vc0  <- ih %*% j0 %*% ih
    tt0  <- diag(vc0)
    
    # Efficiency ratio based on numerical derivatives
    cat('\n\nNumerical results\n')
    cat( (tt[1:2]/tt0) )
  }                                    
}


