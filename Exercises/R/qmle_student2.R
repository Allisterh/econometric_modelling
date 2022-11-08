#================================================================================
#
#   Program to simulate data from a t - distribution and compute
#   scaled covariance matrices for true and misspecified models
#
#================================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")

#
#------------------------- Helper Functions--------------------------------------
#
#load required functions - inv
source("EMTSUtil.R")

#--------------------------------------------------------------------------------
#   Wrapper function      
#--------------------------------------------------------------------------------

lnl <- function( b,y ) {
  logl <- -mean( lnlt( b,y ) )
  return(logl)
  
}

#--------------------------------------------------------------------------------
#   Log-likelihood funciton for a Student t disturbance     
#--------------------------------------------------------------------------------

lnlt <- function( b,y ) {
  u   <- y - b[1]         
  s2  <- abs( b[2] )      
  gam <- abs( b[3] )                                                   

  z     <- u/sqrt(s2)                                                 
  const <- gamma( (gam+1)/2 ) / ( sqrt(pi*(gam-2)) * gamma( gam/2 ) )  
  loglt   <- log(const) - 0.5*log(s2) - 0.5*(gam+1)*log( 1 + (z^2)/(gam-2) )
  return(loglt)  
}

#
#--------------------  QMLE Student t Distribution ----------------------------
#

qmle_student2 <- function() {
  mu  <- 10     # Population mean      
  sig <- 1     # Population standard deviation  
  gam <- 10    # Degrees of freedom
  t   <- 5000

  # Generate data from the true model (Student t)
  v <- qt( runif(t),gam )     
  y <- mu + sig*sqrt( (gam-2)/gam )*v

  # ------------------ CORRECT SPECIFICATION ------------------------
  # Estimate the Student t model  
  bstart <- c(mu,  sig^2,  gam)

  estResults <- optim(bstart, lnl, y=y, method="BFGS")
  bhat <- estResults$par
  

  # Compute gradient and Hessian
  G <- numgrad( lnlt,bhat,y )
  H <- numhess( lnl,bhat,y )  

  iH <- inv( t*H)
    
  cat('\nResults based on the correct specification')
  cat('\n------------------------------------------')
    
  cat('\nCovariance matrix (Hessian)')
  cat('\n---------------------------\n')
  print( t*iH )

    
  cat('\nCovariance matrix (OPG)')
  cat('\n---------------------------\n')
  print( t*inv( t(G) %*% G ) )

    
  cat('\nCovariance matrix (QMLE)')
  cat('\n---------------------------\n')
  print( t*( iH %*% t(G) %*% G %*% iH ) )

  # ------------------ InCORRECT SPECIFICATION ------------------------
  # Idea here is that the QMLE estimates should approximate the 
  # correct estimates for the mean and the variance. So compare this 
  # (2x2) matrix with the upper (2x2) submatrix of the correct model.

  # Estimate the parameters of the normal 
 
  m  <- mean(y)
  s2 <- mean( (y - m)^2 )
 
  # Compute gradients of the misspecified model (normal)
  g1 <- (y - m)/s2
  g2 <- -0.5/s2 + 0.5*(y - m)^2/s2^2
  G  <- cbind(g1, g2)
  J  <- t(G) %*% G

  # Compute hessian of the misspecified model (normal) 
  H <- array(0, c(2,2 ))
  H[1,1] <- -t/s2
  H[1,2] <- -sum( y - m )/s2^2
  H[2,1] <- H[1,2]
  H[2,2] <- t*0.5/s2^2 - sum( (y - m)^2 )/s2^3
 
  I <- -H

  cat('\nResults based on the incorrect specification')
  cat('\n------------------------------------------')
    
  cat('\nCovariance matrix (Hessian)')
  cat('\n---------------------------\n')
  print( t*inv( I ) )

  cat('\nCovariance matrix (OPG)')
  cat('\n---------------------------\n')
  print( unname(t*inv( t(G) %*% G ) ) )

    
  cat('\nCovariance matrix (QMLE)')
  cat('\n---------------------------\n')
  print(  t*( inv(I) %*% t(G) %*% G %*% inv(I) ) )
  
}

  
