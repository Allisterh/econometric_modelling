#============================================================================
#
#   Program to estimate an ARMA model
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Wrapper function for log-liklihood
#----------------------------------------------------------------------------

neglogl <- function(b,y) {  
  logl <- - mean( loglt( b,y ) )
  return(logl)
}


#----------------------------------------------------------------------------
# Log-likelihood function for an ARMA model
#----------------------------------------------------------------------------

loglt <- function(b,y) {
  t <- nrow(y)
  n <- ncol(y)
  
  v       <- rep(0, t)
  lf      <- rep(0, t-1)
  
  # First loop over MA term   
  for (i in 2:t) {
    v[i] <- y[i]-b[1]-b[2]*y[i-1]-b[3]*v[i-1]
  }     
  v <- trimr( v,1,0 )
  omegav <- t(v) %*% v/(t-1)
  
  for (i in seq(t-1)) {
     lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(omegav))- 0.5*v[i] * inv(omegav) * cbind(v[i])    
  }
  return(lf)
}


#
#----------  ARMA Models of  US Macroeconomic variables ---------------------
#
#
stsm_arma <- function() {
  # Read the data: quarterly US data from Jan-1959 to Dec-1998
  ytdata <- as.matrix(read.table("sims_data.dat"))
 
  # Define varaibles
  r    <- ytdata[,1]        
  #lex  <- log( ytdata[,2] )
  #lcp  <- log( ytdata[,3] )
  lm   <- log( ytdata[,4] )
  lp   <- log( ytdata[,5] )
  lo   <- log( ytdata[,6] )
  #sdum <- ytdata[,7:17]

  # Choose y variable  
  y <- cbind(r)                                          # Interest rate
  t <- length(y)
  
  bstart <- 0.1*rep(1, 3)
  estResults <- optim(bstart, neglogl, y=y, method="BFGS", hessian=T)
  theta <- estResults$par
  fval <- estResults$val
  H <- estResults$hessian
  
  lnl1 <- -fval
  vc   <- (1/(t-1))*inv(H)
  
  # Wald test (AR)
  r <- rbind(c(0, 1,0))
  q <- 0
  w <- t( (r %*% theta - q) ) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta - q)
  
  cat('\n ')
  cat('\nWald statistic (AR)     = ', w)
  cat('\np-value                 = ', 1-pchisq(w,1))
  cat('\n')
  
  # Wald test (MA)   
  r <- rbind(c(0, 0,1))  
  q <- 0
  w <- t( (r %*% theta - q) ) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta - q)
  cat('\n ')
  cat('\nWald statistic (MA)     = ', w) 
  cat('\np-value                 = ', 1-pchisq(w,1))
  cat('\n ')
 
  # Wald test (Joint)  
  r <- matrix(c(0, 1,  0,
                0,  0,  1), nrow=2, byrow=T)               
  q <- rbind(0,0)
              
  w <- t( (r %*% theta - q) ) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta - q)
  cat('\n ')
  cat('\nWald statistic (Joint)  = ',w) 
  cat('\np-value                 = ', 1-pchisq(w,2))
  cat('\n ')

  # Likelihood Ratio test  
  lnl0 <- -neglogl( c(mean(y), 0, 0),y ) 
  lr   <- -2*(t-1)*(lnl0 - lnl1)

  cat('\n ')
  cat('\nLR statistic (Joint)    = ', lr) 
  cat('\np-value                 = ', 1-pchisq(lr,2))
  cat('\n ')
}
 
  

  
