#============================================================================
#
#   Program to estimate a one factor model of the term structure
#   Alternative formulation of the Kalman filter
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()


#
#--------------------------- Helper Functions -------------------------------
#

# Load required functions - inv
source("EMTSUtil.R")

#-------------------------------------------------------------------------
#  Log-likelihood function from Kalman filter
#-------------------------------------------------------------------------
neglog <- function ( b,y,flag ){
  # Define matrices for filter
  Lam <- cbind(b[1:6],   diag(abs(b[7:12])) )
  
  if (flag)
    Phi <- diag( c(tanh(b[13]), rep(0, 6)) )
  else
    Phi <- diag( c(b[13], rep(0, 6)) )  
  R     <- matrix(0, 6,6)
  Q     <- diag(7)     
  
  lf    <- -mean( lnlt(y,Phi,Lam,R,Q) ) 
  
  return(lf)
}

#--------------------------------------------------------------------------
# Multivariate Kalman filter
#--------------------------------------------------------------------------
lnlt <- function(y,Phi,Lam,R,Q) {
  # Allocate arrays  
  t <- nrow(y)
  n <- ncol(y)
  k         <- nrow(Q)
  lnl       <- rep(0, t)
  
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- diag(k)
  pt[1,1] <- 0.1    # Impose prior on first factor

  mt <- Lam %*% st
  vt <- Lam %*% pt %*% t(Lam) + R
  ut <- y[1,] - mt

  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
  Kal <- pt %*% t(Lam) %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt
  
  # Main loop over observations
  
  # Main loop over observations
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0     
    pt <- Phi %*% p0 %*% cbind(Phi) + Q          
    
    # Observation
    mt <- Lam %*% st
    vt <- Lam %*% pt %*% t(Lam) + R
    ut <- y[i,] - mt
    
    # Construct log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
    
    # Update    			
    Kal <- pt %*% t(Lam) %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt   
  } 
  return(lnl)
}


#
#--------------- Alternative Formulation of the Kalman Filter ---------------
#

lfac_term_adj <- function( ) {
  # Read the data:
  # US daily zero coupon yields starting 10/04/1988 and ending 12/28/2001
  #   The variables are:
  #
  #  	1.	tcm3m
  #		2.	tcm1y
  #		3.	tcm3y
  #		4.	tcm5y
  #		5.	tcm7y
  #		6.	tcm10y
  
  #load lfac_usdata.mat
  interest <- read.table("lfac_usdata.dat")
  
  #	Choose the variables				
  r <- interest[,c(6:1)]
  
  #	Rescale the variables (used for numerical precision when using mle)		
  y <- r*100   
  y <- t( apply(y, 1, '-', colMeans(y)) )
  t <- nrow(y)
  
  # Estimate parameters with restriction imposed
  flag  <- 1
  start <- c(5.5359, 
             5.8476,   
             6.2928,    
             7.0033,   
             7.6360,   
             7.3535,   
             50.4656,   
             36.6785,  
             26.5888,   
             6.2867,  
             43.6677,  
             65.7365,  
             0.9993)

  
  estResults <- optim(start, neglog, y=y, flag=flag, method="BFGS")
  bhat <- estResults$par
  
  # Restimate without restrictions to get standard errors
  flag  <- 0
  start <- c(bhat[1:12], tanh(bhat[13]))
  
  estResults <- optim(start, neglog, y=y, flag=flag, method="BFGS", hessian=T)
  bhat <- estResults$par
  lf <- estResults$val
  hess <- estResults$hess
  
  vc <- (1/t)*inv(hess)
  
  cat('\nLog-likelihood function = ',-lf)
  cat('\n' )
  print(cbind("Params"=bhat,  "Std. Errors"=sqrt(diag(vc))))
  
}

