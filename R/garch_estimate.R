#=========================================================================
#
#   Estimating GARCH Models by ML
#   This file estimates a variety of GARCH models for daily US 
#   zero-coupon yields.
#
#=========================================================================

rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions  ------------------------------
# 

# Load require functions - trimr, figure
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Wrapper function for an GARCH model
#----------------------------------------------------------------------------
neglog <- function( b,y ) {
  logl <- -mean( lnlt( b,y ) )  
  return(logl)
}

#----------------------------------------------------------------------------
# Likelihood function for a GARCH(1,1) model
#----------------------------------------------------------------------------
lnlt <- function( b,y ) {
  b <- abs(b)
  t <- length( y )
  u <- y   
  h <- sd( y )^2*rep(1, t)
  
  for (i in 2:t) {
    h[i] <- b[1] + b[2]*u[i-1]^2 + b[3]*h[i-1]
  }
  loglt <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u^2/h)
  return(loglt)  
}

#----------------------------------------------------------------------------
# Wrapper function for an ARCH model
#----------------------------------------------------------------------------
neglog0 <- function( b,y ) {
  logl <- -mean( lnlt0( b,y ) )
  return(logl)
}

#----------------------------------------------------------------------------
# Likelihood function for a ARCH(1) model
#----------------------------------------------------------------------------
lnlt0 <- function( b,y ) {
  b <- abs(b)
  t <- length( y )
  u <- y  
  h <- sd(y)^2*rep(1, t)
  
  for (i in 2:t) {
    h[i] <- b[1] + b[2]*u[i-1]^2
  }
  loglt <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u^2/h)
  return(loglt)
}

#----------------------------------------------------------------------------
# Wrapper function for an IGARCH model
#----------------------------------------------------------------------------
neglogi <- function( b,y ) {
  logl <- -mean( lnlti( b,y ) )
  return(logl)
}
#----------------------------------------------------------------------------
# Likelihood function for a IGARCH* model
#----------------------------------------------------------------------------
lnlti <- function( b,y ) {
  b <- abs(b)
  t <- length( y )
  u <- y   
  h <- sd( y )^2*rep(1, t)
  
  for (i in 2:t) {
    h[i] <- b[1] + b[2]*u[i-1]^2 + (1-b[2])*h[i-1]
  }  
  loglt <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u/sqrt( h ))^2
  return(loglt)
}

#
#--------------------------- GARCH models   ---------------------------------
#
garch_estimate <- function() {
  
  # Choose equity index and compute percentage return and subtract mean
  # ftse, dow, nk, dummy vars for day of week (tue, wed, thurs, fri, holiday_ftse, holiday_dow, holiday_nk)
  data <- as.matrix(read.table("equity.dat", 
                               col.names=c('FSTSE', 'DOW', 'NIKKEI', 
                                           'tdum',  'wdum',  'thdum',  'fdum',  'hdum_ftse',	'hdum_dow',	'hdum_nk')))
  equity <- data[,3] # ftse, dow, nikkei
  
  y       <- 100*(trimr( log( equity ),1,0 ) - trimr( log( equity ),0,1 )) 
  y       <- y - mean(y)
  t       <- length( y )
  
  # Estimating the GARCH(1,1) model
  start  <- c(0.05, 0.1, 0.9)
  estResults <- optim(start, neglog, y=y, method="BFGS", hessian=T)
  theta <- estResults$par
  lf1 <- estResults$val
  hess <- estResults$hess
  
  lf1 <- -lf1
  
  cat('\n ')
  cat('\nGARCH(1,1) Results')
  cat('\nalpha_0                                 = ',theta[1])
  cat('\nalpha_1                                 = ',theta[2])
  cat('\nbeta_1                                  = ',theta[3])
  cat('\nLog-likelihood function (unrestricted)  = ',lf1)   
  
  # Compute conditional variance at optimal parameters
  h <- sd( y )^2*rep(1,t)
  u <- y - theta[1]
  
  for (i in 2:t) {
    h[i] <- theta[1] + theta[2]*u[i-1]^2 + theta[3]*h[i-1]
  }
  
  #*********************************************************************
  #**     Generate graph of conditional variance
  #*********************************************************************
  
  figure()
  par(mar=c(5,5,5,5))
  plot(seq(t), h,type="l",
       xlab="t",
       ylab=expression(sigma[t]^2),
       bty="l")
  
  # Comparison of standard errors 
  ih  <- inv(hess)
  seh <- sqrt(diag((1/t)*ih))  
  g   <- numgrad(lnlt,theta,y)
  j   <- t(g) %*% g/t
  sej <- sqrt(diag( (1/t)*inv(j)))
  seq <- sqrt(diag( (1/t)*( ih %*% j %*% ih)))
  
  cat('\n ')
  cat('\nComparison of Standard Errors\n' )
  print( cbind(Hessian=seh, OPG=sej, QMLE=seq))
  
  # Estimating the restricted ARCH(1) model
  start         <- c(0.05, 0.1)
  estResults <- optim(start, neglog0, y=y, method="Nelder-Mead")
  theta0 <- estResults$par
  lf0 <- estResults$val
  
  lf0 <- -lf0
  cat('\n ')
  cat('\nARCH(1) Results')
  cat('\nalpha_0                                 = ',theta0[1])
  cat('\nalpha_1                                 = ',theta0[2])
  cat('\nLog-likelihood function (unrestricted)  = ',lf0)   
  
  # Likelihood Ratio test   
  lr <- -2*t*(lf0-lf1)
  pv <- 1 - pchisq(lr, 1)
  cat('\nLR  test   = ',lr)
  cat('\np-value    = ',pv)
  
  # Estimating the IGARCH
  start         <- c(0.05, 0.1)
  estResults <- optim(start, neglogi, y=y, method="BFGS")
  thetai <- estResults$par
  lfi <- estResults$val
  
  thetai <- abs(thetai)
  lfi <- -lfi
  cat('\n ')
  cat('\nIGARCH(1,1) Results')
  cat('\nalpha_0                           = ',thetai[1])
  cat('\nalpha_1                           = ',thetai[2])
  cat('\nbeta_0                            = ',1-thetai[2])
  cat('\nLog-likelihood function (IGARCH)  = ',lfi)
}

