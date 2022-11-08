#============================================================================
#
#   Estimating GARCH - Student t Model
#  
#============================================================================

rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions  ------------------------------
# 

# Load require functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Likelihood wrapper function
# This function calls the lnlt function and returns the average log-likelihood.
#---------------------------------------------------------------------------
lnl <- function( b,y ) {
  logl <- -mean( lnlt( b,y ) )
  return(logl)
}

#----------------------------------------------------------------------------
# Likelihood function for a GARCH-t(1,1) model
#----------------------------------------------------------------------------
lnlt <- function( b,y ) {
  b <- abs(b)
  t <- length( y )
  u <- y - b[1]  
  h <- sd( y )^2*rep(1, t)
  
  for (i in 2:t) {
    h[i] <- b[2] + b[3]*u[i-1]^2 + b[4]*h[i-1]
  }
  loglt <- - 0.5*log( h ) + log(gamma((b[5]+1)/2)/((pi*(b[5]-2)^0.5)*gamma(b[5]/2))) - ((b[5]+1)/2)*log(1+(u/sqrt( h ))^2/(b[5]-2))
  return(loglt)  
}

#----------------------------------------------------------------------------
# Likelihood wrapper function for restricted model
#----------------------------------------------------------------------------
lnl0 <- function( b,y ) {
  logl <- -mean( lnlt0( b,y ) )
  return(logl)
}

#----------------------------------------------------------------------------
# Likelihood function for a GARCH-t(1,1) restricted model
#----------------------------------------------------------------------------
lnlt0 <- function( b,y ) {
  b <- abs(b)
  t <- length( y )
  u <- y - b[1]  
  h <- b[2]*rep(1, t)
  
  for (i in 2:t) {
    h[i] <- b[2]
  }  
  loglt1 <- - 0.5*log( h ) + log(gamma((b[3]+1)/2)/((pi*(b[3]-2)^0.5)*gamma(b[3]/2))) - ((b[3]+1)/2)*log(1+(u/sqrt( h ))^2/(b[3]-2))
  
  return(loglt1)  
}


#
#--------------------- GARCH Model with Conditional Student t  --------------
#
garch_student <- function( ) {
  # Choose equity index and compute percentage return and subtract mean
  # ftse, dow, nk, dummy vars for day of week (tue, wed, thurs, fri, holiday_ftse, holiday_dow, holiday_nk)
  data <- as.matrix(read.table("equity.dat", 
                               col.names=c('FSTSE', 'DOW', 'NIKKEI', 
                                           'tdum',  'wdum',  'thdum',  'fdum',  'hdum_ftse',  'hdum_dow',	'hdum_nk')))
  equity <- data[,1] # ftse, dow, nikkei
  
  y       <- 100*(trimr( log( equity ),1,0 ) - trimr( log( equity ),0,1 )) 
  y       <- y - mean(y)
  t       <- length( y )

  
  # Estimating the GARCH_t(1,1) model
  start                  <- c(0.1, 0.05, 0.9, 0.1, 8)
  estResults <- optim(start, lnl, y=y, method="BFGS", hessian=T)
  theta <- estResults$par
  lf1 <- estResults$val
  hess <- estResults$hess
  
  lf1 <- -lf1
  vc  <-(1/t)*inv(hess)
  
  
  
  cat('\nGARCH(1,1) unrestricted\n') 
  print(cbind(Estimates=theta, Std.Errors=sqrt(diag(vc)) ))
  cat('\nUnconditional Variance = ',theta[2]/(1-theta[3]-theta[4]))
  cat('\nLog-likelihood (unrestricted) = ',lf1) 
  
  # Wald test of normality: Transform the parameter so under H0 it is zero
  c  <- 1/theta[5] 
  q  <- 0
  dc <- -1/theta[5]^2     # Jacobian  
  wd <- t*(c - q)*inv(dc*vc[5,5]*dc)*(c - q)
  
  cat('\nWald statistic    = ',wd)
  cat('\np-value           = ',1-pchisq(wd,1))
  
  
  
  # Estimating the GARCH_t(1,1) model restricted
  start  <- c(0.1, 0.1, 8)
  estResults <- optim(start, lnl0, y=y, method="BFGS", hessian=T)
  theta0 <- estResults$par
  lf0 <- estResults$va;
  hess <- estResults$hess
  
  lf0 <- -lf0
  vc  <- inv(hess)
  
  cat('\nGARCH(1,1) restricted\n') 
  print(cbind(Estimates=theta0, Std.Errors=sqrt(diag(vc))))
  cat('\nLog-likelihood (unrestricted) = ',lf0)
}


