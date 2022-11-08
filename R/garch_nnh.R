#============================================================================
#
#   Estimating ARCH-NNH Models by ML (Han and Park, J. Ects. (2008))
#
#============================================================================

rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions  ------------------------------
# 
# Load require functions - trimr, figure
source("EMTSUtil.R")


#----------------------------------------------------------------------------
# Likelihood function for a GARCH(1,1) model
#----------------------------------------------------------------------------
neglog <- function( b,y ) {
  b <- abs(b)  
  t <- length(y)
  u = y - b[1]
  h <- sd(y)*rep(1, t)
  
  for(i in 2:t) {
    h[i] <- b[2] + b[3]*u[i-1]^2 + b[4]*h[i-1]
  }
  f  <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u/sqrt( h ))^2
  lf <- -mean( f )  
  return(lf)
}
#----------------------------------------------------------------------------
# Likelihood function for a GARCH-NNH(1,1) model
#----------------------------------------------------------------------------
neglognnh <- function( b,y,w,flag,restrict ) {
  b <- abs(b)  
  t <- length( y )
  u <- y-b[1] 
  h <- sd( y )^2*rep(1, t)
  
  if (flag == 'ARCH') {
    for (i in 2:t) {
      if (restrict) 
        h[i] <- b[2]*u[i-1]^2 + b[3]*abs(w[i])
      else
        h[i] <- b[2]*u[i-1]^2 + b[3]*abs(w[i])^b[4]      
    }    
  }else if (flag =='GARCH') {
    for (i in 2:t) {
      if (restrict)
        h[i] <- b[2]+b[3]*u[i-1]^2 + b[4]*h[i-1]+ b[5]*abs(w[i])
      else
        h[i] <- b[2]+b[3]*u[i-1]^2 + b[4]*h[i-1]+ b[5]*abs(w[i])^b[6]      
    }    
  }         
  else
    cat('\nWrong flag in call to GARCH-NNH')  
  f  <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u/sqrt( h ))^2
  lf <- -mean( f )
  return(lf)  
}




#
#--------------------------- ARCH-NNH Model  --------------------------------
#

garch_nnh <- function( ) {
  # Frequency
  freq <- 'daily' # 'daily','weekly','monthly'
  
  # Load data  
  load("han_park.Rdata")
  
  if (freq == 'daily') {
    spread <- dailydata[,3]
    rets   <- dailydata[,4]
  }else if (freq =='weekly') {
    spread <- weeklydata[,3]
    rets   <- weeklydata[,4]
  }else if (freq == 'monthly'){
    spread <- monthlydata[,3]
    rets   <- monthlydata[,4]  
  } else
    cat('\nWrong frequency....')
  
  # Define variables  
  y <- trimr(rets,1,0)           
  w <- trimr(spread,0,1)   # Lag spreads (not in percentage)  
  t <- length(y)
  
  # Estimate the GARCH model
  start         <- c(0.0006, 0.001, 0.04, 0.9)
  estResults <- optim(start, neglog, y=y, method="BFGS")
  theta <- estResults$par
  lf1 <- estResults$val
  
  lf1   <- -lf1
  theta <- abs(theta)
  
  cat('\n ')
  cat('\nGARCH(1,1) Results')
  cat('\nmu                                      = ',theta[1])
  cat('\nalpha_0                                 = ',theta[2])
  cat('\nalpha_1                                 = ',theta[3])
  cat('\nbeta_1                                  = ',theta[4])
  cat('\nLog-likelihood function (unrestricted)  = ',lf1)   
  
  # Estimate the GARCH-NNH model 
  flag <- 'GARCH' # 'ARCH' or 'GARCH'
  
  if (flag =='ARCH') {
    start         <- c(0.05, 0.2, 0.1, 1.0)
    estResults <- optim(start, neglognnh, y=y, w=w, flag=flag, restrict=F, method="Nelder-Mead")
    theta1 <- estResults$par
    lfn1 <- estResults$val
    
    lfn1   <- -lfn1
    theta1 <- abs(theta1) 
    cat('\n ')
    cat('\nARCH(1)-NNH Results')
    cat('\nmu                                      = ',theta1[1])  
    cat('\nalpha                                   = ',theta1[2]) 
    cat('\nlambda                                  = ',theta1[3])
    cat('\nphi                                     = ',theta1[4])
    
    start         <- c(0.05, 0.2, 0.1)
    estResults <- optim(start, neglognnh, y=y, w=w, flag=flag, restrict=T, method="Nelder-Mead")  
    lfn0 <- estResults$val  
    
    lfn0 <- -lfn0
    
    # LR test
    cat('\nLog-likelihood function (unrestricted)  = ',lfn1)
    cat('\nLog-likelihood function (restricted)    = ',lfn0) 
    
    lr <- -2*t*(lfn0 - lfn1)
    cat('\nLR test        = ',lr)
    cat('\np-value        = ',1-pchisq(lr,1))
  } else if(flag == 'GARCH') {
    start         <- c(0.05, 0.00, 0.2, 0.5, 0.1, 1.0)
    estResults <- optim(start, neglognnh, y=y, w=w, flag=flag, restrict=F, method="BFGS")
    theta1 <- estResults$par
    lfn1 <- estResults$val
    
    lfn1 <- -lfn1
    theta1 <- abs(theta1)
    
    cat('\n ')
    cat('\nGARCH(1)-NNH Results')
    cat('\nmu                                      = ',theta1[1]) 
    cat('\nalpha_0                                 = ',theta1[2])
    cat('\nalpha_1                                 = ',theta1[3])
    cat('\nbeta_1                                  = ',theta1[4])
    cat('\nlambda                                  = ',theta1[5])
    cat('\nphi                                     = ',theta1[6])
    start         <- c(0.05, 0.00, 0.2, 0.5, 0.1)
    estResults <- optim(start, neglognnh, y=y, w=w, flag=flag, restrict=T, method="BFGS")
    lfn0 <- estResults$val
    
    cat('\nLog-likelihood function (unrestricted)  = ',lfn1)
    cat('\nLog-likelihood function (restricted)    = ',lfn0) 
    
    lr <- -2*t*(lfn0 - lfn1)
    cat('\nLR test        = ',lr)
    cat('\np-value        = ',1-pchisq(lr,1))
  }  
}

