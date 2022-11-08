#=========================================================================
#
#      Program to compute the MIC model 
#
#=========================================================================

rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions  ------------------------------
# 

# Load require functions - inv, trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Returns the correct data
#----------------------------------------------------------------------------
getdata <- function(index) {
  data <- as.matrix(read.table("intequity.dat", col.names= c("aord", "hang", "sing"))) # aord, hang, sing  
  
  if (index == 'aust') {
    aord <- data[,'aord']
    # Daily returns
    y <- 100*(trimr(log(aord),1,0) - trimr(log(aord),0,1))
    
    # Delete zero returns corresponding to no trading
    ind <- y != 0 & !is.nan(y) & is.notinf(y)
    y   <- y[ind]
    
    # Delete outliers
    tmp <- c(1373, 1374, 3962, 3963,4145, 4146)
    y <- y[-tmp]
    
  } else if (index == 'hang') {
    hang <- data[,hang]
    # Daily returns
    y <- 100*(trimr(log(hang),1,0) - trimr(log(hang),0,1))
    
    # Delete zero returns corresponding to no trading
    ind <- y != 0 & !is.nan(y) & is.notinf(y)
    y   <- y[ind]
    
    # Delete outliers
    tmp <- c(1939, 1940, 4651, 4652)
    y <- y[-tmp]
  }else if (index =='sing') {
    sing <- data[,'sing']
    # Daily returns
    y <- 100*(trimr(log(sing),1,0) - trimr(log(sing),0,1))
    
    # Delete zero returns corresponding to no trading
    ind <- y != 0 & !is.nan(y) & is.notinf(y)
    y   <- y[ind]    
  }
}


#-------------------------------------------------------------------------
# Likelihood function for GARCH(1,1) 
#-------------------------------------------------------------------------
neglog <- function(b,y) {
  t  <- length(y)
  v  <- rep(0, t)
  s2 <- sd(y)^2*rep(1, t)
  
  for (i in 2:t) {
    s2[i] <- abs(b[5])+abs(b[6])*v[i-1]^2+abs(b[7])*v[i-1]^2*(v[i-1]<0) +abs(b[8])*s2[i-1]
    v[i]  <- y[i] - ( b[1]+b[2]*y[i-1]+b[3]*s2[i]+b[4]*s2[i]*y[i-1])
  }  
  z  <- v/sqrt(s2)
  f  <- -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z^2
  lf <- -mean( f )
  return(lf)
}

#-------------------------------------------------------------------------
# Return mulag and s2lag
#-------------------------------------------------------------------------
parms <- function(b,y) {
  t  <- length(y)
  v  <- rep(0, t)
  s2 <- sd(y)^2*rep(1, t)
  
  for (i in 2:t) {
    s2[i] <- abs(b[5])+abs(b[6])*v[i-1]^2+abs(b[7])*v[i-1]^2*(v[i-1]<0)+abs(b[8])*s2[i-1]
    v[i]  <- y[i] - ( b[1]+b[2]*y[i-1]+b[3]*s2[i]+b[4]*s2[i]*y[i-1])    
  }
  muelag <- mean(y)-mean(v)
  s2lag  <- mean(s2)
  return(list(muelag=muelag, s2lag=s2lag))
}


#
#--------------------------- Mean Impact Curve  -----------------------------
#

garch_mic <- function( ) {
  # Choose country and load data
  index <- 'sing' # 'aust', 'hang', 'sing'
  y <- getdata(index)
  t <- length(y)
  
  # Estimate the model
  start <- c(0.01, 0.1, -0.01, -0.01, 0.01, 0.05, 0.1, 0.9)
  estResults <- optim(start, neglog, y=y, method="BFGS", hessian=T)
  theta <- estResults$par
  hess <- estResults$hess
  
  vc <- (1/t)*inv(hess)
  
  theta[5:8] <- abs(theta[5:8])
  
  # Get muelag and s2lag at the optimal parameters
  results <- parms(theta,y)
  muelag <- results$muelag
  s2lag <- results$s2lag
  
  # Graph MIC and 90# confidence interval
  
  gam0   <- theta[1]
  gam1   <- theta[2]
  gam2   <- theta[3]
  gam3   <- theta[4]
  alpha0 <- theta[5]
  alpha1 <- theta[6]
  phi    <- theta[7]
  beta1  <- theta[8]
  
  # Vector of derivatives of mue with respect to the parameters in theta'
  d   <- rep(0, length(theta))
  mue <- rep(0, 21)
  se  <- rep(0, 21)
  
  for (i in 1:21) {
    v <- -11 + i
    d[1] <- 1
    d[2] <- muelag + v
    d[3] <- alpha0 + beta1*s2lag + (alpha1 + phi*(v<0))*v^2
    d[4] <- (alpha0*muelag + beta1*s2lag*muelag) + (alpha0 + beta1*s2lag)*v + (alpha1*muelag + phi*(v<0)*muelag)*v^2 + (alpha1 + phi*(v<0))*v^3
    d[5] <- (gam2 + gam3*muelag) + gam3*v
    d[6] <- (gam2 + gam3*muelag)*v^2 + gam3*v^3
    d[7] <- (gam2*(v<0) + gam3*(v<0)*muelag)*v^2
    d[8] <- gam2*s2lag + gam3*s2lag*muelag + gam3*s2lag*v
    
    mue[i] <- (gam0 + gam2*alpha0 + gam1*muelag + gam3*alpha0*muelag + gam2*beta1*s2lag + gam3*beta1*s2lag*muelag) + (gam1 + gam3*alpha0 + gam3*beta1*s2lag)*v + (gam2*alpha1 + gam2*phi*(v<0) + gam3*alpha1*muelag + gam3*phi*(v<0)*muelag)*v^2 + (gam3*alpha1 + gam3*phi*(v<0))*v^3
    se[i] <- sqrt(d %*% vc %*% cbind(d))
  }
  
  tt <- -10:10
  print(cbind(v=tt, MIC=mue, se=se))
  
  # Plot MIC
  figure()
  matplot (tt,cbind(mue, mue+1.645*se, mue-1.645*se), type="l",
           xlab='V',
           ylab='Mean Impact Curve',
           bty="l")
}


