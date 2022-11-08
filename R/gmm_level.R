#============================================================================
#
#      Program to estimate level effect in interest rates by GMM
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions - inv, numgrad, figure
source("EMTSUtil.R")
# Load required library - repmat
library("matlab")
#
#----------------------------------------------------------------------------
# Define the moment equations 
#----------------------------------------------------------------------------
meqn <- function(b,drt,r1t) {
   ut <- drt - b[1] - b[2]*r1t
   zt <- cbind(rep(1, length(ut)) ,r1t)
   dt <- repmat(ut,1,2)*zt
   dt <- cbind(dt,repmat((ut^2 - (b[3]^2)*r1t^(2*b[4]) ),1,2)*zt)
   return(dt) 
}

#----------------------------------------------------------------------------
# Defines the mean of the moment conditions  
#----------------------------------------------------------------------------
meaneqn <- function(b,drt,r1t){
  ret <- cbind((colMeans(meqn(b,drt,r1t))))
  return(ret) 
}

#----------------------------------------------------------------------------
# GMM objective function which also computes the optimal w   
#----------------------------------------------------------------------------   
q <- function(b,drt,r1t,lmax,flag) {
  t <- length(drt)
  d <- meqn(b,drt,r1t)
  g <- cbind(colMeans(d))
  if (flag) {
    w   <- t(d) %*% d
    tau <- 1
    while (tau <= lmax) {
      wtau <- t( d[(tau+1):nrow(d),] ) %*% d[1:(nrow(d)-tau),]
      w    <- w + (1.0-tau/(lmax+1))*(wtau + t(wtau))
      tau  <- tau + 1      
    }  
    # Store globally
    w <<- w/t    
  } else 
     w <<- diag(ncol(d))
   ret <- t(g) %*% inv(w) %*% g 
  return (ret)
}



#
#--------------- Level Effects in U.S. Interest Rates -----------------------
#
                                           
gmm_level <- function( ) {
      
  # Load data --- monthly December 1946 to February 1991                                          
  #     0,1,3,6,9 mth and 10 yrs                                        
  load('level.Rdata')
  
  #     Choose the interest rate series         
  rt  <- xdata[,3]
  drt <- trimr(rt,1,0) - trimr(rt,0,1)
  r1t <- trimr(rt,0,1)
  t   <- length(drt)

  lmax <- 5    #  Length of the Newey-West Covariance matrix estimation                            
  
  # Estimate the model to obtain consistent estimates
  flag <- FALSE  
  b0   <- c(0.1,0.1,0.1,1.0)
  estResults <- optim(b0, q, drt=drt, r1t=r1t, lmax=lmax, flag=flag, method="BFGS")
  bgmm <- estResults$par

  # Now estimate with optimal weighting matrix
  flag <- TRUE
  estResults <- optim(bgmm, q, drt=drt, r1t=r1t, lmax=lmax, flag=flag, method="BFGS")
  bgmm <- estResults$par

  # Compute optimal weigthing matrix at GMM estimates  
  obj <- q(bgmm,drt,r1t,lmax,flag) #(obj, w)

  # Compute standard errors of GMM estimates
  dg <- numgrad(meaneqn,bgmm,drt,r1t)
  v  <- t(dg) %*% inv(w) %*% dg
  cov <- inv(v)/t
  se <- sqrt(diag(cov))

  cat('\n')
  cat('\nThe value of the objective function  = ', obj)
  cat('\nJ-test                               = ', t*obj, '\n') 
  print(cbind(Estimates=bgmm, se=se, t_stats=bgmm/se ))
  cat('\nNewey-West estimator with max lag    = ', lmax) 
  cat('\n')

  # Test of gam = 0.0,0.5,1.0,1.5 
  stat <- (bgmm[4] - 0.0)^2/cov[4,4]
  cat('\nTest of (gam=0.0) = ', stat) 
  cat('\np-value           = ', 1-pchisq(stat,1))
  cat('\n')

  stat <- (bgmm[4] - 0.5)^2/cov[4,4]
  cat('\nTest of (gam=0.5) = ', stat)
  cat('\np-value           = ', 1-pchisq(stat,1))
  cat('\n')

  stat <- (bgmm[4] - 1.0)^2/cov[4,4]
  cat('\nTest of (gam=1.0) = ', stat)
  cat('\np-value           = ', 1-pchisq(stat,1))
  cat('\n')

  stat <- (bgmm[4] - 1.5)^2/cov[4,4]
  cat('\nTest of (gam=1.5) = ', stat)
  cat('\np-value           = ', 1-pchisq(stat,1))
  cat('\n')

  #**************************************************************************
  #  
  #                   Generate Graphs
  #
  #**************************************************************************
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(2,2), bty="l")

  # Plot volatility function for alternative values of gam
  tt <- seq(from=1946+12/12,by=1/12,length.out=t)

  plot(tt,drt/r1t^0.0, type="l", xlab="", ylab="",
       main = expression(paste(gamma, " = 0.0")))
  plot(tt,drt/r1t^0.5, type="l",xlab="", ylab="",
       main = expression(paste(gamma, " = 0.5")))
  plot(tt,drt/r1t^1.0, type="l", xlab="", ylab="",
       main = expression(paste(gamma, " = 1.0")))
  plot(tt,drt/r1t^1.5, type="l", xlab="", ylab="",
       main = expression(paste(gamma, " = 1.5")))
}
