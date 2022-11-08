#============================================================================
#
#   Estimate an LSTAR model of Australian unemployment rate based on the 
#   simple specification to be found in Skalin and Terasvirta ( )
#
#   The delay variable is Delta_4 u
#============================================================================

rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions  ------------------------------
# 

# Load require functions - inv, trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#
#   Wrapper function
#
#----------------------------------------------------------------------------
logl <- function( p,x ) {
  f <- - mean( loglt( p,x ) )
  return(f)
}

#----------------------------------------------------------------------------
#
#   Returns concentrated log-likelihood at each observation
#
#----------------------------------------------------------------------------
loglt <- function( p,x ) {
  y    <- x[,1] 
  ylag <- x[,2]
  st   <- x[,3]
  
  # Transition function
  gam <- abs( p[length(p)-1] )
  thr <- abs(p[length(p)])
  tmp <- (st - thr)/sd(st)
  Gt  <- 1/( 1+exp( -gam*tmp ) )
  
  # Compute errors
  tmp1 <- p[1]+p[2]*ylag 
  tmp2 <- p[3]+p[4]*ylag
  ut   <- y - (tmp1+tmp2*Gt)
  
  # Concentrate sigma
  s2 <- sd(ut)     
  
  # Likelihood function
  f <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*ut^2/s2 
  return(f)
}

#
#---------------------- LSTAR Model of unemployment -------------------------
#
nlm_ozurate <- function() {
  # Load data and set up the required lags
  data <- scan("ausu.dat", quiet=T)
  
  # Generate LSTAR variables    
  dy      <- trimr(data,1,0) - trimr(data,0,1)    # dy(t)
  y_lag   <- trimr(data,0,1)                      # y(t-1)
  dy12    <- trimr(data,12,0) - trimr(data,0,12)  # y(t)-y(t-12)
  dy12lag <- trimr(dy12,0,1)
  
  # Adjust variables to have same length as dy12lag
  dy    <- trimr(dy,12,0)
  y_lag <- trimr(y_lag,12,0)
  
  # Create data matrix 
  t <- length(dy)     
  
  # Estimate the linear model
  xlin <- cbind(rep(1, t), y_lag)
  bols <- lm(dy ~ xlin - 1)$coef
  lr_linear <- -bols[1]/bols[2]
  
  # Optimization
  x       <- cbind(dy, y_lag, dy12lag)
  pstart  <- 0.1*rep(1, 6)
  
  estResults <- optim(pstart, logl, x=x, method="BFGS", hessian=T)
  phat <- estResults$par
  hess <- estResults$hess

  ih <- 1/t * inv(hess)
  se <- sqrt(diag(ih))
  
  print( cbind(phat, se=se, t_stat=phat/se))
  
  lr_low <- -phat[1]/phat[2]
  lr_high <- -(phat[1]+phat[3])/(phat[2]+phat[4])
  
  cat('\nLong-run mean unemployment in linear state')
  cat('\n', lr_linear )
  cat('\nLong-run mean unemployment in low state')
  cat('\n', lr_low )
  cat('\nLong-run mean unemployment in high state')
  cat('\n', lr_high )  
}


