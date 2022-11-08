#============================================================================
#
#   Program to estimate the Campbell and Shiller present value model
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr
source("EMTSUtil.R")

 
#----------------------------------------------------------------------------
# Wrapper function for unrestricted log-liklihood
#----------------------------------------------------------------------------
neglogl <- function( b,y ) {
  lf <- - mean( loglt( b,y ) )
  return(lf)  
}

#----------------------------------------------------------------------------
# Log-likelihood function for an unrestricted Campbell-Shiller model
#----------------------------------------------------------------------------

loglt <- function (b,y) {
  nobs <- nrow(y)
  n <- ncol(y)
  
  e  <- array(0, c(nobs,2))
  lf <- rep(0, nobs-1)
  
  for (t in 2:nobs) {
    e[t,1] <- y[t,1] - b[1]- b[2]*y[t-1,1] - b[3]*y[t-1,2]
    e[t,2] <- y[t,2] - b[4] - b[5]*y[t-1,1] - b[6]*y[t-1,2]    
  }
  e <- trimr(e,1,0)
  omega <- t(e) %*% e/(nobs-1)
  for (t in seq(nobs-1)) {
    lf[t] <- - n*0.5*log(2*pi) - 0.5*log(det(omega)) - 0.5*e[t,] %*% inv(omega) %*% cbind(e[t,])
  }
  return(lf)
}

    

#----------------------------------------------------------------------------
# Wrapper function for restricted log-liklihood
#----------------------------------------------------------------------------
negloglr <- function(b,y,alpha) {
  lf <- - mean( logltr(b,y,alpha) )
  return(lf)
}


#----------------------------------------------------------------------------
# Log-likelihood function for the restricted Campbell-Shiller model
#----------------------------------------------------------------------------

logltr <- function( b,y,alpha ) {
  nobs <- nrow(y)
  n <- ncol(y)
  
  e  <- array(0, c(nobs,2))
  lf <- rep(0, nobs-1)
  
  for (t in 2:nobs) {
    e[t,1] <- y[t,1] - b[1] - b[2]*y[t-1,1] - b[3]*y[t-1,2]
    e[t,2] <- y[t,2] - b[4] - (-alpha*b[2])*y[t-1,1] - (1/b[5] - alpha*b[3])*y[t-1,2]
  }  
  e <- trimr(e,1,0)
  omega <- t(e) %*% e/(nobs-1)
  
   for (t in seq(nobs-1)) {
      lf[t] <- -n*0.5*log(2*pi) - 0.5*log(det(omega)) -0.5 * e[t,] %*% inv(omega) %*% cbind(e[t,])     
   }
  return(lf)
}

#
#------------------------ Campbell-Shiller PV Model -------------------------
#
stsm_camshiller <- function() {
    # Load data from 1933:1 to 1990:12   
  #  1. Real equity price (US)
  #	2. Real dividend (US)

  load('campbell_shiller.Rdata')
  
  
  # Define variables
  s <- ytdata[,1]		
  d <- ytdata[,2]
  t <- nrow(ytdata)

  # Estimate alpha from the present value cointegrating equation
  x     <- cbind(rep(1, length(d)), d)
  b_pv  <- lm(s ~ x - 1)$coef
  alpha <- b_pv[2]

  # Construct variables for use in VAR
  y    <- cbind(100*(trimr(d,1,0) - trimr(d,0,1)),  trimr(s - x %*% b_pv,1,0))
  nobs <- nrow(y)

  # Estimate VAR(1)
  xx <- cbind(rep(1, nobs-1), trimr(y,0,1))
  b_var <- lm(trimr(y,1,0) ~ xx - 1)$coef

  estResults <- optim(c(b_var), neglogl, y=y, method="BFGS")
  theta0 <- estResults$par
  fvalu <- estResults$val

  # Estimate the restricted model
  bstart <- c(theta0[1:4],  0.95)  
  estResults <- optim(bstart, negloglr, alpha=alpha, y=y, method="BFGS")
  theta1 <- estResults$par
  fvalr <- estResults$val
  
  cat('Estimated discount rate       = ', 1/theta1[length(theta1)]-1 )
  
  # Likelihood Ratio test 
  fvalr <- -fvalr
  fvalu <- -fvalu
  lr <- -2*t*(fvalr - fvalu)
  dof <- length(theta0) - length(theta1)

  cat('\nLR test                       = ', lr)
  cat('\nNumber of degrees of freedom  = ', dof)
  cat('\np-value                       = ', 1 - pchisq(lr, 4))
  
}

