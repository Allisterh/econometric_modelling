#============================================================================
#
#      Over-identification test of the MA(1) model using the EMM estimates
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - trimr
source("EMTSUtil.R")

#Load required library - repmat
library("matlab")

#-------------------------------------------------------------------------#
# Objective function to compute the EMM estimator 
#-------------------------------------------------------------------------#
q <- function(b,e,lag,bhat,iinv) {
  ys     <- trimr(e,1,0) - tanh(b)*trimr(e,0,1)    
  ys     <- ys - mean(ys)
  
  if (lag == 1) {
    xs    <- trimr(ys,0,1)
    ehats <- trimr(ys,1,0) - xs*bhat
    gs    <- cbind(mean( ehats*xs ))
  }
  else if (lag == 2) {
    xs    <- cbind(trimr(ys,1,1), trimr(ys,0,2))
    ehats <- trimr(ys,2,0) - xs %*% bhat
    gs    <- cbind(colMeans( repmat(ehats,1,2)*xs))
  }
  else if (lag == 3) {
    xs    <- cbind(trimr(ys,2,1), trimr(ys,1,2), trimr(ys,0,3))
    ehats <- trimr(ys,3,0) - xs %*% bhat
    gs    <- cbind(colMeans( repmat(ehats,1,3)*xs ))    
  }
  retp <- t(gs) %*% iinv %*% gs  
  return(retp)
}


#
# ------------- Over-identification Test of a Moving Average Model ----------
#
sim_ma1overid <- function ()
{
  # Parameters of MA(1) process
  t     <- 250        # Sample size                         
  theta <- 0.5        # MA1 Population parameter            
  lag   <- 3          # Choose AR lag for auxiliary model   
  p     <- 0          # Used to construct weighting matrix

  # Simulation settings 
  # Generate the actual data for the MA(1) process       
  u  <- rnorm(t)
  y  <- trimr(u,1,0) - theta*trimr(u,0,1) 
  
  
  # Estimate the the auxiliary model using actual data  
  y   <- y - mean(y)  
  if (lag == 1) {
    x    <- trimr(y,0,1)
    bhat <- lm(trimr(y,1,0) ~ x - 1)$coef
    ehat <- trimr(y,1,0) - x*bhat
    g    <- cbind(ehat*x)
  } else if (lag == 2) {
    x    <- cbind(trimr(y,1,1), trimr(y,0,2))
    bhat <- lm(trimr(y,2,0) ~ x - 1)$coef
    ehat <- trimr(y,2,0) - x %*% bhat
    g    <- repmat(ehat,1,2)*x    
  } else if (lag == 3 ) {
    x    <- cbind(trimr(y,2,1), trimr(y,1,2), trimr(y,0,3))
    bhat <- lm(trimr(y,3,0) ~ x - 1)$coef
    ehat <- trimr(y,3,0) - x %*% bhat
    g    <- repmat(ehat,1,3)*x        
  }
  
  # Compute the optimal weighting matrix
  i <- t(g) %*% g
  if (p > 0) {
    for (l in seq(p)) {
      gam <- t( g[(l+1):nrow(g),] ) %*% g[1:(nrow(g)-l),]
      i   <- i + (1.0 - l/(p+1))*(gam + t(gam))        
    }      
  }      
  i    <- i/length(g)									
  iinv <- inv(i)		

  # Compute EMM estimator (could use a line search algorithm)
  e        <- rnorm(t)
  estResults <- optim(theta, q, e=e, lag=lag, bhat=bhat, iinv=iinv, method="BFGS")
  b <- estResults$par
  qmin <- estResults$value
  
  # Number of restrictions  
  r <- length(bhat) - length(b)       

  cat('\n')
  cat('\nSample size                         =  ', t)
  cat('\nTrue parameter value                =  ', theta)
  cat('\nLags in auxiliary model             =  ', lag)
  cat('\nOver-identification test            =  ', t*qmin)
  cat('\nNumber of restrictions              =  ', r)  
  if (r > 0)
        cat('\np value                             =  ', 1-pchisq(t*qmin,r))
}
