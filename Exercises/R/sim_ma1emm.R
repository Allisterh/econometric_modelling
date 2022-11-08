#============================================================================
#
#      Monte Carlo analysis to investigate the sampling properties
#      of the EMM estimator of a MA(1) model.
#
#      Gourieroux et. al. (1993) J of Appl. Eco.
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - trimr
source("EMTSUtil.R")

#Load required library - repmat
library("matlab")

#----------------------------------------------------------------------------
# Objective function to compute the EMM estimator 
#----------------------------------------------------------------------------
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
#-------------------- EMM Estimation of the AR(1) Model ---------------------
#
sim_ma1emm <- function()
{
  # Parameters of MA(1) process
  t     <- 250        # Sample size                         
  theta <- 0.5        # MA1 Population parameter            
  lag   <- 1          # Choose AR lag for auxiliary model   
  p     <- 0          # Used to construct weighting matrix

  # Simulation settings
  nreps  <- 1000

  b      <- rep(0, nreps)
  # Main DO LOOP to generate sampling disribution
  pb <- txtProgressBar(min=0, max=nreps, style=3)
  for (j in seq(nreps)) {
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
    l <- 1
    while (l <= p) {
      gam <- t( g[(l+1):nrow(g),] ) %*% g[1:(nrow(g)-l),] 
      i   <- i + (1.0 - l/(p+1))*(gam + t(gam))  			
      l   <- l + 1      
    }
    i    <- i/nrow(g)	
    iinv <- inv(i)
    
    # Compute EMM estimator (could use a line search algorithm)
    e    <- rnorm(t)    
    estResutls <- optim(theta, q, e=e, lag=lag, bhat=bhat, iinv=iinv, method="BFGS")
    b[j] <- estResutls$par    
    setTxtProgressBar(pb, j)
  }
  close(pb)

  # Generate statistics on the sampling distribution
  b     <- tanh(b)
  m     <- mean(b)
  stdev <- sd(b)
  rmse  <- sqrt(mean(b-theta)^2)

  cat('\n ')
  cat('\nNumber of replications              =  ', nreps)
  cat('\nSample size                         =  ', t, '\n')      
  print(cbind(True=theta, Mean=m, "Std err"=stdev, RMSE=rmse))
}
