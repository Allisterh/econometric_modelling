# ===========================================================================
#
#   Estimate a multifactor model of the Stock-Watson business cycle model
#   using indirect estimation based on Gallant and Tauchen.
#   Auxiliary model is a VAR.
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - trimr, recserar, seqa, inv
source("EMTSUtil.R")

# Load required library - repmat
library("matlab")

#----------------------------------------------------------------------------
# Objective function to compute the EMM estimator 
#----------------------------------------------------------------------------
q <- function(b,bhat,etas,zs,iinv) {
  lam <- b[1:6]
  sig <- b[7:12]
  phi <- b[13]
  
  ss <- recserar( cbind(etas), cbind(0.0), cbind(phi))  
  ys <- (lam * repmat(ss, 1, 6)) + t(apply(zs, 1, "*", sig))  
  gs <- c()
  for (i in 1:6) {
    # Estimate auxiliary model evaluated at bhat  
    yvar  <- trimr(ys[,i],2,0)          
    xvar  <- cbind(rep(1, length(yvar)), trimr(ys[,i],1,1), trimr(ys[,i],0,2))
    vhats <- yvar - xvar %*% bhat[c(1:3),i]   
    
    gs    <- cbind(gs,  (c(vhats) * xvar)/bhat[4,i],  (0.5*vhats^2/bhat[4,i]^2 - 0.5/bhat[4,i]) )
  }
  
  ret <- colMeans(gs) %*% iinv %*% cbind(colMeans(gs))
  return(ret)  
}


#
# ------------------------------- Business Cycles ---------------------------
#
sim_bcycle <- function()
{
  # Parameter values
  t <- 600                                
  lam <- c(1,0.8,0.7,0.5, -0.5, -1)
  
  sig <- 0.2*seqa(1,1,6)
  phi <- 0.8 
  
  # Generate data
  eta <- rnorm(t)
  z   <- matrix(rnorm(t*6), nrow=t)

  s   <- recserar( cbind(eta), cbind(0.0), cbind(phi))
  y   <- (lam * repmat(s, 1, 6)) + t(apply(z, 1, "*", sig))

  # Estimate the auxiliary model by a sequence of ar(1) regressions  
  g <- c()
  bhat <- c()
  for (i in 1:6) {    
    yvar <- trimr(y[,i],2,0)
    xvar <- cbind(rep(1, length(yvar),1),  trimr(y[,i],1,1), trimr(y[,i],0,2) )
    b    <- lm(yvar ~ xvar - 1)$coef
    vhat <- yvar - xvar %*% b
    s2   <- mean(vhat^2)    
    bhat <- cbind(bhat, c(b, s2))
    
    # Moment conditions based on the data    
    g <- cbind(g,  cbind((c(vhat) * xvar)/s2,   (0.5*vhat^2/s2^2 - 0.5/s2)))    
  } 

  iinv <- inv(t(g) %*% g/nrow(g))
  # Simulation Estimation  
  n    <- 10*t                  # Length of simulation run 
  etas <- rnorm(n)            # Fix simulated disturbances   
  zs   <- matrix(rnorm(n*6), nrow=n)

  # Call optimiser  
  theta0 <- c(lam,  sig,  phi)
  estResults <- optim(theta0, q, bhat=bhat, etas=etas, zs=zs, iinv=iinv)
  theta <- estResults$par
  qmin <- estResults$value

  dof  <- length(bhat)-length(theta)
  jstat <- t*qmin

  cat('\nParameter estimates')
  cat('\n True   Estimated \n')
  print(cbind(theta0, theta))
  cat('\nValue of the objective function (Q) = ', qmin)
  cat('\nValue of the J-test (TQ)            = ', jstat) 
  cat('\nDegrees of freedom                  = ', dof)
  cat('\nP-value                             = ', 1-pchisq(jstat,dof))
}
