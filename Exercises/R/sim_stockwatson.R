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
  phi <- tanh(b[13])
  
  ss <- recserar( cbind(etas), cbind(0.0), cbind(phi))  
  ys <- (lam * repmat(ss, 1, 6)) + t(apply(zs, 1, "*", sig))  
  ys <-t(apply(ys, 1, "-", colMeans(ys)))
 
  gs <- c()
  for (i in 1:6) {
    # Estimate auxiliary model evaluated at bhat    
    yvar  <- trimr(ys[,i],2,0) 
    xvar  <- cbind(trimr(ys[,i],1,1), trimr(ys[,i],0,2))
    vhats <- yvar - xvar %*%  bhat[c(1:2),i]
    gs    <- cbind(gs,  (c(vhats) * xvar)/bhat[3,i],  (0.5*vhats^2/bhat[3,i]^2 - 0.5/bhat[3,i]) )
  }  
  ret <- colMeans(gs) %*% iinv %*% cbind(colMeans(gs))
  return(ret)  
}


#
# ------------------------------- Business Cycles ---------------------------
#
sim_stockwatson <- function ()
{
    # Load monthly data September 1959 to September 2009 for Australia
  load('bcycle.Rdata')
  employed <- data[,1]    # Employed     
  gdp      <- data[,2]    # GDP          
  hincome  <- data[,3]    # Household income   
  indprod  <- data[,4]    # Industrial production     
  retail   <- data[,5]    # Retail sales  
  unemp    <- data[,6]    # Unemployment rate
  index    <- data[,7]    # Coincident index 
 
  y <- log( cbind(employed, gdp, hincome, indprod, retail, 1/unemp))
  y <- 100*(trimr(y,3,0) - trimr(y,0,3))
  y <- t(apply(y, 1, "-", colMeans(y)))
  stdy <- apply(y, 2, sd)
  y <- t(apply(y, 1, "/", stdy))
  
  t <- nrow(y)
  # Estimate the auxiliary model by a sequence of ar(1) regressions  
  g <- c()
  bhat <- c()

  for (i in 1:6) {
     yvar <- trimr(y[,i],2,0)
     xvar <- cbind(trimr(y[,i],1,1),  trimr(y[,i],0,2) )
     b    <- lm(yvar ~ xvar - 1)$coef
     uhat <- yvar - xvar %*% b
     s2   <- mean(uhat^2)
     bhat <- cbind(bhat, c(b, s2) )
     
     # Moment conditions based on the data
     g <- cbind(g,  cbind((c(uhat) * xvar)/s2,   (0.5*uhat^2/s2^2 - 0.5/s2)))    
  }

  iinv <- inv(t(g) %*% g/nrow(g))

  # Simulation Estimation  
  n    <- 30*t                  # Length of simulation run 
  etas <- rnorm(n)            # Fix simulated disturbances   
  zs   <- matrix(rnorm(n*6), nrow=n)

  # Call optimiser  
  theta0       <- c(runif(12), 1)
  estResults <- optim(theta0, q, bhat=bhat, etas=etas, zs=zs, iinv=iinv)
  theta <- estResults$par
  qmin <- estResults$value
  theta[length(theta)]   <- tanh(theta[length(theta)])
  dof  <- length(bhat)-length(theta)
  jstat <- t*qmin

 
  cat('\nParameter estimates\n')
  print(cbind(theta))
  cat('\nValue of the objective function (Q) = ', qmin)
  cat('\nValue of the J-test (TQ)            = ', jstat) 
  cat('\nDegrees of freedom                  = ', dof)
  cat('\nP-value                             = ', 1-pchisq(jstat,dof))
}

