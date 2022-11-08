# ===========================================================================
#
#   Monte Carlo analysis to investigate the sampling properties
#   of the indirect estimator of geometric Brownian motion.
#
#   Gourieroux et. al. (1993) J of Appl. Eco.
# ===========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - trimr, recserar
source("EMTSUtil.R")

#Load required library - repmat
library("matlab")

#----------------------------------------------------------------------------
#     The objective function to compute the indirect estimator  
#----------------------------------------------------------------------------
fobj <- function(b,bhat,h,dt,dw,y0,t) {
  btilda <- array(0, c(length(bhat),h))
  ysim   <- exp(recserar( cbind((b[1] - 0.5*b[2]^2)*dt + b[2]*dw),
                          cbind(log(y0)*rep(1, h)),
                          cbind(rep(1,h))))
  
  # Generate discrete time data by choosing every t-th observation  
  nys  <- rep(0, t)
  for  (i in seq(t)) {
     temp   <- i*(1/dt)
     nys[i] <- ysim[temp]
  }
  yt <- nys
  btilda[1,] <- mean(trimr(yt,1,0)/trimr(yt,0,1)) - 1
  btilda[2,] <- sqrt(mean((trimr(yt,1,0)/trimr(yt,0,1) - 1 - btilda[1,])^2))
  w    <- diag(length(bhat))
  retp <- t( (bhat-colMeans(btilda)) ) %*% w %*% (bhat - colMeans(btilda))
  return(retp)
}
        

#
# ------------------------- Geometric Brownian Motion -----------------------
#
sim_geobrind <- function ()
{
  t     <- 100            #  Sample size                 
  h     <- 1
  mu    <- 1.0
  sig   <- 0.5
  theta <- c(mu,sig)
  y0    <- 1
  dt    <- 0.01

  #    Main DO LOOP to generate sampling disribution
  ndraws <- 100
  bmsm   <- array(0, c(ndraws,2))
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (j in seq(ndraws)) {    
    # Generate the actual data                        
    yact <- exp(recserar(cbind(mu - 0.5*sig^2 + sig*rnorm(t)),cbind(log(y0)),cbind(1.0)))
    
    # Estimate the auxiliary model using actual data   
    bhat <- rep(0, 2)
    bhat[1] <- mean(trimr(yact,1,0)/trimr(yact,0,1)) - 1
    bhat[2] <- sqrt(mean((trimr(yact,1,0)/trimr(yact,0,1) - 1 - bhat[1])^2))
    
    # Generate errors to simulate the true model      
    dw <- matrix(rnorm(t/dt*h), ncol=h)*sqrt(dt)
    
    # Estimate the indirect model using simulated data   
    b0  <- c(mu,sig)
    estResults <- optim(b0, fobj, bhat=bhat, h=h, dt=dt, dw=dw, y0=y0, t=t, method="BFGS")
    b <- estResults$par
    bmsm[j,] <- b      
    setTxtProgressBar(pb, j)
  }
  close(pb)

  # Generate statistics on the sampling distribution    
  m     <- colMeans(bmsm)
  stdev <- sqrt(colMeans((bmsm - repmat(m,nrow(bmsm),1))^2))
  rmse  <- sqrt(colMeans((bmsm - repmat(theta,nrow(bmsm),1))^2))

  cat('\n')
  cat('\nNumber of replications              = ', ndraws)
  cat('\nSample size                         = ', t)
  cat('\n')
  cat('\nTrue population parameter           = ', theta)
  cat('\nMean of estimates                   = ', m)
  cat('\nStandard deviation of estimates     = ', stdev)
  cat('\nRMSE of Theta                       = ', rmse)
}
