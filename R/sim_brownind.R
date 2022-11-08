# ===========================================================================
#
#      This program performs a Monte Carlo analysis to investigate the
#      sampling properties of the Indirect estimator of Brownian motion with 
#      drift.
#
#      Gourieroux et. al. (1993) J of Appl. Eco.
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, recserar, repmat
source("EMTSUtil.R")
library("matlab")

#----------------------------------------------------------------------------
#     The objective function to compute the indirect estimator  
#----------------------------------------------------------------------------
q <- function(b,dt,v,y0,bhat,n) {
  b[2]<- abs(b[2])                                
  
  # Generate continuous time data and choose every t-th observation  
  ys  <- recserar( cbind(b[1]*dt + sqrt(b[2])*v) , cbind(y0) , cbind(1) )                                               
  nys <- rep(0, n*dt)
       
  for (i in seq(n*dt)) {
    temp   <- i*(1/dt) 
    nys[i] <- ys[temp]    
  }       
  ys    <- nys
  dys   <- trimr(ys,1,0) - trimr(ys,0,1)
  bhats <- c(mean(dys),  mean((dys - mean(dys))^2))
  val   <- t( (bhat - bhats) ) %*% (bhat - bhats)
  return(val)  
}

#
#------------------------- Brownian Motion (indirect) -----------------------
#

sim_brownind <- function() {
  t     <- 500            #     Sample size                         
  mu    <- 0.5
  sig2  <- 0.25
  theta <- c(mu, sig2)      #     Parameters of Brownian motion       
  y0    <- 1              #     Initial value of Brownian process   
  dt    <- 0.1            #     Continuous time step interval       
  h     <- 1/dt           #     Scalar to control the simulation run
  n     <- t*h            #     Length of the simulated series      


  # Main DO LOOP to generate sampling disribution
  ndraws     <- 200
  #ndraws     <- 1000
  theta_ind  <- array(0, c(ndraws,2))
  pb <- txtProgressBar(min=0, max=ndraws, style=3)

  for (j in seq(ndraws)) { 
    #  Simulate sample data using an exact discretisation      
    u  <- rnorm(t)
    y  <- recserar( cbind(mu + sqrt(sig2)*u) , cbind(y0) , cbind(1.0) )                            
    
    # Estimate the auxiliary model using actual data    
    dy   <- trimr(y,1,0) - trimr(y,0,1)                
    bhat <- c(mean(dy),  mean((dy - mean(dy))^2))
    
    # Generate errors to be used to compute the emm estimator.
    # Note that these errors need to be generated outside of the procedure.
    v <- sqrt(dt)*rnorm(n)
    
    # Estimate model  
    bstart  <- c(mu,  sig2)
    estResults <- optim(bstart, q, dt=dt, v=v, y0=y0, bhat=bhat, n=n, method="BFGS")
    b <- estResults$par
    
    theta_ind[j,] <- b    
    setTxtProgressBar(pb, j)    
  } 
  close(pb)

  # Generate statistics on the sampling distribution    
  mean_ <- colMeans(theta_ind)
  stdev <- sqrt(colMeans((theta_ind - repmat(mean_,nrow(theta_ind),1))^2))
  rmse  <- sqrt(colMeans((theta_ind - repmat(theta,nrow(theta_ind),1))^2))

  cat('\n')
  cat('\nNumber of replications              = ', ndraws)
  cat('\nSample size                         = ', t)
  cat('\n ')
  cat('\nTrue population parameter           = ', theta)
  cat('\nMean of estimates                   = ', mean_)  
  cat('\nStandard deviation of estimates     = ', stdev)
  cat('\nRMSE of Theta                       = ', rmse)
}
