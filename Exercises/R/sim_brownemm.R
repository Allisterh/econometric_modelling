# ===========================================================================
#
#    This program performs a Monte Carlo analysis to investigate the
#    sampling properties of the EMM estimator of Brownian motion with drift.
#
#    Gourieroux et. al. (1993) J of Appl. Eco.
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, recserar, repmat, zeros, inv
source("EMTSUtil.R")
library("matlab")

#----------------------------------------------------------------------------
#     The objective function to compute the emm estimator  
#----------------------------------------------------------------------------

q <- function(b,dt,v,y0,bhat,n,iinv) {
  b[2]<- abs(b[2])   
  
  # Generate continuous time data and choose every t-th observation  
  ys  <- recserar( cbind(b[1]*dt + sqrt(b[2])*v) , cbind(y0) , cbind(1) )                                               
  nys <- zeros(n*dt,1)
  
  for  (i in seq(n*dt)) {
    temp<- i*(1/dt) 
    nys[i] <- ys[temp]    
  }
  dys <- trimr(nys,1,0) - trimr(nys,0,1)
  g1s <- (dys - bhat[1])/bhat[2]
  g2s <- ((dys - bhat[1])^2/bhat[2] - 1)*0.5/bhat[2]
  gs  <- colMeans(cbind(g1s, g2s))
  val <- t(gs) %*% iinv %*% gs  
}

#
#------------------------- Brownian Motion (emm) -----------------------
#

sim_brownemm <- function() {
  t     <- 500            #     Sample size                         
  mu    <- 0.5
  sig2  <- 0.25
  theta <- c(mu, sig2)      #     Parameters of Brownian motion       
  y0    <- 1              #     Initial value of Brownian process   
  dt    <- 0.1            #     Continuous time step interval       
  h     <- 1/dt           #     Scalar to control the simulation run
  n     <- t*h            #     Length of the simulated series      


  #    Main DO LOOP to generate sampling disribution
  ndraws <- 200
  #ndraws <- 1000
  theta_emm  <- zeros(ndraws,2)

  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (j in seq(ndraws)) {
    #  Simulate Brownian motion using an exact discretisation      
    u  <- rnorm(t)
    y  <- recserar(cbind(mu + sqrt(sig2)*u) , cbind(y0) , cbind(1.0) )                            

    # Estimate the auxiliary model using actual data    
    dy   <- trimr(y,1,0) - trimr(y,0,1)                
    bhat <- c(mean(dy),  mean((dy - mean(dy))^2))
    g1   <- (dy - bhat[1])/bhat[2]
    g2   <- ((dy - bhat[1])^2/bhat[2] - 1)*0.5/bhat[2]
    g    <- cbind(g1, g2)
    

    # Compute the optimal weighting matrix           
    i <- t(g) %*% g
    p <- nrow(g)-1
    
    l <- 1
    while (l < p) {
      gam <- t( g[(l+1):nrow(g), ] ) %*% g[1:(nrow(g)-l), ]    
      i   <- i + (1.0 - l/(p+1))*(gam + t(gam))
      l   <- l + 1        
    }    
    i    <- i/nrow(g)  							
    iinv <- inv(i)
    # Generate errors to be used to compute the emm estimator.
    # Note that these errors need to be generated outside of the procedure.
    v <- sqrt(dt)*rnorm(n)
    
    # Estimate model  
    bstart   <- c(mu,  sig2)
    estResults <- optim(bstart, q, dt=dt, v=v, y0=y0, bhat=bhat, n=n, iinv=iinv, method="BFGS")
    b <- estResults$par
    theta_emm[j,] <- b   
    setTxtProgressBar(pb, j)
  }
  close(pb)

  # Generate statistics on the sampling distribution    
  mean_ <- colMeans(theta_emm)
  stdev <- sqrt(colMeans((theta_emm - repmat(mean_,nrow(theta_emm),1))^2))
  rmse  <- sqrt(colMeans((theta_emm - repmat(theta,nrow(theta_emm),1))^2))

  cat('\n')
  cat('\nNumber of replications              = ', ndraws)
  cat('\nSample size                         = ', t)
  cat('\n')
  cat('\nTrue population parameter           = ', theta)
  cat('\nMean of estimates                   = ', mean_)
  cat('\nStandard deviation of estimates     = ', stdev)
  cat('\nRMSE of Theta                       = ', rmse)  
}
