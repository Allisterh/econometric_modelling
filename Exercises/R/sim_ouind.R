#============================================================================
#
#   Monte Carlo analysis to investigate the sampling properties 
#   of the indirect estimator of Ornstein-Uhlenbech process.
#
#      Gourieroux et. al. (1993) J of Appl. Eco.
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - recserar, trimr
source("EMTSUtil.R")

#-------------------------------------------------------------------------#
# Objective function to compute the indirect estimator 
#-------------------------------------------------------------------------#
q <- function(b,dt,v,y0,bhat,n) {
   b[3] <- abs(b[3])                                                    
   ys   <- recserar(cbind(b[1]*b[2]*dt + sqrt(b[3])*v),cbind(y0),cbind(1-b[2]*dt))

   # Generate discrete time data by choosing every t-th observation  
   nys  <- rep(0, n*dt)
   for  (i in seq(n*dt)) {
     temp   <- i*(1/dt) 
     nys[i] <- ys[temp]    
   }
   
   ys <- nys
   n <- length(ys) - 1   
   xs  <- cbind(rep(1, n) ,trimr(ys,0,1))
   ys  <- trimr(ys,1,0)
   
   b   <- lm (ys ~ xs - 1)$coef
   bhats    <- rep(0, 3)
   bhats[1] <- b[1]/(1 - b[2])
   bhats[2] <- 1 - b[2]
   bhats[3] <- mean((ys - xs %*% b)^2)
 
   retp <- ( t( (bhat - bhats) ) %*% (bhat - bhats) )
   return(retp)
}    

# 
#--------------------- Ornstein-Uhlenbeck Process ---------------------------
#
sim_ouind <- function ( )
{
  # Parameters of Ornstein-Uhlenbech process
  t     <- 250                                    
  kappa <- 0.8
  alpha <- 0.1
  sig2  <- 0.06^2
                 
  y0    <- 0.1            #     Initial value of Ornstein-Uhlenbech process  
  dt    <- 0.1            #     Continuous time step interval       
  h     <- 10/dt          #     Scalar to control the simulation run
  n     <- t*h            #     Length of the simulated series      

  # Simulation settings  
  nreps <- 200
  theta <- array(0, c(nreps,3))
  
  # Main DO LOOP to generate sampling disribution
  pb <- txtProgressBar(min=0, max=nreps, style=3)
  for (j in seq(nreps)) {
    u  <- rnorm(t)
    y  <- recserar(cbind(alpha*(1-exp(-kappa)) + sqrt(sig2)*sqrt((1-exp(-2*kappa))/(2*kappa))*u),cbind(y0),cbind(exp(-kappa)))
    
    # Estimate the auxiliary model using actual data  
    x  <- cbind(rep(1, nrow(y)-1) ,trimr(y,0,1))
    y  <- trimr(y,1,0)
    b  <- lm(y ~ x - 1)$coef
    
    bhat    <- rep(0, 3)
    bhat[1] <- b[1]/(1 - b[2])
    bhat[2] <- 1 - b[2]
    bhat[3] <- mean((y - x %*% b)^2)
    
    # Compute the indirect estimator.
    v          <- sqrt(dt)*rnorm(n)
    b0         <- c(alpha,  kappa, sig2)
    estResults <- optim(b0, q, dt=dt, v=v, y0=y0, bhat=bhat, n=n, method="BFGS" )
    b <- estResults$par
    
    b[3]       <- abs(b[3])
    theta[j,] <- b   
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
   
  # Generate the actual data for the Ornstein-Uhlenbech process       
  # Generate statistics on the sampling distribution
  m     <- colMeans(theta)
  stdev <- apply(theta, 2, sd)
  diff <- t(apply(theta, 1, '-', cbind(alpha, kappa, sig2)^2))
  rmse  <- sqrt(colMeans(diff))
  
  cat('\n ')
  cat('\nNumber of replications              =  ', nreps)
  cat('\nSample size                         =  ', t, '\n')
  true <- cbind(c(alpha, kappa, sig2))
  print(cbind("True"=true, "Mean"=m, "Std err"=stdev, "RMSE"=rmse))
}

