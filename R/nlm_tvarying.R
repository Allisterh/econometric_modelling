#=========================================================================
#
#   Estimate Markov Switching model (Hamilton, Econometrica, 1989, 357-384)
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr, inv
source("EMTSUtil.R")

# Load matlab function - ones, zeros, repmat
library("matlab")

#-------------------------------------------------------------------------
#  Negative log-likelihood function
#-------------------------------------------------------------------------
neglog <- function(b,y,smooth) {
  logl <- loglt(b,y,smooth)$logl
  lf <- -mean( logl )
  return(lf)  
}

#-------------------------------------------------------------------------
#  Kalman filter
#-------------------------------------------------------------------------

loglt <- function(b,y,smooth) {
  
  t <- length(y)  
  n <- 1
  logl <- zeros(t,1)
  
  # Define the matrices for the filter   
  lam <- 0.0
  r   <- (b[1]^2)
  phi <- tanh(b[2])
  q   <- 1 
  k  <- length(q)
  
  s_predict <- zeros(t,k)     #     Declare predicted factor estimates     
  s_update  <- zeros(t,k)     #     Declare updated factor estimates       
  s_smooth  <- zeros(t,k)     #     Declare smoothed factor estimates      
  p_predict <- zeros(t,k^2)   #     Declare updated factor estimates vcov  
  p_update  <- zeros(t,k^2)   #     Declare smoothed factor estimates vcov 
  p_smooth  <- zeros(t,k^2)   #     Declare smoothed factor estimates vcov   
  
  # Run through the filter and construct the likelihood 
  st <- zeros(k,1)
  pt <- eye(k)*0.1    #     Initialization based on the diffuse prior distribution   
  
  s0 <- st
  p0 <- pt
  
  for (i in 2:t) {
    # Prediction      
    st <- phi %*% s0
    pt <- phi %*% p0 %*% t(phi) + q
    
    # Observation     
    lam <- y[i-1]     
    mt  <- lam %*% st
    vt  <- lam %*% pt %*% t(lam) + r
    ut  <- y[i] - mt
    
    # Updating        
    s0 <- st + pt %*% t(lam) %*% inv(vt) %*% ut
    p0 <- pt - pt %*% t(lam) %*% inv(vt) %*% lam %*% pt
    
    # Log-likelihood  
    logl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*cbind(ut) %*% inv(vt) %*% ut
    
    if (smooth == 1) {
      s_predict[i,] <- t(st)       #     Store predicted factor estimates        
      s_update[i,]  <- t(s0)       #     Store updated factor estimates          
      p_predict[i,] <- t(c(pt))  #     Store predicted factor estimates vcov   
      p_update[i,]  <- t(c(p0))  #     Store updated factor estimates vcov  
    }    
  }   
  
  # Generate smoothed factor estimates
  if (smooth == 1) {
    s_smooth <- s_update
    p_smooth <- p_update    
    
    for (i in seq(t-1)) {
      j <- t(matrix((p_update[t-i, ]),nrow = k,ncol = k)) %*% t(phi) %*% inv( t(matrix((p_predict[t-i + 1, ]),nrow = k,ncol = k)) ) 
      s_smooth[t-i, ] <-  t(t(s_update[t-i, ]) + j*(t(s_smooth[t-i + 1, ])-t(s_predict[t-i + 1, ]))) 
      p_smooth[t-i,] <- t( c( t(matrix((p_update[t-i, ]),nrow = k,ncol = k)) ) + j %*% (t (matrix(p_smooth[t-i+1,],k,k) ) - t(matrix(p_predict[t-i+1,],k,k) ) ) %*% t(j) )
    }    
  } 
  return(list(logl=logl, s_smooth=s_smooth))
}




#
#--------------------------- Time-varying models ----------------------------
#

nlm_tvarying <- function( )
{
  # Simulate data  
  t <- 200
  y <- rep(0, t+100)
  v <- sqrt(3)*rnorm(t+100)
  
  # Parameters
  # Choose gam = 0.0 for first experiment and gam = 0.8 for the second experiment
  phi0 <- 10.0
  phi1 <- 0.6
  phi2 <- 0.3
  gam  <- 0.8            
  
  # for (k in 3:(t+100)) {
  #   y[k] <- phi0 + phi1*y[k-1] + phi2*y[k-2] + gam * y[k-1] * v[t-1] + v[k]
  # }
  
  # Use matlab values
  y <- as.matrix(read.table("nlm_tvarying.dat"))
  
  y <- trimr(y,100,0)
  
  
  # True conditional expectation
  my_true <- phi0 + phi1*trimr(y,1,1) + phi2*trimr(y,0,2)                 
  
  # Nonparametric approximation     
  yt <- trimr(y,1,0)
  x  <- trimr(y,0,1)
  xt <- x
  h  <- as.numeric( 1.06*t(sd(yt))*t^(-1/5) )
  
  
  fx  <- t( colMeans( dnorm( t( (repmat(x,1,length(x)) - repmat(xt,length(x),1))/h ) )/h ) )
  fyx <- t( colMeans( dnorm(t( (repmat(x,1,length(x)) - repmat(xt,length(x),1))/h ) ) * repmat(yt,1,length(yt)) ) )/h
  
  my_npr <- cbind(c(fyx/fx))
  my_npr <- trimr(my_npr,1,0)
  
  # Estimate the model as a Kalman filter with time-varying parameters      
  theta_0 <- c(0.1,  0.1)
  smooth  <- 0             # Toggle to compute the smoothed factor (smooth = 1).
  estResults <- optim(theta_0, neglog, y=y, smooth=smooth, method="BFGS")
  theta <- estResults$par
  
  # Compute the smoothed factor estimates using the constrained MLEs  
  smooth <- 1
  s_smooth <- loglt(theta,y,smooth)$s_smooth
  
  my_tvar <- trimr(s_smooth,1,0) * trimr(y,0,1)                   
  my_tvar <- trimr(my_tvar,1,0)
  
  # *********************************************************************
  #
  #  Graph conditional expectations
  #
  # *********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i")
  
  tt <- 1:length(my_true)
  matplot(tt, cbind(my_true, my_tvar, my_npr), type="l",
          xlab = "y",
          ylab = "Conditional mean",     
          lty = c(1,1,1),
          col = c('red', 'green', 'blue')
          )
  
  legend("topright",                       
         legend=c('True','Time Varying','Nonparametric'),
         col = c('red', 'green', 'blue'),
         lwd=c(1,1,1,1))
}




