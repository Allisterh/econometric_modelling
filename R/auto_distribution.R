#==========================================================================
#
#   Simulation example to reproduce the asymptotic distribution of the 
#   MLE estimator for the regression model with autocorrelation.
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions ----------------------------------
#
#load required functions - figure, inv, recserar
source("EMTSUtil.R")

#---------------------------------------------------------------------------
# Conditional log-likelihood function
#---------------------------------------------------------------------------
neglogc <- function (b,y,x) {
  beta0 <- b[1]                                              
  beta1 <- b[2]
  rho1  <- tanh(b[3])    # Stay in the unit circle    
  sig2  <- abs(b[4])     # Variance is positive
      

  u    <- y - beta0 - beta1*x
  v    <- u[-1] - rho1*u[-length(u)] 
  tmp  <- - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2
     
  lf <- -mean( tmp )
   return(lf)
}

#-------------------------------------------------------------------------
# Exact log-likelihood function
#-------------------------------------------------------------------------
negloge <- function(b,y,x) {
  beta0 <- b[1]                                              
  beta1 <- b[2]
  rho1  <- tanh(b[3])    # Stay in the unit circle    
  sig2  <- abs(b[4])     # Variance is positive
  fac   <- 1 - rho1^2

  u    <- y - beta0 - beta1*x
  v    <- u[-1] - rho1*u[-length(u)]
  tmp  <- -0.5*log(2*pi)-0.5*log(sig2)+0.5*log(fac)-0.5*(u[1])^2/(sig2/fac) 
  tmp1 <- - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2
     
  lf <- -mean( c(tmp, tmp1) )
  return (lf) 
}

#
#---------------------- Asymptotic Distribution --------------------------
#


auto_distribution <- function () {
  # Set parameter values
  beta0  <- 1.0                                           
  beta1  <- 1.0                                                     
  rho1   <- 0.6                                                       
  sig2   <- 10                 
  theta0 <- c(beta0, beta1, rho1, sig2)  # Start values

  t       <- 500                                                                     
  ndraws  <- 5000

  # Simulate a regression model with an AR(1) disturbance term     
  x <- runif(t) - 0.5  # x fixed in repeated samples
  cond   <- array(0, c(ndraws,1))
  
  pb <- txtProgressBar(min=0, max=ndraws, style = 3)
  for (k in 2:ndraws) {    
    v <- sqrt(sig2)*rnorm(t)
    u <- recserar(cbind(v) , cbind(sqrt(1/(1-rho1^2))*v[1]) , cbind(rho1))             
    y <- beta0 + beta1*x + u
    
    # Exact MLE
    # tmp   <-  optim(theta0, negloge, y=y, x=x)$par        
       
    # Conditional MLE     
    estResults <- optim(theta0, neglogc, y=y, x=x, method="BFGS")
    tmp <- estResults$par     
    cond[k] <- tmp[2] 
    setTxtProgressBar(pb, k)
  }
  close(pb)
  mse       <- mean( (cond - beta1)^2 )
  ssq       <- sum( ( x[-1] - rho1*x[-length(x)] )^2 ) 
  i_beta    <- ssq/sig2             # Inforation matrix of beta hat    
  var_beta1 <- inv(i_beta)          # Asymptotic variance              

  cat('\nSample size                       = ', t )
  cat('\nSum of squares                    = ', ssq )
  cat('\nAsymptotic variance (theoretical) = ', var_beta1)
  cat('\nAsymptotic variance (simulated)   = ', mse)
    
  #********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************
  figure()
  par(xaxs="i", yaxs="i", yaxt='n')

  z  <- ( cond - beta1 )/sqrt( mse )
  tt <- seq(-5, 5, 0.1)

  hist(z, freq=FALSE, breaks = 31,       
       ylab = expression(f(z)),
       xlab = expression(z))
  npdf <- dnorm(tt)
  lines(tt, npdf, type="l", lwd=2)
  box(bty = "l")
}
    
