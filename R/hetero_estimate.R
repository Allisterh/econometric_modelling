#==============================================================================
#    Estimating a model of heteroskedasticity 
#==============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")

#
#------------------------- Helper Functions------------------------------------
#
#load required functions - inv
source("EMTSUtil.R")

#-----------------------------------------------------------------------------
#   Simulate the data  
#-----------------------------------------------------------------------------
simulatedata <- function(t) {
  beta0 <- 1
  beta1 <- 2
  gam0  <- 0.1
  gam1  <- 0.1

  x <- rnorm(t)                                            
  w <- seq(0.0,0.1*(t-1), 0.1)                                               

  u <- sqrt(exp(gam0 + gam1*w))*rnorm(t)    #  ut is N(0,sig2)                      
  y <- beta0 + beta1*x + u                   #  yt
  
  return(list(y=y, x=x, w=w))
}

#-----------------------------------------------------------------------------
# Negative unconstrained log-likelihood  
#-----------------------------------------------------------------------------
neglog1 <- function(theta,y,x,w) {
  lf <- -mean( lnlt1(theta,y,x,w) )
  return(lf)
}

#-----------------------------------------------------------------------------
# Unconstrained log-likelihood function at each observation
#-----------------------------------------------------------------------------
lnlt1 <- function(b,y,x,w) {
  mu  <- b[1] + b[2]*x
  sig <- sqrt( exp(b[3] + b[4]*w) )
  lnl <- -(1/2)*log(2*pi*sig^2) - (y - mu)^2 /(2*sig^2)
  return(lnl)
}

#-----------------------------------------------------------------------------
# Negative constrained log-likelihood  
#-----------------------------------------------------------------------------
neglog0 <- function(theta,y,x,w) {
  lf <- -mean( lnlt0(theta,y,x,w) )
  return(lf)
}

#-----------------------------------------------------------------------------
# Constrained log-likelihood function at each observation
#-----------------------------------------------------------------------------
lnlt0 <- function(b,y,x,w) {
  mu   <- b[1] + b[2]*x
  sig  <- sqrt( exp(b[3] + 0*w) )
  lnl  <- -(1/2)*log(2*pi*sig^2) - (y - mu)^2 /(2*sig^2)
  return(lnl)
}

#
#------------------------- Model Estimation ---------------------------------
#

hetero_estimate <- function() {
  
  # Simulate the model    
  make.new.random <- TRUE
  if(make.new.random) {
    t <- 500
    simResults <- simulatedata(t)
    y <- simResults$y
    x <- simResults$x    
    w <- simResults$w   
  } else {
    simResults <- read.table("estimate.dat")
    y <- as.matrix(simResults[,1])
    x <- as.matrix(simResults[,2])
    w <- as.matrix(simResults[,3])
    t <- nrow(x)    
  }  
  
  # Estimate the unconstrained model by MLE 
  theta <- c(1, 2, 0.1, 0.1)  
  estResults <- optim(theta, neglog1, y=y, x=x, w=w, method="BFGS", hessian=T)
  theta1 <- estResults$par
  lf1 <- estResults$value
  H1 <- estResults$hessian

  cat('\nLog-likelihood (unconstrained) = ', -lf1)
  cat('\nUnconstrained parameter estimates\n' )
  sterr <- (1/t)*sqrt( diag(inv(H1) ))
  print( cbind(theta1, sterr) ) 
  
  # Estimate the constrained model by MLE 
  theta <- c(1, 2, 0.1)
  estResults <- optim(theta, neglog0, y=y, x=x, w=w, method="BFGS", hessian=T)
  theta0 <- estResults$par
  lf0 <- estResults$value
  H0 <- estResults$hessian
  
  cat('\nLog-likelihood (constrained) = ', -lf0)
  cat('\nConstrained parameter estimates\n' )
  sterr <- (1/t)*sqrt( diag(inv(H0) ))
  print( cbind(theta0, sterr) )  
}





