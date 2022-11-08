#=========================================================================
#
#      Program to estimate a simultaneous model with vector hetero-
#      skedasticity and first order autocorrelation
#      The system is defined as yt*b + xt*a = u      
#      where
#              u = ru(-1) + v
#      and yt is a (1xn) set of dependent variables at time t
#          xt is a (1xk) set of explanatory variables at time t
#          wt is a (1xs) set of variance explanatory variables at time t
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")
#
# -------------------------- Helper Functions ---------------------------------
#
#load required functions - inv
source("EMTSUtil.R")

#------------------------------------------------------------------------------
#   Simulate the data and returns y, x, w as a list 
#------------------------------------------------------------------------------
simulatedata <- function(t) {
  beta1  <- 0.6 
  alpha1 <- 0.4
  beta2  <- 0.2 
  alpha2 <- -0.5 

  c11 <- 1.0
  c21 <- 0.5
  c22 <- 2.0

  d11  <- 0.5
  d21  <- 0.2
  d22  <- 0.2
  
  rho11 <- 0.8 
  rho12 <- 0.1
  rho21 <- -0.2 
  rho22 <- 0.6
  
  b  <-  matrix(c(1, -beta2,
                  -beta1, 1), nrow=2, byrow=T)
  a  <- matrix(c(-alpha1, 0,
                 0, -alpha2), nrow=2, byrow=T)
  c  <-  matrix(c(c11,  0,
                  c21, c22), nrow=2, byrow=T)
  d  <-  matrix(c(d11,  0,
                  d21,  d22), nrow=2, byrow=T)
  # Exogenous variables                      
  x <- cbind(10*runif(t), 3*rnorm(t))
  w <- runif(t)
  
  # Disturbances 
  zeros <- array(0, c(t,2))
  u <- zeros
  v <- zeros
  for (i in 2:t) {
    l      <- c + d * w[i]    
    v[i,] <- rnorm(2) %*% t(l)
    u[i,1] <- rho11*u[i-1,1] + rho12*u[i-1,2] + v[i,1]
    u[i,2] <- rho21*u[i-1,1] + rho22*u[i-1,2] + v[i,2]  
  }   
  # Simulate the reduced form   
  y <- zeros
  for (i in seq(t)) {
    y[i,] <- -x[i,] %*% a %*% inv(b) + u[i,] %*% inv(b)    
  }

  
  return(list(y=y, x=x, w=w))
}

#------------------------------------------------------------------------------
# Negative unconstrained log-likelihood  
#------------------------------------------------------------------------------
neglog <- function(theta,y,x,w) {
  lf <- -mean( lnlt(theta,y,x,w) )
  # print estimates to show progress
  cat('\n theta = [', theta, "], fn = ", -lf)
  return(lf)
}

#------------------------------------------------------------------------------
# Unconstrained log-likelihood function at each observation
#------------------------------------------------------------------------------
lnlt <- function(theta,y,x,w) {
  t <- nrow(y)
  n <- ncol(y)
  b     <- matrix(c(1, -theta[3],
                  -theta[1], 1), nrow=2, byrow=T)
  a     <- matrix(c(-theta[2], 0,
                    0, -theta[4]), nrow=2, byrow=T)
  c   <- matrix(c(theta[5], 0,
                  theta[7],  theta[9]), nrow=2, byrow=T)
  d   <- matrix(c(theta[6], 0,
                  theta[8],  theta[10]), nrow=2, byrow=T) 
  rho <- matrix(c(theta[11], theta[13],
                  theta[12], theta[14]), nrow=2, byrow=T)
  
  zeros <- array(0, c(t,n))
  u   <- zeros
  v <- zeros  
  lnl <- array(0, c(t,1))
  
  for (i in 2:t) {
    u[i,] <- y[i,] %*% b + x[i,] %*% a
    v[i,] <- u[i,] - u[i-1,] %*% rho
    l     <- c + d * w[i]
    V     <- l %*% t(l)    
    lnl[i] <- - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(V)) - 0.5*v[i,] %*% inv(V) %*% cbind(v[i,])
  }  
  return(lnl)  
}

hetero_general <- function() {
  # Simulate the model   
  t <- 2000
  simResults <- simulatedata(t)
  x <- simResults$x
  y <- simResults$y
  w <- simResults$w

  # Estimate the model
  theta0 <- c(0.6, 0.4, 0.2, -0.5, 1.0, 0.5, 0.5, 0.2, 2.0, 0.2, 0.8, 0.1,-0.2, 0.6)
  estResults <- optim(theta0, neglog, y=y, x=x, w=w, method="BFGS", hessian=T)
  theta <- estResults$par
  a <- estResults$value
  H <- estResults$hessian

  vcov <- inv(H)
  cat('\nLog-likelihood function     = ', -a)
  cat('\nParameter estimates and standard errors\n' )
  sterr <- sqrt(diag(vcov))
  print(cbind(theta0, theta, sterr))  

  # Wald test of no vector heteroskedasticty and no autocorrelation             
  r   <- matrix(c(0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0,
                  0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0,
                  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0,
                  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0,
                  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0,
                  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0,
                  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1), nrow=7, byrow=T)
  q   <- cbind( rep(0, nrow(r)) )
  wd  <- t* t( (r %*% theta - q) ) %*% inv(r %*% vcov %*% t(r)) %*% (r %*% theta - q)
  dof <- nrow(r)
  cat('\nWald statistic          = ', wd) 
  cat('\nDegrees of freedom      = ', dof)
  cat('\np-value                 = ',1-pchisq(wd, dof))  
}

