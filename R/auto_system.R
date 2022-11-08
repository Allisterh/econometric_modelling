# ========================================================================
#
#      Program to estimate a simultaneous model with first order vector
#      autocorrelation. The set of equations is defined as yt*b + xt*a = u
#      where
#              u = ru(-1) + v
#          where yt is a (1xn) set of dependent variables at time t
#                xt is a (1xk) set of explanatory variables at time t
#
# =======================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
# -------------------------- Helper Functions --------------------------
#

#load required functions - inv, trimr
source("EMTSUtil.R")


#-----------------------------------------------------------------------
# Negative unconstrained log-likelihood  
#-----------------------------------------------------------------------
neglog1 <- function(theta,y,x) {
  lf <- -mean( lnlt1(theta,y,x) )
  return(lf)
}

#-----------------------------------------------------------------------
# Unconstrained log-likelihood function
#-----------------------------------------------------------------------
lnlt1 <- function(theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)
  b     <- matrix(c(1, -theta[3],
                  -theta[1], 1), nrow=2, byrow=T)
  a     <- matrix(c(-theta[2], 0,
                    0, -theta[4]), nrow=2, byrow=T)
  rho   <- matrix(c(theta[5], theta[7],
                    theta[6],  theta[8]), nrow=2, byrow=T)
  
  # Construct residuals and concentrate the covariance matrix  
  u <- array(0, c(t,n))
  v <- array(0, c(t,n))
  for (i in 2:t) {
    u[i,] <- y[i,] %*% b + x[i,] %*% a     
    v[i,] <- u[i,] - u[i-1,] %*% rho 
  }
  
  omega <- t(v) %*% v/t  
    
  lnl <- array(0, c(t,1))
  for (i in 2:t) {
    lnl[i] <- - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*v[i,] %*% inv(omega) %*% cbind(v[i,])
  }
  lnl <-  trimr(lnl,1,0)
  return(lnl)  
}

#-----------------------------------------------------------------------
# Negative constrained log-likelihood function
#-----------------------------------------------------------------------
neglog0 <- function(theta,y,x) {
  lf <- -mean( lnlt0(theta,y,x) )
  return(lf)
}

#-----------------------------------------------------------------------
# Constrained log-likelihood function
#-----------------------------------------------------------------------  
lnlt0 <- function(theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)  
  b     <- matrix(c(1, -theta[3],
                  -theta[1], 1), nrow=2, byrow=T)
  a     <- matrix(c(-theta[2], 0,
                    0, -theta[4]), nrow=2, byrow=T)
  rho   <- matrix(0, nrow=2, ncol=2, byrow=T)
  
  # Construct residuals and concentrate the covariance matrix 
  u <- array(0, c(t,n))
  v <- array(0, c(t,n))
  for (i in 2:t) {
    u[i,] <- y[i,] %*% b + x[i,] %*% a     
    v[i,] <- u[i,] - u[i-1,] %*% rho 
  }
  
  omega <- t(v) %*% v/t  
    
  lnl <- array(0, c(t,1))
  for (i in 2:t) {
    lnl[i] <- - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*v[i,] %*% inv(omega) %*% cbind(v[i,])
  }
  lnl <-  trimr(lnl,1,0)
  return(lnl)
}

    
#-----------------------------------------------------------------------
# Negative constrained log-likelihood function (independent auto)
#-----------------------------------------------------------------------
neglog2 <- function(theta,y,x) {
  lf <- -mean( lnlt2(theta,y,x) )
  return(lf)  
}

#-----------------------------------------------------------------------
#   Constrained log-likelihood function at each observation
#   with independent autocorrelation
#-----------------------------------------------------------------------
lnlt2 <- function(theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)
  b     <- matrix(c(1, -theta[3],
                  -theta[1], 1), nrow=2, byrow=T)
  a     <- matrix(c(-theta[2], 0,
                    0, -theta[4]), nrow=2, byrow=T)
  rho   <- matrix(c(theta[5], 0,
                    0,  theta[6]), nrow=2, byrow=T)
  
  # Construct residuals and concentrate the covariance matrix  
  u <- array(0, c(t,n))
  v <- array(0, c(t,n))
  for (i in 2:t) {
    u[i,] <- y[i,] %*% b + x[i,] %*% a     
    v[i,] <- u[i,] - u[i-1,] %*% rho 
  }
  
  omega <- t(v) %*% v/t  
    
  lnl <- array(0, c(t,1))
  for (i in 2:t) {
    lnl[i] <- - n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*v[i,] %*% inv(omega) %*% cbind(v[i,])
  }
  lnl <-  trimr(lnl,1,0)
  return(lnl)
}

#
# -------------------- Autocorrealted Systems ----------------------
#

auto_system <- function() {
  t <- 500

  # Load GAUSS data to reproduce results
  gaussdata <- read.table("system.dat")
    
  y <- as.matrix(gaussdata[, 1:2])
  x <- as.matrix(gaussdata[, 3:4])
    
  # Estimate the unconstrained model  
  theta <-  c(0.6, 0.4, 0.2, -0.5, 0.8, 0.1,-0.2, 0.6)
  estResults <- optim(theta, neglog1, y=y, x=x, method="BFGS", hessian=T)
  theta1 <- estResults$par
  a1     <- estResults$value
  H      <- estResults$hessian

  lnl1  <- -(t-1)*a1                           # Unconstrained log-likelihood                       
  vcov1 <- inv(H)

  # Estimate the constrained model 
  theta <-  c(0.6, 0.4, 0.2, -0.5)
  estResults <- optim(theta, neglog0, y=y, x=x, method="BFGS")
  theta0 <- estResults$par
  a0     <- estResults$value

  lnl0 <- -(t-1)*a0                           # Constrained log-likelihood      

  # Estimate the constrained model with independent disturbances
  theta <-  c(0.6, 0.4, 0.2, -0.5, 0.8, 0.1)
  estResults <- optim(theta, neglog2, y=y, x=x, method="BFGS")
  theta2 <- estResults$par
  a2     <- estResults$value

  lnl2 <- -(t-1)*a2                            # Constrained log-likelihood 

  cat('\n')
  cat('\nUnconstrained log-likelihood function                  = ', lnl1)
  cat('\nConstrained log-likelihood function                    = ', lnl0)
  cat('\nConstrained log-likelihood function (independent auto) = ', t*lnl2 )
  
  cat('\n')
  # LR test of no autocorrelation     
  lr  <- -2*(lnl0 - lnl1)
  dof <- length(theta1) - length(theta0)
  cat('\nLR test (no autocorrelation)       = ', lr)
  cat('\nDegrees of freedom                 = ', dof)
  cat('\np-value                            = ', 1-pchisq(lr,dof))
  
  cat('\n')
  # LR  test of independent autocorrelation
  lr  <- -2*(lnl2 - lnl1)
  dof <- length(theta1) - length(theta2)
  cat('\nLR test (indep autocorrelation)       = ', lr)
  cat('\nDegrees of freedom                    = ', dof)
  cat('\np-value                               = ', 1-pchisq(lr,dof))
	  
  #     Wald test of no autocorrelation
  r  <- matrix(c(0 , 0 , 0 , 0 , 1 , 0 , 0 , 0,
                 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0,
                 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0,
                 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1), nrow=4, byrow=T)
  q  <- cbind(c(0,0,0,0))
    
  wd <- t* t( (r %*% theta1 - q) ) %*% inv(r %*% vcov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof <- nrow(r)
  cat('\n')
  cat('\nWald test (no autocorrelation)       = ', wd)
  cat('\nNumber of degrees of freedom         = ', dof)
  cat('\np-value                              = ', 1-pchisq(wd,dof))

  #     Wald test of common autocorrelation                   
  r  <- matrix(c(0 , 0 , 0 , 0 , 1 , 0 ,  0 , -1,
                 0 , 0 , 0 , 0 , 0 , 1 , -1 ,  0), nrow=2, byrow=T)

  q  <- cbind(c(0,0))
  wd <- t* t( (r %*% theta1 - q) ) %*% inv(r %*% vcov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof <- nrow(r)

  cat('\n')
  cat('\nWald test (indep autocorrelation)    = ', wd)
  cat('\nNumber of degrees of freedom         = ', dof)
  cat('\np-value                              = ', 1-pchisq(wd,dof))


  # Lagrange Multiplier test (based on numerical opg matix)              
  dof    <- length(theta1) - length(theta0)
  theta  <- cbind(theta0,  rep(0, 4))                
  gmat   <- numgrad(lnlt1,theta,y,x)
  g      <- cbind(colMeans(gmat))
  j      <- t(gmat) %*% gmat/t
  lm     <- t* t( g ) %*% inv(j) %*% g

  cat('\n')
  cat('\nLM test (no autocorrelation)         = ', lm)
  cat('\nNumber of degrees of freedom         = ', dof)
  cat('\np-value                              = ', 1-pchisq(lm,dof))
}    



