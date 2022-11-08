#============================================================================
#
#   Program to investigate alternative identification strategies
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1245)

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr
source("EMTSUtil.R")
library(MASS)

#----------------------------------------------------------------------------
#  Reduced rank log-likelihood function: Model 2 and p=1
#----------------------------------------------------------------------------
neglog <- function(b,z0,z1,z2) {
  nobs <- nrow(z0)
  f   <- rep(0, nobs)
  
  m1  <- b[1] + b[4]*(z1[,1] - b[10]*z1[,3]) + b[7]*(z1[,2] - b[10]*z1[,3]) + z2[,c(2, 3, 4)] %*% b[c(11, 12, 13)]
  m2  <- b[2] + b[5]*(z1[,1] - b[10]*z1[,3]) + b[8]*(z1[,2] - b[10]*z1[,3]) + z2[,c(2, 3, 4)] %*% b[c(14, 15, 16)]
  m3  <- b[3] + b[6]*(z1[,1] - b[10]*z1[,3]) + b[9]*(z1[,2]) - b[10]*z1[,3] + z2[,c(2, 3, 4)] %*% b[c(17, 18, 19)]
  
  v1 <- z0[,1] - m1
  v2 <- z0[,2] - m2
  v3 <- z0[,3] - m3
  
  v <- cbind(v1,  v2,  v3)
  k <- nrow(v)  
  n <- ncol(v)
  
  omegav <- t(v) %*% v/k
  for (t in seq(nobs)) {
    f[t] <- -n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[t,] %*% inv(omegav) %*% cbind(v[t,])
  }  
  lf <- -mean( f )  
  return(lf)
}


# 
#------------------------- Identification of VECMs --------------------------
#
coint_ident <- function( ) {
  t <- 200 
  n <- 3
  
  # Parameters of true DGP based on Model 3   
  beta   <- matrix(c(1,  0,
                     0,   1,
                     -0.8,   -0.8), nrow=3, byrow=T)
  beta0  <- c(4, 1)
  beta1  <- c(0,  0)
  alpha  <- matrix(c(-0.2,  0,
                     0,  -0.3,
                     0.1, 0.1), nrow=3, byrow=T)
  
  alpha0 <- Null(alpha)
  delta0 <- alpha0*0.1 
  delta1 <- alpha0*0
  psi1   <- matrix(c(-0.2, 0,  0,
                     0,  0.4, 0,
                     0, 0,  0.1), nrow=3, byrow=T)
  # Simulate the model: DGP is a vecm based on Model 3  
  v <- matrix(rnorm(t*n), nrow=t, ncol=n)
  y <- array(0, c(t+2,n))
  
  for (j in 3:(t+2)) {
    u <- beta0 + beta1*j + y[j-1,] %*% beta        
    y[j,] <- y[j-1,] + t(delta0 + delta1*j + alpha %*% t(u)) + (y[j-1,] - y[j-2,]) %*% psi1 + v[j-2,]
  }
  
  y   <- trimr(y,2,0)
  
  
  # Estimate the vecm by maximum likelihood using Johansen estimator 
  p <- 2      # Number of lags in VAR      
  r <- 1      # Number of cointegrating equations    
  model <- 3   
  
  # Estimate the vecm by maximum likelihood using the iterative estimator for triangular cross-equation normalizations  **/
  dy <- trimr(y,1,0)-trimr(y,0,1) 
  z0 <- trimr(dy,p-1,0)
  z1 <- trimr(y,p-1,1)
  
  z2 <- c()
  
  for (j in 1:(p-1)) {
    z2 <- cbind(z2, trimr(dy,p-1-j,j))
  }
  
  
  # Model 3    
  z2 <- cbind(rep(1, nrow(y)-p), z2)
  nobs <- nrow(z0)
  
  theta_0 <- 0.1*rep(1, 19)
  estResults <- optim(theta_0, neglog, z0=z0, z1=z1, z2=z2, method="BFGS", hessian=T)
  theta <- estResults$par
  logl <- estResults$val
  hess <- estResults$hess
  
  
  cat('\nCointegrating vectors (triangularization)')
  cat('\nwith cross-equation restrictions\n')
  print( matrix(c(diag(2),-theta[10], -theta[10]), nrow=3, byrow=T))
  
}


