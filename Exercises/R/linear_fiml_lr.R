#=========================================================================
#
#  Program to estimate model by full information maximum likelihood 
#  and do a LR test. 
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#------------------------ Helper Functions -------------------------------------#

# load utility functions - inv and numhess functions
source("EMTSUtil.R")

#--------------------------------------------------------------------------------
# Simulate data an return the results as a list of x and y
#--------------------------------------------------------------------------------
 
simulatedata <- function(t) {
  beta1  <- 0.6 
  beta2  <- 0.2 
  alpha1 <- 0.4
  alpha2 <- -0.5
  
  omega <- matrix(c(1, 0.5, 0.5, 1.0), nrow = 2, byrow=T)
  
  # Construct population parameter matrices     
  b <- matrix(c(1, -beta2, -beta1, 1.0), nrow = 2, byrow=T) 
  a <- matrix(c(-alpha1, 0, 0, -alpha2), nrow = 2, byrow=T)
  
  # Construct exogenous variables from a uniform distribution                     
  x <- cbind( 10*rnorm(t), 3*rnorm(t) )
  
  # Construct disturbances                   
  u <- matrix(rnorm(t*2), ncol = 2) %*% chol(omega)                   
  invB <- inv(b)
  
  # Simulate the model
  y <- array(0, c(t,2))
  
  for (i in seq(t)) {
    y[i,] <- -x[i,] %*% a %*% invB + u[i,] %*% invB
  }
  return(list(x=x, y=y))
}

#--------------------------------------------------------------------------------
# Constrained log-likelihood function
#--------------------------------------------------------------------------------
neglog0 <- function(theta,y,x){
  lf <- -mean(lnlt0(theta,y,x))
  return(lf)
}

#--------------------------------------------------------------------------------
# Constrained log-likelihood function at each observation
#-------------------------------------------------------------------------------- 
lnlt0 <- function(theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)
  
  b <- matrix(c(1, -theta[3], -theta[1], 1), nrow = 2, byrow=T)
  a <- matrix(c(-theta[2], 0, 0, theta[2]), byrow = T, nrow = 2)
  
  u <- array(0, c(t,n))  
  for (i in seq(t)) {
     u[i,] <- y[i,] %*% b + x[i,] %*% a    
  }  
  omega <- t(u) %*% u/t
  lf    <- array(0, c(t,1))
  for (i in seq(t)) {
    lf[i] <- -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u[i,] %*% inv(omega) %*% cbind(u[i,])
  }
  return(lf)
}
   

#-----------------------------------------------------------------
# Unconstrained Log-likelihood function
#-----------------------------------------------------------------
        
neglog1 <- function(theta,y,x) {
   lf <- -mean(lnlt1(theta,y,x))
   return(lf)
}

#-----------------------------------------------------------------
# Unconstrained Log-likelihood function at each observation 
#-----------------------------------------------------------------
lnlt1 <- function (theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)  
    
  b <- matrix(c(1, -theta[3], -theta[1], 1), nrow = 2, byrow=T)
  a <- matrix(c(-theta[2], 0, 0, -theta[4]), nrow = 2, byrow=T)
  u <- array(0, c(t,n))
  for (i in seq(t)) {
     u[i,] <- y[i,] %*% b + x[i,] %*% a    
  }
  omega <- t(u) %*% u/t                         # Concentrate out resid var-covar matrix  
  lnl <- array (0, c(t,1))

  for (i in seq(t)) {
    lnl[i] <- -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u[i,] %*% inv(omega) %*% cbind(u[i,])
  }
  return(lnl)  
}

#------------------------------------- FIML LR Test  -------------------------------------#
  
linear_fiml_lr <- function() {
  # Simulate
  t <- 500
  make.new.random <- FALSE
  if (make.new.random) {
    simResults <- simulatedata(t)
    x <- simResults$x
    y <- simResults$y    
  } else {
    simResults <- read.table("linear_fimltestdata.dat")
    y <- as.matrix(simResults[,1:2])
    x <- as.matrix(simResults[,3:4])
  }
  # Estimate the unconstrained model  
  start <- runif(4) 
  estResults <- optim(start, neglog1, y = y, x = x, method = "BFGS")
  theta1 <- estResults$par
  lf1    <- estResults$value
  lf1 <- -lf1
  
  u <- array(0, c(t,2))
  b0 <- matrix(c(1, -theta1[3], -theta1[1], 1), nrow = 2, byrow = T)
  a0 <- matrix(c(-theta1[2], 0, 0, -theta1[4]), nrow = 2, byrow = T)
  for (i in seq(t)) {
    u[i,] <- y[i,] %*% b0 + x[i,] %*% a0
  }
  cat('\nResidual variance-covariance matrix (unrestricted)\n')
  omega1 <- t(u) %*% u/t
  print(omega1)

  
  # Estimate the constrained model  
  start <- runif(3)
  estResults <- optim(start, neglog0, y = y, x = x, method = "BFGS")
  theta0 <- estResults$par
  lf0 <- estResults$value
  lf0 <- -lf0  

  b0 <- matrix(c(1, -theta0[3], -theta0[1], 1), nrow = 2, byrow = T)
  a0 <- matrix(c(-theta0[2], 0, 0, theta0[2]), nrow = 2, byrow = T)

  u <- array(0, c(t,2))
  for (i in seq(t)) {
    u[i,] = y[i,] %*% b0 + x[i,] %*% a0
  }
  cat('\nResidual variance-covariance matrix (restricted)\n')
  omega0 <- t(u) %*% u/t
  print(omega0)

    
  cat('\nParameter estimates (unconstrained)   = ', theta1)
  cat('\nLog of the likelihood (unconstrained) = ', lf1)
  cat('\n\n')

  cat('\nParameter estimates (constrained)    = ', theta0)
  cat('\nLog of the likelihood (constrained)  = ', lf0)
  cat('\n\n')

  # Likelihood ratio test       
  lr <- -2*(t*lf0 - t*lf1)    
  cat('\nLikelihood ratio test                = ', lr)
  cat('\np-value                              = ', 1-pchisq(lr, 1) )

  # Alternative form of the LR test
  b1 <- matrix(c(1, -theta1[3], -theta1[1], 1), nrow = 2, byrow = T)
  b0 <- matrix(c(1, -theta0[3], -theta0[1], 1), nrow = 2, byrow = T)
    
  lr_alternative <- t*( log(det(omega0)) - log(det(omega1)) ) - 2*t*( log(abs(det(b0))) - log(abs(det(b1))) )
 
  cat('\nLikelihood ratio test (alternative)  = ', lr_alternative)
  cat('\np-value                              = ', 1-pchisq(lr_alternative,1) )  
}








