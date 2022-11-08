#=========================================================================
#
#    Program to estimate model by full information maximum likelihood.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
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

#-----------------------------------------------------------------
# Log-likelihood function
#-----------------------------------------------------------------
        
neglog <- function(theta,y,x) {
   lf <- -mean(lnlt(theta,y,x))
   return(lf)
}

#-----------------------------------------------------------------
# Log-likelihood function at each observation 
#-----------------------------------------------------------------
lnlt <- function (theta,y,x) {
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

#-----------------------------------------------------------------
# Instrumental variable estimation
#-----------------------------------------------------------------
iv <- function(y,w,x) {
  
  tmp0   <- inv(t(w) %*% x %*% inv(t(x) %*% x) %*% t(x) %*% w)
  tmp1   <- (t(w) %*% x %*% inv( t(x) %*% x) %*% t(x) %*% y)
  b      <- tmp0 %*% tmp1                                             # IV estimates                            
  
  # Standard error of regression
  e       <- y - w %*% b  
  t       <- nrow(y)
  sigma   <- sqrt(t(e) %*% e/t)
  
  # Variance-covariance matrix
  vcov    <- c(sigma^2) *tmp0
  sterr   <- sqrt(diag(vcov))
  tstats  <- b/sterr  
  return(b)  
}

#
#--------------------------- FIML Model -----------------------------------
#

linear_fiml <- function() {
  # Simulate
  t <- 500
  make.new.random <- FALSE
  if (make.new.random) {
    simResults <- simulatedata(t)
    x <- simResults$x
    y <- simResults$y    
  } else {
    simResults <- read.table("linear_fimldata.dat")
    y <- as.matrix(simResults[,1:2])
    x <- as.matrix(simResults[,3:4])
  }
  # Estimate the model
  theta_0 <- runif(4)
  results <- optim(theta_0, neglog, y = y, x = x, method="BFGS")
  theta <- results$par

  ht      <- numhess(neglog,theta,y,x)
  cov     <- (1/t)*inv(ht)  
  u       <- array(0, c(t,2))
  
  b0 <- matrix(c(1, -theta[3], -theta[1], 1), nrow = 2, byrow = T)
  a0 <- matrix(c(-theta[2], 0, 0, -theta[4]), nrow = 2, byrow = T)
  for (i in seq(t)) {
    u[i,] <- y[i,] %*% b0 + x[i,] %*% a0    
  }
  
  cat('\n\nBeta 1 and se = ', theta[1],'   ', sqrt(cov[1,1]))
  cat('\nBeta 2 and se = ', theta[2],'   ', sqrt(cov[2,2]))
  cat('\nBeta 3 and se = ', theta[3],'   ', sqrt(cov[3,3]))
  cat('\nBeta 4 and se = ', theta[4],'   ', sqrt(cov[4,4]))

  cat('\nResidual variance-covariance matrix\n')
  print(t(u) %*% u/t)
  
  #Instrumental variable estimation
  tmpy1 <- as.matrix(y[,1])
  tmpy2x1 <- matrix(c(y[,2], x[,1]), ncol=2)
  tmpx1x2 <- x[,1:2]
  beta1_iv <- iv( tmpy1, tmpy2x1, tmpx1x2)

  cat('\nIV estimate of beta1   = ', beta1_iv[1])
  cat('\nIV estimate of alpha1  = ', beta1_iv[2])
  
  tmpy2 <- as.matrix(y[,2])
  tmpy1x2 <- matrix(c(y[,1], x[,2]), ncol=2)
  beta2_iv <- iv( tmpy2, tmpy1x2, tmpx1x2)  

  cat('\nIV estimate of beta2   = ', beta2_iv[1])
  cat('\nIV estimate of alpha2  = ', beta2_iv[2])
}



