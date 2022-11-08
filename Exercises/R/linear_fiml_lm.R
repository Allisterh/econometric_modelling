#=========================================================================
#
#  Program to estimate model by full information maximum likelihood 
#  and do a LM test. 
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

#------------------------------------- FIML LM Test  -------------------------------------#
 

linear_fiml_lm <- function() {  
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
  
  # Estimate the model  
  start <- runif(3) 
  estResults <- optim(start, neglog0, y = y, x = x, method = "BFGS")
  theta <- estResults$par

  theta0 <- c(theta[1:3], -theta[2])
  gmat <- numgrad(lnlt1,theta0,y,x)                               
  g    <- cbind(colMeans(gmat))
  j    <- t(gmat) %*% gmat/t
  lm   <- t* t( g ) %*% inv(j) %*% g
  vcov <- (1/t) * inv(j)

  cat('\nGradient evaluated at contrained estimates\n')
  print(g)

	cat('\nOuter product of gradients matrix\n')
  print(j)
    
  cat('\nCovariance matrix\n')
  print(vcov)

  cat('\nLM test = ', lm)
  cat('\np-value = ', 1-pchisq(lm,1))  
}





