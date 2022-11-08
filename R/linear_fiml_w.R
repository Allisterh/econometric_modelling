#=========================================================================
#
#  Program to estimate model by full information maximum likelihood 
#  and do a Wald test. 
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

#------------------------------------- FIML Wald Test  -------------------------------------#
 
linear_fiml_w <- function() {
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

  ht    <- numhess(neglog1,theta1,y,x)             
  vcov1 <- (1/t)*inv(ht)

  cat('\nParameter estimates (unconstrained)\n')
  print(theta1)
  cat('\nLog of the likelihood (unconstrained) = ', lf1)

  # Wald test
  q <- 0
  r <- rbind(c(0, 1, 0, 1))
  w <- t ((r %*% theta1 - q)) %*% inv(r %*% vcov1 %*% t(r)) %*% (r %*% theta1 - q)
  pv <- 1 - pchisq(w,1)
  
  cat('\nCovariance matrix\n')
  print(vcov1)

  cat('\nWald test     = ', w)
  cat('\np-value       = ', pv)  
}  	 






