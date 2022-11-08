#=========================================================================
#
#    Program to perform a Wald test of a regression model
#
#=========================================================================


rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
source("EMTSUtil.R")

#--------------------------------------------------------------------------------
# Simulate data an return the results as a list of x and y
#--------------------------------------------------------------------------------
simulatedata <- function(t)  {
  beta0  <- 1.0                               # Parameter values                            
  beta1  <- 0.7
  beta2  <- 0.3
  sig   <- sqrt(4.0)                              

  x1 <- rnorm(t)                                                                 
  x2 <- rnorm(t)
  u  <- sig*rnorm(t)                        # Generate disturbances          
  y  <- beta0 + beta1*x1 + beta2*x2 + u
    
  ones <- rep(1, t)    
  x  <- matrix(c(ones, x1, x2), ncol = 3)
  return(list(x=x, y=y))  
}  

#--------------------------------------------------------------------------------
# Unconstrained log-likelihood function
#--------------------------------------------------------------------------------
neglog1 <- function(theta,y,x1,x2) {
  lf <- -mean(lnl1(theta,y,x1,x2))
  return(lf)
}

#--------------------------------------------------------------------------------
# Unconstrained log-likelihood function at each observation
#--------------------------------------------------------------------------------     
lnl1 <- function (theta,y,x1,x2) {
  m    <- theta[1] + theta[2]*x1 + theta[3]*x2       # Mean                        
  s2   <- theta[4]                                   # Variance                    
  z    <- (y - m)/sqrt(s2)                           # Standardised residual       
    
  rval <- -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z^2  
  return(rval)
}

#
#--------------------------- Wald Test Functions -----------------------------------
#

linear_w <- function() {
  t <- 200
  make.new.random <- FALSE
  if (make.new.random) {
    simResults <- simulatedata(t)
    x <- simResults$x
    y <- simResults$y    
  } else {
    simResults <- read.table("linear_testdata.dat")
    y <- as.matrix(simResults[,1])
    x <- as.matrix(simResults[,2:4])
  }
  x1 <- x[,2]
  x2 <- x[,3]
 

  # Estimate the unconstrained model  
  theta_0 <- runif(4) 
  estResults <- optim(theta_0, neglog1, y = y, x1 = x1, x2 = x2, method="BFGS") 
  
  theta <- estResults$par
  ht    <- numhess(neglog1, theta, y, x1, x2)               # As f = negative log likelihood
  vcov1 <- (1/t)*inv(ht)
 
  # Wald test       
  r <- rbind(c(0, 1, 1, 0))
  q <- 1
  w <- t ((r %*% theta - q)) %*% inv( r %*% vcov1 %*% t(r)) %*% (r %*% theta - q)
  pv <- 1 - pchisq(w,1)
  
 cat('\nvcov based on numerical derivatives\n')
 print(vcov1) 
 cat('\nWald test = ', w)
 cat('\np-value   = ', pv)  
} 
 
    

