#=========================================================================
#
#    Program to perform a likelihood ratio test of a regression model
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
# 

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
# Constrained log-likelihood function
#--------------------------------------------------------------------------------
neglog0 <- function(theta,y,x1,x2){
  lf <- -mean(lnl0(theta,y,x1,x2))
  return(lf)
}

#--------------------------------------------------------------------------------
# Constrained log-likelihood function at each observation
#-------------------------------------------------------------------------------- 
lnl0 <- function(theta,y,x1,x2) {
  m <-  theta[1] + theta[2]*x1 + (1-theta[2])*x2   # Mean                        
  s2 <- theta[3]                                   # Variance                    
  z <-  (y - m)/sqrt(s2)                           # Standardised residual       
    
  rval <- -0.5*log(2*pi) - 0.5*log(s2) - 0.5*z^2
  return(rval)
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
#--------------------------- LR Test Functions -----------------------------------
#

linear_lr <- function() {
  
  # Simulate Data
  t <- 200
  make.new.random <- FALSE
  if (make.new.random) {
    simResults <- simulatedata(t)
    x <- simResults$x
    y <- simResults$y    
  }
  else {
    simResults <- read.table("linear_testdata.dat")
    y <- as.matrix(simResults[,1])
    x <- as.matrix(simResults[,2:4])
  }   
  x1 <- x[,2]
  x2 <- x[,3]
    
  # Estimate the unconstrained model  
  theta_0 <- runif(4) 
  estResults <- optim(theta_0, neglog1, y = y, x1 = x1, x2 = x2, method="BFGS") 
  
  theta1 <- estResults$par
  lf1 <- estResults$value    
  lf1 <- -lf1

  # Estimate the constrained model  
  theta_0 <- runif(3)
  estResults <- optim(theta_0, neglog0, y = y, x1 = x1, x2 = x2, method="BFGS")
  
  theta0 <- estResults$par
  lf0 <- estResults$value  
  lf0 <- - lf0

  cat('\nParameter estimates (unconstrained)   = ', theta1)
  cat('\nLog of the likelihood (unconstrained) = ', lf1)
  cat('\n')

    
  cat('\nParameter estimates (constrained)     = ', theta0)
  cat('\nLog of the likelihood (constrained)   = ', lf0)
  cat('\n')

  # Likelihood ratio test       
  lr <- -2*(t*lf0 - t*lf1)
  pv <- 1 - pchisq(lr, 1)
    
  cat('\nLikelihood ratio test                 = ', lr)
  cat('\np-value                               = ', pv )
  cat('\n') 
}


