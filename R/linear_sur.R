#=============================================================================
# Program to estimate a SUR model by full information maximum likelihood.
# The set of equations is defined as yt*b + xt*a = u, with b = i
#          where yt is a (1xn) set of dependent variables at time t
#           xt is a (1xk) set of explanatory variables at time t
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")
    
# load utility functions - inv 
source("EMTSUtil.R")

#
#--------------------------- Helper Functions -----------------------------------
#

#--------------------------------------------------------------------------------
# Simulate data an return the results as a list of x, y, alpha and omega
#-------------------------------------------------------------------------------- 
simulatedata <- function(t) {
  alpha1 <-  0.4
  alpha2 <- -0.5
  alpha3 <-  1.0
  omega <- matrix(c(1, 0.5, -0.1,
                  0.5, 	1.0, 	 0.2,
                  -0.1, 0.2, 	 1.0), nrow = 3, byrow = T)
  b <- diag(1, nrow(omega)) 
  a <- matrix(c(-alpha1, 	0,	0,
              0,  -alpha2,  0,
              0,  	0,	-alpha3), nrow = 3, byrow = T)
  # Exogenous variables                       
  x <- cbind(1*rnorm(t), 2*rnorm(t), 3*rnorm(t))
  
  
  # Construct disturbances                   
  u <- matrix(rnorm(t*3), ncol = 3) %*% chol(omega)                   
  invB <- inv(b)
  
  # Simulate the model by simulating the reduced form   
  # Note that the SUR system is the reduced form
  y <- array(0, c(t,3))
  
  for (i in seq(t)) {
    y[i,] <- -x[i,] %*% a %*% invB + u[i,] %*% invB
  }
  return(list(x=x, y=y, alpha=c(alpha1,alpha2,alpha3), omega=omega))
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
    
  b <- diag(1, n)
  a <- matrix(c( -theta[1], 0, 0,
                0, -theta[2], 0,
                0,    0, -theta[3]), 
                nrow = n, byrow=T)  
  
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

#
#--------------------------- SUR Model -----------------------------------
#
linear_sur <- function() {
  # simulate data
  t <- 500
  simResults <- simulatedata(t) 
  y <- simResults$y
  x <- simResults$x
  
  # Estimate the model  
  theta0 <- runif(3)
  estResults <- optim(theta0, neglog, 
                      y=y, x=x, method="BFGS")
  theta <- estResults$par
  
  # Compute OLS estimates for comparative purposes
  alpha1ols <-  lm(y[,1] ~ x[,1] - 1)$coef
  alpha2ols <-  lm(y[,2] ~ x[,2] - 1)$coef
  alpha3ols <- lm(y[,3] ~ x[,3] - 1)$coef
  
  alphaols <- c(alpha1ols, alpha2ols, alpha3ols)
    
  cat('\nComparing true and estimated parameter values')
  cat('\n        Actual     MLE       OLS\n')
  print( cbind(simResults$alpha, theta, alphaols) )

  # Compute residuals at optimal parameters
  n <- ncol(y)
  b <- diag(1, n)
  a <- matrix(c( -theta[1], 0, 0,
                0, -theta[2], 0,
                0,    0, -theta[3]), 
                nrow = n, byrow=T)
  u <- array(0, c(t,n))
  for (i in seq(t)) {
     u[i,] <- y[i,] %*% b + x[i,] %*% a    
  }
  omegahat <- t(u) %*% u/t

  cat('\nComparing true and estimated elements of Omega\n') 
  print( cbind(c(simResults$omega), c(omegahat)) )
}

 
   

