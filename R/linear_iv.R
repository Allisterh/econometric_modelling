#==========================================================================
#  Program to estimate a model by full information maximum likelihood 
#  and instrumental variables.
#
#  The set of equations is defined as yt*b + xt*a = u
#  where yt is a (1xn) set of dependent variables at time t
#  xt is a (1xk) set of explanatory variables at time t
#==========================================================================

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
# The set of equations is given by the bivariate system
#         y1t = beta*y2t + u1t
#         y2t = gam*y1t + alpha*xt + u2t
#         where E(ut'*ut) = omega
#--------------------------------------------------------------------------------
 
simulatedata <- function(t) {  
  beta <- 0.6 
  gam <- 0.4 
  alpha <- -0.5 
  
  omega <- matrix(c(2.0, 0.0, 0.0, 1.0), nrow = 2, byrow=T)
  
  # Construct population parameter matrices     
  b <- matrix(c(1, -gam, -beta, 1.0), nrow = 2, byrow=T) 
  a <- matrix(c(0, -alpha), nrow = 1, byrow=T)
  
  # Construct exogenous variables from a uniform distribution                     
  x <- cbind(10*rnorm(t))  
  
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
        
neglog <- function(theta,y,x) {
   lf <- -mean(lnlt(theta,y,x))
   return(lf)
}

#-----------------------------------------------------------------
# Unconstrained Log-likelihood function at each observation 
#-----------------------------------------------------------------
lnlt <- function (theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)  
    
  b <- matrix(c(1, -theta[2], -theta[1], 1), nrow = 2, byrow=T)
  a <- matrix(c(0, -theta[3]), nrow = 1, byrow=T)
  u <- array(0, c(t,n))
  for (i in seq(t)) {
     u[i,] <- y[i,] %*% b + x[i,] %*% a    
  }
  i<-1
  me2 <- colMeans(u^2)
	omega <- diag(1,n)
	
	for (j in seq( ncol(omega) ) ) {
    omega[i,j] <- me2[i]
  	i<-i+1   
	}  
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
  b      <- tmp0 %*% tmp1                                           # IV estimates                          
  return(b)  
}

#
#--------------------------- FIML and IV Estimation -----------------------------------
#

linear_iv <- function() {  
  t <- 500
  theta_0 <- runif(3)
  
  simResults <- simulatedata(t)
  x <- simResults$x
  y <- simResults$y
  u <- simResults$u
  
  # Estimate the model using BFGS algorithm and compute QMLE se  
  estResults <- optim(theta_0, neglog, y = y, x = x, method = "BFGS")
  theta <- estResults$par
  
  
  cat('\nFIML (iterative) estimate of beta  = ', theta[1])
  cat('\nFIML (iteratuve) estimate of gam   = ', theta[2])
  cat('\nFIML (iterative) estimate of alpha = ', theta[3])
  cat('\n\n')
  
  t <- nrow(y)
  n <- ncol(y)
  b <- matrix(c(1, -theta[2], -theta[1], 1), nrow = 2, byrow=T)
  a <- matrix(c(0, -theta[3]), nrow = 1, byrow=T)
    
  u <- array(0, c(t,n))
  for (i in seq(t)) {
     u[i,] <- y[i,] %*% b + x[i,] %*% a    
  }
  # Concentrate out resid var-covar matrix and restrict it to be diagonal
  i<-1
  me2 <- colMeans(u^2)
  drv <- diag(1,n)
  
  for (j in seq(length(me2)) ) {
    drv[i,j] <- me2[i]
  	i<-i+1   
  }      
  cat('\nResidual variance-covariance matrix\n')
  print(drv)  

  # FIML analytical solution    
  beta_fiml <- sum(y[,1]*x) / sum(y[,2]*x)

  e1 <- y[,1] - beta_fiml * y[,2]

  num <- sum(y[,2]*e1)*sum(x^2) - sum(x*e1)*sum(y[,2]*x)

  den <- sum(y[,1]*e1)*sum(x^2) - sum(x*e1)*sum(y[,1]*x)
  
  gam_fiml <- num/den
  
  num <- sum(y[,1]*e1)*sum(y[,2]*x) - sum(y[,1]*x)*sum(y[,2]*e1)
  
  alpha_fiml <- num/den
  
  cat('\nFIML (analytical) estimate of beta   = ', beta_fiml)
  cat('\nFIML (analytical) estimate of gam    = ', gam_fiml)
  cat('\nFIML (analytical) estimate of alpha  = ', alpha_fiml)
  cat('\n\n')  
  
  # Instrumental variable estimation
  beta1_iv <- iv(y[,1], y[,2], x)                 		# First equation
  
  tmpy2 <- as.matrix(y[,2])
  tmpy1x <- matrix(c(y[,1], x), ncol=2)
  tmpe1x <- matrix(c(e1, x), ncol=2)
  
  beta2_iv <- iv(tmpy2, tmpy1x, tmpe1x)            	# Second equation 
  
  cat('\nIV estimate of beta  = ', beta1_iv)
  cat('\nIV estimate of gam   = ', beta2_iv[1])
  cat('\nIV estimate of alpha = ', beta2_iv[2])
  cat('\n\n')  
}


