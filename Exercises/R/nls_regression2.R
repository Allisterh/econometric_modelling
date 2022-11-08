#============================================================================
#
#   Program to estimate a nonlinear regression model
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

set.seed(12, kind="Mersenne-Twister")

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Log-likelihood function
#----------------------------------------------------------------------------
neglog <- function(b,y,x) {
  t <- length(y)
  m <- b[1] + b[2]*x
  u <- y^b[3] - m
  
  # Concentrate the likelihood
  s2 <- t(u) %*% u/t             
  lt <- - 0.5*log(2*pi) - 0.5*log(s2) + log(abs(b[3])) + (b[3]-1)*log(y) - 0.5*(u^2)/s2
  lf <- -mean(as.real(lt))  
  return(lf)
}


# 
#--------------------------------- Regression Model ---------------------
#
nls_regression2 <- function( ) {
  t <- 100
  
  # Parameters
  beta0 <- 10
  beta1 <- 2
  beta2 <- 0.5
  sig2  <- 0.1
  
  # Generate the data                   
  x <- runif(t)^2
  u <- rnorm(t)
  y <- (beta0 + beta1*x + sqrt(sig2)*u)^(1/beta2)
  
  # Estimate the model using BGS and compute Hessian se   
  start <- c(beta0,  beta1, beta2)
  estResults <- optim(start, neglog, y=y, x=x, method="BFGS", hessian=T)
  bhat <- estResults$par
  hess <- estResults$hess
  
  vc <- (1/t)*inv(hess)
  
  cat('\n  Estimates  \n')
  print(cbind(bhat))
  
}

