#=========================================================================
#
#     Program to do LR, Wald and LM test of Weibull Distribution
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(13, kind="Mersenne-Twister")

#
#--------------------------- Functions -----------------------------------
# 

#-------------------------------------------------------------------------
# Wrapper function to calculate inverse of a given matrix
#-------------------------------------------------------------------------
inv <- function (M) {
  return(solve(M))
}

#-------------------------------------------------------------------------
#  Unconstrained likelihood function 
#-------------------------------------------------------------------------
loglu <- function (b,y) {
  logl <- -mean( logltu(b,y) )
  return(logl)
}

#-------------------------------------------------------------------------
#  Unconstrained likelihood function at each observation
#-------------------------------------------------------------------------
logltu <- function(b,y) {
  alpha <- b[1]
  beta  <- b[2]
  lt <- log(alpha)+log(beta)+(beta-1)*log(y)- alpha*y^(beta) 
  return(lt)
}

#-------------------------------------------------------------------------
#  Constrained likelihood function
#-------------------------------------------------------------------------
loglc <- function(b,y) {
  logl <- -mean(logltc(b,y))
  return(logl)
}

#-------------------------------------------------------------------------
#  Constrained likelihood function at each observation
#-------------------------------------------------------------------------
logltc <- function(alpha,y) {
  beta <- 1
  lt <- log(alpha)+log(beta)+(beta-1)*log(y)- alpha*y^(beta)
  return(lt)
}

# load numeric lib required for this function
source("EMTSUtil.R")

#
#--------------------------- Exercise Functions -----------------------------------
#
test_weibull <- function() {
  # Simulate data
  t     <- 20
  alpha <- 1
  beta  <- 2
  
  # Generate Weibull random numbers
  y <- rweibull(t, scale = (1/alpha)^(1/beta), shape = beta)
  
   # or load data
#     y <- c(0.293, 0.589, 1.374, 0.954, 0.608, 1.199, 1.464, 
#         0.383, 1.743, 0.022, 0.719, 0.949, 1.888, 0.754, 
#         0.873, 0.515, 1.049, 1.506, 1.090, 1.644)
  
  bstart <- c(alpha, beta)
  
  # Unconstrained estimation
  results <- optim(bstart, loglu, y =y, method="BFGS")
  bu <- results$par
  logL1 <- results$value
  
  logL1 <- -logL1
  
  cat('\n')
  cat('\nUnconstrained estimation results')
  cat('\nalpha  = ', bu[1])
  cat('\nbeta   = ' , bu[2] )
  cat('\nlog L  = ' , t*logL1)
  
  # Constrained estimation
  results <- optim(1, loglc, y = y, method="BFGS")
  b0 <- results$par
  logL0 <- results$value
  logL0 <- -logL0
  
  cat('\n ')
  cat('\nConstrained estimation results')
  cat('\nalpha  = ', b0)
  cat('\nbeta   = ' , 1)
  cat('\nlog L  = ', t*logL0)
  
  # Likelihood ratio test
  lr <- -2*t*(logL0 - logL1)
  pv  <- 1 - pchisq(lr,1)

  cat('\nLikelihood ratio test')
  cat('\nLR stat      = ', lr)
  cat('\np-value      = ', pv)
  cat('\n')
  
  # Wald test  
  h <- numhess(loglu, bu, y)
    
  cat('\nHessian evaluated at unconstrained estimates\n')
  print(h)
  cat('\n')
  
  r <- c(0, 1)
  q <- 1
  w <- t*cbind((r%*%bu - q)) %*% inv(r %*% h %*% cbind(r)) %*% (r %*% bu - q)
  w <- as.numeric(w)
  pv <- 1 - pchisq(w,1)
    
  cat('\nWald test')
  cat('\nWald stat    = ', w )
  cat('\np-value      = ', pv )
  cat('\n' )

  #LM test
  th0  <- c(b0, 1)
  gmat <- numgrad(logltu, th0, y)  
  g    <- cbind((colMeans(gmat)))
  j    <- t(gmat) %*% gmat/t
  lm   <- t * t(g) %*% inv(j) %*% g
  pv   <- 1 - pchisq(lm, 1)

  cat('\nLM test')
  cat('\nLM stat    = ', lm )
  cat('\np-value    = ', pv)  
}




