#=========================================================================
#
#     Program to do LR, Wald and LM tests of a simple regression model
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(1234, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -----------------------------------
# 

#-------------------------------------------------------------------------
# Wrapper function to calculate inverse of a given matrix
#-------------------------------------------------------------------------
inv <- function (M) {
  return(solve(M))
}

#
#--------------------------- Regression Model Tests -----------------------------------
#
test_regress <- function () {
  x    <- c(1,2,4,5,8)
  t    <- length(x)
  beta <- 1
  sig  <- 2
  
  # Generate data
  u <- sig*rnorm(t)
  y <- beta*x + u
  
  # Constrained estimation
  bhat0 <- 0.0
  sig20 <- mean((y - x*bhat0)^2)
  
  logL0  <- - 1/2*log(2*pi) - 1/2*log(sig20) - 1/(2*sig20*t)*sum((y - x*bhat0)^2)
  
  cat('\nConstrained estimation results')
  cat('\nbeta_hat_0   = ',bhat0)
  cat('\nsigma2_hat_0 = ',sig20)
  cat('\nlogL_0       = ',logL0)
  cat('\n\n')
             
  # Unconstrained estimation
  bhat1 <- lm(y~x - 1)$coef
  err2  <- (y - x*bhat1)^2
  sig21 <- mean(err2)
   
  
  logL1 <- - 1/2*log(2*pi) - 1/2*log(sig21) - 1/(2*sig21*t)*sum((y - x*bhat1)^2)
  
  cat('\nUnconstrained estimation results \n')
  cat('\nbeta_hat_1   = ',bhat1)
  cat('\nsigma2_hat_1 = ',sig21)
  cat('\nlogL_1       = ',logL1)
  cat('\n\n')
  
  # Likelihood ratio test
  lr  <- -2*t*(logL0 - logL1)
  pv  <- 1-pchisq(lr,1)
  
  cat('\nLikelihood ratio test')
  cat('\nLR stat      =  ' , lr)
  cat('\np-value      =  ', pv )
  cat('\n')
  
  # Wald test
  r      <- c(1, 0)
  q      <-  0 
  theta1 <- cbind(c(bhat1, sig21))
  
  i <- matrix(c( (x %*% x)/(t*sig21),  0.0, 0.0, 1/(2*(sig21)^2) ), nrow=2)
  w      <- t*(r %*% theta1 - q) %*% inv(r %*% i %*% r) %*% (r %*% theta1 - q)
  pv     <- pchisq(w,1)
  
  cat('\nWald test')
  cat('\nUnconstrained estimates = ', bhat1, sig21 )
  cat('\nWald stat               =  ', w )
  cat('\np-value                 =  ', pv)
  cat('\n')
  
  # LM test analytic derivatives
  g <- c( x %*% (y-x*bhat0)/(t*sig20),
        -(1/(2*sig20))+(1/(2*(sig20)^2*t))*sum( (y - x*bhat0)^2 ) )
  i <- matrix( c( (x %*% x)/(t*sig20),  0.0, 0.0, 1/(2*(sig20)^2) ), nrow=2 )
      
  lm <- t*g %*% inv(i) %*% g  
  pv <- pchisq(lm,1)
  
  
  cat('\nLM test')
  cat('\nConstrained estimates   = ', 0, sig20 )
  cat('\nLM stat                 =  ',lm)
  cat('\np-value                 =  ', pv)
  cat('\n')  
}
