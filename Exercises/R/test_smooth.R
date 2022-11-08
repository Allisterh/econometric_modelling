#=========================================================================
# 
# Program for Neyman's Smooth Goodness of Fit Test   
#   
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

set.seed(123457, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions ----------------------------
# 
# load required functions - inv
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Program used to compute normalizing constant numerically
#-------------------------------------------------------------------------
normconst <- function(u,theta) {
  z <- u - 0.5
  t <- theta[1]*sqrt(3)*2*z + theta[2]*sqrt(5)*(6*z^2 - 0.5) + theta[3]*sqrt(7)*(20*z^3 - 3*z) + theta[4]*3*(70*z^4 - 15*z^2 + 3/8)
  y <- exp(1 + t)
  return (y)
}


#-------------------------------------------------------------------------
# Log-likelihood function
#-------------------------------------------------------------------------
neglog <- function(b,u) {
  z <- u - 0.5
  c  <- integrate(function(u) normconst(u,b),0,1)$value
  f  <- 1 + b[1]*sqrt(3)*2*z + b[2]*sqrt(5)*(6*z^2 - 0.5) + b[3]*sqrt(7)*(20*z^3 - 3*z) + b[4]*3*(70*z^4 - 15*z^2 + 3/8) - log(c)
  lf <- -mean(f)
  return(lf)  
}

#
#--------------------------- Tests --------------------------------------
#

test_smooth <- function() {
  # Simulate data
  t <- 1000
  y <- rnorm(t)        # normal (simulate under H0)

  # Transform the data to uniform values and estimate my maximum likelihood
  u <- pnorm(y)

  theta0 <- 0.1*rep(1, 4)
 
  # Estimate model
  estResults <- optim(theta0, neglog, u=u, method="BFGS", hessian=T)
  theta <- estResults$par
  lnlu  <- estResults$value
  H     <- estResults$hessian

  # LR test
  lnl1 <- -lnlu         
  lnl0 <- 0.0     # Under the null lnl is zero as f under H0 is 1 so lnl=0
  lr   <- -2*t*(lnl0 - lnl1)
  dof <- 4
  cat('\nLR test statistic    = ', lr)
  cat('\nDegrees of freedom   = ', dof)
  cat('\nP-value              = ', 1-pchisq(lr,dof))
  
  cat('\n')
  # Wald test
  wd <- t*theta %*% (inv(H)) %*% theta
  cat('\nWald statistic       = ', wd)
  cat('\nP-value              = ', 1-pchisq(wd,dof))
  cat('\n')

  # LM test
  z <- u - 0.5
  phi1 <- sqrt(3)*2*z
  phi2 <- sqrt(5)*(6*z^2 - 0.5)
  phi3 <- sqrt(7)*(20*z^3 - 3*z)
  phi4 <- 3*(70*z^4 - 15*z^2 + 3/8)

  lm <- (sum(phi1)/sqrt(t))^2 + (sum(phi2)/sqrt(t))^2 + (sum(phi3)/sqrt(t))^2 + (sum(phi4)/sqrt(t))^2
 
  cat('\nLM statistic         = ', lm)
  cat('\nP-value              = ', 1-pchisq(lm,dof))  
}
 

