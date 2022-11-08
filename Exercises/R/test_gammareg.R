#=========================================================================
#
# Program to do LR, Wald and LM tests of a gamma regression model
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

#-------------------------------------------------------------------------
# Define the unconstrained log-likelihood at each observation
#-------------------------------------------------------------------------
loglt_u <- function(theta,y,x) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  rho   <- theta[3]
  mu_t  <- beta0 + beta1*x
  
  LogLt <- -log(gamma(rho)) - rho*log(mu_t) + (rho-1)*log(y) -y/mu_t
  return(LogLt)  
}
  
#-------------------------------------------------------------------------
# -log likelihood summed over all observations
#-------------------------------------------------------------------------
logl_u <- function(theta,y,x) {
  LogLt <- loglt_u(theta,y,x)
  LogL  <- -sum(LogLt)
  return(LogL)
}

#-------------------------------------------------------------------------
# constrained log-likelihood at each observation
#-------------------------------------------------------------------------
loglt_c <- function(theta,y,x) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  rho   <- 1
  mu_t  <- beta0 + beta1*x
  
  LogLt <- -log(gamma(rho)) - rho*log(mu_t) + (rho-1)*log(y) - y/mu_t
  return(LogLt)  
}

#-------------------------------------------------------------------------
# -log likelihood summed over all observations
#-------------------------------------------------------------------------
logl_c <- function(theta,y,x) {
  LogLt <- loglt_c(theta,y,x)
  LogL  <- -sum(LogLt)
  return(LogL)
}

# load required functions
source("EMTSUtil.R")

#
#--------------------------- Gamma Model Tests -----------------------------------
#

test_gammareg <- function() {
  
  # Simulation  
  beta0 <- 1
  beta1 <- 2
  T     <- 2000
  
  x  <- runif(T)
  mu_t <- beta0 + beta1*x
  u_t  <- runif(T)
  y  <- -mu_t*log(1 - u_t)
  
  # Unconstrained estimation  
  th0 <- c(1,1,1)
  results <- optim(th0, function(theta) {logl_u(theta, y, x)})
  theta_1 <- results$par
  logL_1  <- results$value
  
  cat('\n')
  cat('\nUnconstrained estimation results')
  cat('\nbeta0_hat_1  = ',theta_1[1])
  cat('\nbeta1_hat_1  = ',theta_1[2])
  cat('\nrho_hat_1    =', theta_1[3])
  cat('\nlogL_1       = ',-logL_1)
  cat('\n') 

  # LR Test  
  th0 <- c(1,1,1)
  results <- optim(th0, function(theta) {logl_c(theta, y, x)})
  theta_0 <- results$par
  logL_0  <- results$value
  theta_0[3] <- 1
  
  cat('\n')
  cat('\nConstrained estimation results')
  cat('\nbeta0_hat_0  = ',theta_0[1])
  cat('\nbeta1_hat_0  = ',theta_0[2])
  cat('\nrho_hat_1    = ',theta_0[3])
  cat('\nlogL_0       = ',-logL_0)
  cat('\n')
  
  LR <- 2*(logL_0 - logL_1)
  p  <- 1 - pchisq(LR,1)
  
  cat('\nLikelihood ratio test')
  cat('\nLR stat      = ',LR)
  cat('\np-value      = ',p)
  cat('\n\n')
 
  
  # Wald test with Hessian computed from numerical derivatives  
  H_1 <- -numhess(logl_u,theta_1,y,x)
  
  R <- rbind(c(0,0,1))
  Q <- 1
  
  cat('\nHessian evaluated at theta_hat_1, determined numerically')
  cat('\nH(th_hat_1)\n')
  print(H_1)
  cat('\n')
  
  WH <- t( (R %*% theta_1 - Q) ) %*% inv( R %*% inv(-H_1) %*% t(R) ) %*% (R %*% theta_1 - Q)
  pH <- 1 - pchisq(WH,1)
  
  cat('\nWald tests with numerical derivatives')
  cat('\nUsing Hessian')
  cat('\nWald stat    = ',WH)
  cat('\np value      = ',pH)
  cat('\n')
  
  #LM test with OPG computed from numerical derivatives  
  Gt  <- numgrad(loglt_u,theta_0,y,x)
  J_0 <- t(Gt) %*% Gt
  
  cat('\nOPG evaluated at theta_hat_0, determined numerically')
  cat('\nJ(th_hat_0)\n')
  print(J_0)
  cat('\n')
  
  G_0 <- as.matrix(colSums(Gt))
  
  cat('\ngradient evaluated at theta_hat_0')
  cat('\nG(th_hat_0)\n')
  print(G_0)
  cat('\n')
  
  LMJ      <- t(G_0)  %*% inv(J_0) %*% G_0
  pJ       <- 1 - pchisq(LMJ,1)
  
  cat('\nLM test with OPG matrix computed from numerical derivatives')
  cat('\nLM stat      = ',LMJ)
  cat('\np value      = ',pJ)
  cat('\n')  
}







