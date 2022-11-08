#=========================================================================
#
# Program to do LR, Wald and LM tests of a exponential regression model
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(1234, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions ----------------------------
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
loglt_u <-function(theta,y,x) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu_t  <- beta0 + beta1*x

  LogLt <- -log(mu_t) - y/mu_t
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
# analytic gradient at each observation
#-------------------------------------------------------------------------
GtExp <-function (theta,y,x){
  T <- length(y)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu_t  <- beta0 + beta1*x

  Gt      <- array(0, dim = c(T,2))  
  Gt[,1] <- -mu_t^(-1) + y*mu_t^(-2)
  Gt[,2] <- x*Gt[,1]
  return(Gt)
}

#-------------------------------------------------------------------------
# analytic gradient summed over all observations
#-------------------------------------------------------------------------
GExp <- function(theta,y,x){
  Gt <- GtExp(theta,y,x)  
  G  <- as.matrix(colSums(Gt))
  return(G)
}

#-------------------------------------------------------------------------
# analytic hessian
#-------------------------------------------------------------------------
HExp <- function(theta,y,x){
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu_t  <- beta0 + beta1*x

  H      <- array(0, dim = c(2,2) )
  H[1,1] <- sum(mu_t^(-2) - 2*y*mu_t^(-3))
  H[1,2] <- sum( (mu_t^(-2) - 2*y*mu_t^(-3))*x )
  H[2,1] <- H[1,2]
  H[2,2] <- sum( (mu_t^(-2) - 2*y*mu_t^(-3))*x^2 )
  return(H)  
}


#-------------------------------------------------------------------------
# analytic information matrix
#-------------------------------------------------------------------------
IExp <- function(theta,y,x){
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu_t  <- beta0 + beta1*x
  
  I      <- array(0, dim = c(2,2))
  I[1,1] <- sum(mu_t^(-2))
  I[1,2] <- sum( (mu_t^(-2))*x )
  I[2,1] <- I[1,2]
  I[2,2] <- sum( (mu_t^(-2))*x^2 )
  return(I)  
}

#-------------------------------------------------------------------------
# constrained log-likelihood at each observation
#-------------------------------------------------------------------------
loglt_c <-function(theta,y,x){
  beta0 <- theta[1]
  beta1 <- 0
  mu_t  <- beta0 + beta1*x
  
  LogLt <- -log(mu_t) - y/mu_t
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

# load required functions - numhess, numgrad
source("EMTSUtil.R")


#
#--------------------------- Exp Model Tests ----------------------------
#

test_expreg <- function() {

  # Simulate
  beta0 <- 1
  beta1 <- 2
  T     <- 2000
  
  x  <- runif(T)
  mu_t <- beta0 + beta1*x
  u_t  <- runif(T)
  y  <- -mu_t*log(1 - u_t)
  
  # Unconstrained estimation  
  th0 <- c(1,1)
  results <- optim(th0, function(theta) {logl_u(theta, y, x)})
  theta_1 <- results$par
  logL_1  <- results$value
  
  cat('\n')
  cat('\nUnconstrained estimation results')
  cat('\nbeta0_hat_1  = ',theta_1[1])
  cat('\nbeta1_hat_1  = ',theta_1[2])
  cat('\nlogL_1       = ',-logL_1)
  cat('\n')
  
  #   Constrained estimation  
  th0 <- c(1,1)
  results <- optim(th0, function(theta) {logl_c(theta, y, x)})
  theta_0 <- results$par
  logL_0  <- results$value
  theta_0[2] <- 0
  
  cat('\n')
  cat('\nConstrained estimation results')
  cat('\nbeta0_hat_0  = ',theta_0[1])
  cat('\nbeta1_hat_0  = ',theta_0[2])
  cat('\nlogL_0       = ',-logL_0)
  cat('\n')
  
  
  # Perform Tests
  
  # LR Test  
  lr <- 2*(logL_0 - logL_1)
  p  <- 1 - pchisq(lr,1)
  
  cat('\nLikelihood ratio test')
  cat('\nLR stat      = ',lr)
  cat('\np-value      = ',p)
  cat('\n\n')
  
  # Wald test with analytic derivatives  
  H_1 <- HExp(theta_1,y,x)
  I_1 <- IExp(theta_1,y,x)
  Gt  <- GtExp(theta_1,y,x)
  J_1 <- t(Gt) %*% Gt
  
  cat('\nAnalytically determined matrices for covariance')
  cat('\nHessian evaluated at theta_hat_1')
  cat('\nH(th_hat_1)  = \n')
  print(H_1)
  cat('\nInformation matrix evaluated at theta_hat_1')
  cat('\nI(th_hat_1)  = \n')
  print(I_1)
  cat('\nOPG evaluated at theta_hat_1')
  cat('\nJ(th_hat_1)  = \n')
  print(J_1)
  cat('\n')
  
  R <- rbind(c(0, 1))
  Q <- 0
  
  WH <- t( (R %*% theta_1 - Q) ) %*% inv( R %*% inv(-H_1) %*% t(R) ) %*% (R %*% theta_1 - Q)
  WI <- t( (R %*% theta_1 - Q) ) %*% inv( R %*% inv(I_1) %*% t(R) ) %*% (R %*% theta_1 - Q)
  WJ <- t( (R %*% theta_1 - Q) ) %*% inv( R %*% inv(J_1) %*% t(R) ) %*% (R %*% theta_1 - Q)
  
  pH <- 1 - pchisq(WH,1)
  pI <- 1 - pchisq(WI,1)
  pJ <- 1 - pchisq(WJ,1)
  
  cat('\nWald tests with analytical derivatives')
  cat('\nUsing Hessian')
  cat('\nWald stat    = ',WH)
  cat('\np value      = ',pH)
  cat('\nUsing information matrix')
  cat('\nWald stat    = ',WI)
  cat('\np value      = ',pI)
  cat('\nUsing OPG')
  cat('\nWald stat    = ',WJ)
  cat('\np value      = ',pJ)
  cat('\n')
   
  #   Wald test with numerical derivatives  
  Gt  <- numgrad(loglt_u,theta_1,y,x)
  J_1 <- t(Gt) %*% Gt
  H_1 <- -numhess(logl_u,theta_1,y,x)
  
  
  cat('\nNumerically determined matrices for covariance')
  cat('\nHessian evaluated at theta_hat_1')
  cat('\nH(th_hat_1)  =\n')
  print(H_1)
  cat('\nOPG evaluated at theta_hat_1')
  cat('\nJ(th_hat_1)\n')
  print(J_1)
  cat('\n')
  
  WH <- t ( (R %*% theta_1 - Q) ) %*% inv( R %*% inv(-H_1) %*% t(R) ) %*% (R %*% theta_1 - Q)
  WJ <- t ( (R %*% theta_1 - Q) ) %*% inv( R %*% inv(J_1) %*% t(R) ) %*% (R %*% theta_1 - Q)
  
  pH <- 1 - pchisq(WH,1)
  pJ <- 1 - pchisq(WJ,1)
  
  cat('\nWald tests with numerical derivatives')
  cat('\nUsing Hessian')
  cat('\nWald stat    = ',WH)
  cat('\np value      = ',pH)
  cat('\nUsing OPG')
  cat('\nWald stat    = ',WJ)
  cat('\np value      = ',pJ)
  cat('\n')
  
  # LR tests with analytic derivatives  
  H_0 <- HExp(theta_0,y,x)
  I_0 <- IExp(theta_0,y,x)
  Gt  <- GtExp(theta_0,y,x)
  J_0 <- t(Gt) %*% Gt
  G_0 <- GExp(theta_0,y,x) 
  
  cat('\nAnalytically determined results')
  cat('\ngradient evaluated at theta_hat_0')
  cat('\nG(th_hat_0)\n')
  print(G_0)
  cat('\n')
  
  cat('\nHessian evaluated at theta_hat_0')
  cat('\nH(th_hat_0)\n')
  print(H_0)
  cat('\nInformation matrix evaluated at theta_hat_0')
  cat('\nI(th_hat_0)\n')
  print(I_0)
  cat('\nOPG evaluated at theta_hat_0')
  cat('\nJ(th_hat_0)\n')
  print(J_0)
  cat('\n')  
  
  LMH      <- t(G_0) %*% inv(-H_0) %*% G_0
  LMI      <- t(G_0) %*% inv(I_0)%*% G_0
  LMJ      <- t(G_0) %*% inv(J_0)%*% G_0
  
  pH       <- 1 - pchisq(LMH,1)
  pI       <- 1 - pchisq(LMI,1)
  pJ       <- 1 - pchisq(LMJ,1)
  
  cat('\nLM tests with analytical derivatives')
  cat('\nUsing Hessian')
  cat('\nLM stat      = ',LMH)
  cat('\np value      = ',pH)
  cat('\nUsing information matrix')
  cat('\nLM stat      = ',LMI)
  cat('\np value      = ',pI)
  cat('\nUsing OPG')
  cat('\nLM stat      = ',LMJ)
  cat('\np value      = ',pJ)
  cat('\n')
  
  # LR tests with numerical derivatives  
  Gt  <- numgrad(loglt_u,theta_0,y,x)
  J_0 <- t(Gt) %*% Gt
  H_0 <- -numhess(logl_u,theta_0,y,x)
  G_0 <- as.matrix(colSums(Gt))
  
  cat('\nNumerically determined results')
  cat('\ngradient evaluated at theta_hat_0')
  cat('\nG(th_hat_0)\n')
  print(G_0)
  cat('\n')
  cat('\nHessian evaluated at theta_hat_0 ')
  cat('\nH(th_hat_0)\n')
  print(H_0)
  cat('\nOPG evaluated at theta_hat_0')
  cat('\nJ(th_hat_0)\n')
  print(J_0)
  cat('\n')
  
  LMH      <- t(G_0) %*% inv(-H_0) %*% G_0
  LMJ      <- t (G_0) %*% inv(J_0) %*% G_0
  
  pH       <- 1 - pchisq(LMH,1)
  pJ       <- 1 - pchisq(LMJ,1)
  
  cat('\nLM tests with numerical derivatives')
  cat('\nUsing Hessian\n')
  cat('\nLM stat      = ',LMH)
  cat('\np value      = ',pH)
  cat('\nUsing OPG')
  cat('\nLM stat      = ',LMJ)
  cat('\np value      = ',pJ)
  cat('\n')   
}



