#==========================================================================
#
#      Program to estimate a model by full information maximum likelihood.
#      This program demonstrates the equivalence of OLS and MLE in this
#      special case of the SUR model.
#
#      The set of equations is defined as yt*b + xt*a = u
#      where yt is a (1xn) set of dependent variables at time t
#            xt is a (1xk) set of explanatory variables at time t
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#

# load utility functions - inv and numhess functions
source("EMTSUtil.R")

#--------------------------------------------------------------------------------
# Generate simulated data by simulating the reduced form.
#   The set of equations is given by the bivariate system
#         y1t = beta1*y2t + alpha1*x1t + u1t
#         y2t = beta2*y1t + alpha2*x2t + u2t
#         where E(ut'*ut) = omega
#--------------------------------------------------------------------------------

simulatedata <- function(t, beta1, alpha1, beta2, alpha2, beta3, alpha3) { 
  omega <- matrix(c(2, 0, 0,
                    0, 1, 0,
                    0, 0, 5), nrow = 3, byrow = T)


  # Construct population parameter matrices 
  b <- matrix(c(1, -beta1, -beta2, 
                0, 1,   -beta3,
                0, 0,      1), nrow = 3, byrow = T)
   

  a <- matrix(c(-alpha1, -alpha2, -alpha3,
           0,       0,       0,
           0,       0,       0),  nrow = 3, byrow = T)
  
  
  # Construct exogenous variables 
  x <- cbind( 1*rnorm(t), 2*rnorm(t), 3*rnorm(t) )

  u <- matrix(rnorm(t*3), ncol = 3) %*% chol(omega)                   
  invB <- inv(b)
  
  # Simulate the model
  y <- array(0, c(t,3))
  
  for (i in seq(t)) {
    y[i,] <- -x[i,] %*% a %*% invB + u[i,] %*% invB
  }
  return(list(x=x, y=y, u=u, omega=omega)) 
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
lnlt <- function(theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)
  
  b <- matrix(c(1, -theta[2], -theta[4], 
              0,  1,  -theta[5],
              0,  0,    1), nrow = 3, byrow = T)
  
  a <- matrix (c(-theta[1], -theta[3], -theta[6],
                 0,        0,         0,
                 0,        0,         0), nrow = 3, byrow = T)
  
  u <- array(0, c(t,n))
  for (i in seq(t)) {
    u[i,] <- y[i,] %*% b + x[i,] %*% a       # Construct residuals 
  }
  
  omega <- t(u) %*% u/t                         # Concentrate out resid var-covar matrix  
  lnl <- array (0, c(t,1))  
  for (i in seq(t)) {
    lnl[i] <- -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u[i,] %*% inv(omega) %*% cbind(u[i,])
  }
  
  return (lnl)
}

#-----------------------------------------------------------------
# Create a results table
#-----------------------------------------------------------------
create.empty.table <- function() {
  rnames <- c("alpha1", "beta1", "alpha2", "beta2", "beta3", "alpha3")
  cnames <- c("Actual", "MLE", "OLS")
  results <- matrix(nrow = length(rnames), ncol = length(cnames), 
                    dimnames = list(rnames, cnames))
  return(results)
}


#------------------------------------------ Recursive method ------------------------------------------	

linear_recursive <- function() {
  t <- 200  
  beta1 <- 0.6 
  alpha1 <-  0.4
  beta2 <- 0.2
  alpha2 <- -0.5 
  beta3 <- 1.0 
  alpha3 <-  0.2
  
  simResults <- simulatedata(t, beta1, alpha1, beta2, alpha2, beta3, alpha3)
  x <- simResults$x
  y <- simResults$y
  u <- simResults$u
  omega <- simResults$omega
  
  theta_0 <- runif(6)
  
  # Estimate the model using BFGS algorithm and compute QMLE se
  estResults <- optim(theta_0, neglog, y = y, x = x, method="BFGS")
  theta <- estResults$par
   
  comparison <- create.empty.table()
  # Actual
  comparison[,1] <- c(alpha1, beta1, alpha2, beta2, beta3, alpha3)
  
  # MLE Estimates
  comparison[,2] <- c(theta)
  
  # OLS Estimates
  eq1 <- lm(y[,1] ~ x[,1] - 1)$coef
  eq2 <- lm(y[,2] ~ y[,1] + x[,1] - 1)$coef
  eq3 <- lm(y[,3] ~ y[,1] + y[,2] + x[,1] - 1)$coef
  comparison[,"OLS"] <- c(eq1, eq2, eq3)
  
  cat('\nComparison of true and estimated parameters\n') 
  print(comparison)
  cat('\n\n')
  
  omegahat <- t(u) %*% u/t              # Note that u now represents the residuals from estimating the model    
  
  comparison <- cbind(diag(omega), diag(omegahat))
  
  cat('\nComparison of true and estimated omega\n')
  print(comparison)
}

