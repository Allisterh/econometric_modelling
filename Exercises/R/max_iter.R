#==========================================================================
#
#   Program to find the MLEs using the Newton-Raphson algorithm
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(66, kind="Mersenne-Twister")

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
#--------------------------- Algorithms -----------------------------------
#

# Simulate the model 
x <- c(1, 2, 4, 5, 8)

beta <- 1.0
sig2 <- 4.0
t <- length(x)
y <- beta*x + sqrt(sig2)*rnorm(t)


# Evaluate gradient and Hessian at theta0
beta  <- 1.0
sig2  <- 4.0
theta <- beta | sig2

g       <- cbind(rep(0, 2))
g[1,1]  <- sum( (y - beta*x) * x ) / sig2
g[2,1]  <- -0.5*t/sig2 + 0.5*sum( (y - beta*x)^2 )/sig2^2

h       <- matrix(0, c(2,2), nrow = 2, ncol = 2)
h[1,1]  <- -sum(x^2)/sig2
h[1,2]  <- -sum( (y - beta*x) * x )/sig2^2
h[2,1]  <- -sum( (y - beta*x) * x )/sig2^2
h[2,2]  <- 0.5*t/sig2^2 - sum( (y - beta*x)^2 )/sig2^3

# Average log-likelihood at theta0
a_0     <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( (y - beta*x)*2 )/sig2

# Newton-Raphson
tol <- 0.001
while (t(g) %*% g > tol) 
{
  theta <- theta - inv(h) %*% g
  beta  <- theta[1]
  sig2  <- theta[2]
  
  g      <- cbind(rep(0, 2))
  g[1,1]  <- sum( (y - beta*x)*x )/sig2
  g[2,1]  <- -0.5*t/sig2 + 0.5*sum( (y - beta*x)^2 )/sig2^2
  
  h       <- matrix(0, c(2,2), nrow = 2, ncol = 2)
  h[1,1]  <- -sum(x^2)/sig2
  h[1,2]  <- -sum( (y - beta*x)*x )/sig2^2
  h[2,1]  <- -sum( (y - beta*x)*x )/sig2^2
  h[2,2]  <- 0.5*t/sig2^2 - sum( (y - beta*x)^2 )/sig2^3 
}


# Optimum Average log-likelihood at final theta

a_1     <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( (y - beta*x)*2 )/sig2

cat('\nAverage log-likelihood at theta0 ', a_0, '\n')
cat('\nAfter applying newton raphson')
cat('\nFinal value of beta              = ', beta)
cat('\nFinal value of variance          = ', sig2)
cat('\nFinal average log-likelihood     = ', a_1, '\n')



