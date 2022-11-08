#==========================================================================
#
#    Program to find the MLEs using one iteration of Newton-Raphson
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
# (1) - Simulate the model 

x <- c(1, 2, 4, 5, 8)

beta <- 1.0
sig2 <- 4.0
t <- length(x)
y <- beta*x + sqrt(sig2)*rnorm(t)


# (2) - Evaluate gradient and Hessian at theta0 

beta_0  <- 1.0
sig2_0  <- 4.0
theta_0 <- beta_0 | sig2_0

g       <- cbind(rep(0, 2))
g[1,1]  <- sum( (y - beta_0*x) * x ) / sig2_0
g[2,1]  <- -0.5*t/sig2_0 + 0.5*sum( (y - beta_0*x)^2 )/sig2_0^2

h       <- matrix(0, c(2,2), nrow = 2, ncol = 2)
h[1,1]  <- -sum(x^2)/sig2_0
h[1,2]  <- -sum( (y - beta_0*x) * x )/sig2_0^2
h[2,1]  <- -sum( (y - beta_0*x) * x )/sig2_0^2
h[2,2]  <- 0.5*t/sig2_0^2 - sum( (y - beta_0*x)^2 )/sig2_0^3


# (3) - Average log-likelihood at theta0

a_0     <- -0.5*log(2*pi) - 0.5*log(sig2_0) - 0.5*mean( (y - beta_0*x)*2 )/sig2_0

# (4) - Newton-Raphson update

theta_1 <- theta_0 - inv(h) %*% g
beta_1  <- theta_1[1]
sig2_1  <- theta_1[2]

# (5) - Average log-likelihood at theta1

a_1     <- -0.5*log(2*pi) - 0.5*log(sig2_1) - 0.5*mean( (y - beta_1*x) * 2 )/sig2_1


cat('\nAverage log-likelihood at theta0 = ', a_0)
cat('\nAverage log-likelihood at theta1 = ', a_1, '\n')

