#============================================================================
#
#   Program to demonstrate multiple roots of a regression model
#   with nonlinear parameterisation
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

set.seed(12, kind="Mersenne-Twister")


# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - figure, seqa
source("EMTSUtil.R")

# 
#------------------------- Nonlinear Regression -----------------------------
#
t <- 100                            

# Parameters
beta0 <- 2.5
beta1 <- -1.5
sig2  <- 1.0

# Generate data
u <- rnorm(t)
x <- matrix(runif(t*2), t, 2)
y <- beta0 + beta1*x[,1] + beta1^2*x[,2] + sqrt(sig2)*u


ngrid <- 299
lnl   <- rep(0, ngrid)
g     <- array(0, c(ngrid,2))
b1    <- seqa(-2.1,0.01,ngrid)         


for (i in seq(ngrid)) {
  b0 <- beta0
  z      <- (y - b0 - b1[i]*x[,1] - b1[i]^2*x[,2])/sig2
  lnl[i] <- - log(2*pi) - 0.5*log(sig2) - 0.5*mean(z^2)
  g[i,1] <- mean( (y - b0 - b1[i]*x[,1] - b1[i]^2*x[,2]) )
  g[i,2] <- mean( (y - b0 - b1[i]*x[,1] - b1[i]^2*x[,2])*(x[,1] + 2*b1[i]*x[,2]) )
}


figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))
matplot(b1,cbind(g[,2], rep(0, ngrid)), type='l',
     main = '(a) Gradient',
     ylab = expression(G[T](theta)),
     xlab = expression(theta),
     bty = 'l')

plot(b1,lnl, type='l',
     main = '(b) Log-likelihood',
     ylab = expression(L[T](theta)),
     xlab = expression(theta),
     bty = 'l')

