#=========================================================================
#
#     Program to find the MLEs using graphical methods
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(123, kind="Mersenne-Twister")


# Simulate the model   

x    <- c(1, 2, 4, 5, 8)
beta <- 1.0
sig2 <- 4.0
t    <- length( x )
y    <- beta*x + sqrt(sig2)*rnorm(t)


# Plot average log-likelihood sig2 = 4
sig2 <- 4.0
beta <- seq(0.0, 2.0, 0.1)
a1   <- rep(0, length(beta))

for (i in seq(beta)) {
    z <- ( y - beta[i]*x )/sqrt(sig2)
    a1[i] <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( z^2 )
}

figure()
par(mfrow=c(2,2), xaxs = "i", yaxs="i")
plot(beta,a1, type="l",
     main = "Average log-likelihood: sig2 = 4.0",
     xlab = expression(beta),
     ylab = expression(paste("ln L"[T]*"","(", beta, ")")),
     ylim = c(-5.5, -1.5))


# Plot average log-likelihood sig2 = 3.5

sig2 <- 3.5
beta <- seq(0.0, 2.0, 0.1)
a2   <- rep(0, length(beta))

for (i in seq(beta)) {
    z <- ( y - beta[i]*x )/sqrt(sig2)
    a2[i] <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( z^2 )
}

plot(beta,a2, type="l",
     main = "Average log-likelihood: sig2 = 3.5",
     xlab = expression(beta),
     ylab = expression(paste("ln L"[T]*"","(", beta, ")")),
     ylim = c(-5.5, -1.5))



# Plot average log-likelihood beta = 1.0
beta <- 1.0
sig2 <- seq(1, 11, 0.5)
a3   <- rep(0, length(sig2))

for (i in seq(sig2)) {
    z <- ( y - beta*x)/sqrt(sig2[i])
    a3[i] <- -0.5*log(2*pi) - 0.5*log(sig2[i]) - 0.5*mean( z^2 )
}

plot(sig2, a3, type="l",
     main = "Average log-likelihood: beta = 1.0",
     xlab = expression(sigma^2),
     ylab = expression(paste("ln L"[T]*"","(", sigma^2, ")")))



# Plot average log-likelihood beta = 0.9
beta <- 0.9
sig2 <- seq(1, 11, 0.5)
a4   <- rep(0, length(sig2))

for (i in seq(sig2)) {
    z <- ( y - beta*x)/sqrt(sig2[i])
    a4[i] <- -0.5*log(2*pi) - 0.5*log(sig2[i]) - 0.5*mean( z^2 )
}

plot(sig2, a4, type="l",
     main = "Average log-likelihood: beta = 0.9",
     xlab = expression(sigma^2),
     ylab = expression(paste("ln L"[T]*"","(", sigma^2, ")")))
