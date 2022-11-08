#=========================================================================
#
#     Program to find the MLEs using grid search methods
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(123, kind="Mersenne-Twister")


# Simulate the model

x <- c(1, 2, 4, 5, 8)

beta <- 1.0
sig2 <- 4.0

t <- length(x)
y <- beta*x + sqrt(sig2)*rnorm(t)


# Grid search on gradient sig2 = 4      

sig2 <- 4.0
beta <- seq(0.5, 1.5, 0.1)
g1   <- rep(0, length(beta))

for (i in seq(beta)) {
    g1[i] <- sum( (y - beta[i]*x) * x) / sig2
}

figure()
par(mfrow=c(1,2), xaxs = "i", yaxs="i")
plot(beta,g1, type="l",
     main = "Gradient: sig2 = 4.0",
     xlab = expression(beta),
     ylab = expression(paste("G"[T]*"","(", beta, ")")))
lines(beta, rep(0, length(beta)))

# Grid search on gradient sig2 = 3.5
sig2 <- 3.5
beta <- seq(0.5, 1.5, 0.1)
g2   <- rep(0, length(beta))

for (i in seq(beta)) {
    g2[i] <- sum( (y - beta[i]*x) * x) / sig2
}

plot(beta,g2, type="l",
     main = "Gradient: sig2 = 3.5",
     xlab = expression(beta),
     ylab = expression(paste("G"[T]*"","(", beta, ")")))
lines(beta, rep(0, length(beta)))


