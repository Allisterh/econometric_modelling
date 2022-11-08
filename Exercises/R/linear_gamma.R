#============================================================================
#
#   Program to demonstrate the distribution of the 
#   estimates of a gamma regression model
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv
source("EMTSUtil.R")


# 
#--------------------------------- Gamma Regression Model -------------------
#
ndraws <- 5000                                                         

# Parameters
b0    <- c(1, 2)
rho   <- 0.25                           
alpha <- 0.1

# For the sample size of T = 10     
set.seed(1, kind="Mersenne-Twister")

t  <- 10                           
x  <- cbind(rep(1, t), rnorm(t))
z1 <- rep(0, ndraws)
z2 <- rep(0, ndraws)

pb <- txtProgressBar(min=0, max=ndraws, style=3)
for (i in seq(ndraws)) {
  u  <- alpha*rgamma(t, shape=rho,scale=1) - rho*alpha                      
  y  <- x %*% cbind(b0) + u
  b  <- lm(y ~ x - 1)$coef
  e  <- y - x %*% b
  
  s2 <- t(e) %*% e/t
  vc <- c(s2)*inv(t(x) %*% x)
  z1[i] <- (b[1] - b0[1])/sqrt(vc[1,1])
  z2[i] <- (b[2] - b0[2])/sqrt(vc[2,2])
  setTxtProgressBar(pb, i)
}       
close(pb)
# For the sample size of T = 100 
set.seed(1, kind="Mersenne-Twister")

t  <- 100                              
x  <- cbind(rep(1, t), rnorm(t))
z3 <- rep(0, ndraws)
z4 <- rep(0,ndraws)
pb <- txtProgressBar(min=0, max=ndraws, style=3)
for (i in seq(ndraws)) {
  
  u  <- alpha*rgamma(t, shape=rho,scale=1) - rho*alpha                      
  y  <- x %*% cbind(b0) + u
  b  <- lm(y ~ x - 1)$coef
  e  <- y - x %*% b
  s2 <- t(e) %*% e/t
  vc <- c(s2)*inv(t(x) %*% x)
  z3[i] <- (b[1] - b0[1])/sqrt(vc[1,1])
  z4[i] <- (b[2] - b0[2])/sqrt(vc[2,2])
  setTxtProgressBar(pb, i)
}
close(pb)

# Plot the results

figure()
par(mfrow=c(2,2))
hist(z1, col='blue', breaks=40,
     main = expression(paste("(a) T=10, ", beta[0])),
     xlab = '',
     ylab = '', xlim=c(-10,10))

hist(z2, col = 'blue', breaks=40,
     main = expression(paste("(a) T=10, ", beta[1])),
     xlab = '',
     ylab = '')

hist(z3, col = 'blue', breaks=40,
     main = expression(paste("(a) T=100, ", beta[0])),
     xlab = '',
     ylab = '')

hist(z4, col = 'blue', breaks=40,
     main = expression(paste("(a) T=100, ", beta[1])),
     xlab = '',
     ylab = '')
