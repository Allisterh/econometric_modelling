#============================================================================
#
#    Program to reproduce the results in Chapman and Pearson, 
#    Journal of Finance (2000), p.355-388. 
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure, trimr
source("EMTSUtil.R")

#
#---------------  Estimating Drift and Diffusion Functions ------------------
#
 
t <- 7500

kappa <- 0.21459
theta <- 0.085711
sigma <- 0.07830
delta <- 1/250

# Simulate the interest rate 
r <- theta + rep(0, t)

for (i in 2:t) {
  r[i] <- r[i-1] + kappa*(theta - r[i-1])*delta + sigma*sqrt(r[i-1])*sqrt(delta)*rnorm(1)
}

   
# Nonparametric regression of the mean
dr <- trimr(r, 1, 0) - trimr(r, 0, 1)
yt <- dr/delta
xt <- trimr(r, 0, 1)
    

h <- 0.023
x <- seq(0.0, 0.2,0.002)

fx <- rep(0, length(x))
fxy <- rep(0, length(x))
          
for (i in seq(x)) {
  z      <- ((x[i] - xt)/h)    
  fx[i]  <- mean( dnorm(z)/h )
  fxy[i] <- mean( dnorm(z)*yt/h )  
}
mx  <- fxy/fx

yt_mean   <- yt
mx_mean   <- mx
true_mean <- kappa*(theta - x)

# Nonparametric regression of the variance 
dr <- trimr(r, 1, 0) - trimr(r, 0, 1)
yt <- dr^2/delta
xt <- trimr(r, 0, 1)
    
h <- 0.023

fx <- rep(0, length(x))
fxy <- rep(0, length(x))
for (i in seq(x)) {
  z      <- ((x[i] - xt)/h)    
  fx[i]  <- mean( dnorm(z)/h )
  fxy[i] <- mean( dnorm(z)*yt/h )  
}
mx  <- fxy/fx

yt_var   <- yt
mx_var   <- mx
true_var <- sigma^2*x


#**************************************************************************
#**
#**     Generate graphs
#**
#**************************************************************************

figure()
par(mfrow=c(2,2))


#--------------------------------------------------------#
# Panel (a)
plot(yt_mean, type="l",
     main = expression(bold( paste("(a) y"[t]*" =(r"[t]*" - r"[t-1]*")/", Delta))),
     xlab = expression(t),
     ylab = expression(y[t]),
     bty="l")

#--------------------------------------------------------#
# Panel (b)
plot(x,true_mean,type="l",
     main = "(b) Estimated Mean",
     xlab = expression(x[t]),
     ylab = expression(m(x[t])),
     bty="l")
lines(x,mx_mean, lty=2, col="blue")

#--------------------------------------------------------#
# Panel (c)
plot(yt_var,type="l",
     main = expression(bold(paste("(c) y"[t]*" =(r"[t]*" - r"[t-1]*")^2/", ,Delta))),
     xlab = expression(t),
     ylab = expression(y[t]),     
     bty="l")
#--------------------------------------------------------#
# Panel (d)
plot(x,true_var,type="l",
     main = "(c) Estimated Variance",
     xlab = expression(x[t]),
     ylab = expression(m(x[t])),     
     bty="l")
lines(x,mx_var,lty=2, col="blue")
