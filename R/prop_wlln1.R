#=========================================================================
#
#   Program to demonstrate the Law of Large Numbers 
#   (Exponential distribution example) 
#
#=========================================================================

rm(list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(4, kind="Mersenne-Twister")

mu <- 5      
t  <- 500

# Generate exponential random numbers from a Gamma distribution
e <- rgamma(t, 1)
y <- mu * e           

# Generate sample means from sample of size t = 1,2,3,...,tmax
ybar <- array(0, c(t, 1))

for(i in seq(t))
{
  ybar[i] <- mean(y[1:i])
}


#**************************************************************************
#**
#**     Generate graph
#**
#**************************************************************************

figure()

tt <- seq(1, t, 1)
plot(tt, ybar, type="l",
     col="darkblue",
     main = "Weak Law of Large Numbers (Necessary Condition)",
     ylab = expression(bar(y)),
     ylim = c(2, 8),     
     xlab = "T",
     xlim = c(0, t),
     bty="l")
lines(tt, mu * rep(1, t))
lines(tt, 4.80 * rep(1, t), lty=2)
lines(tt, 5.20 * rep(1, t), lty=2)
