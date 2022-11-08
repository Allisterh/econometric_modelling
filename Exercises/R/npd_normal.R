#============================================================================
#
#     Normal distribution example using different bandwidths.
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")

#
#-----------------------  Normal Distribution -------------------------------
#

t <- 201
mue <- 0
sig <- 3

yt <- mue + sig*rnorm(t)

# Generate the population (normal) distribution
y <- seq(-10, 10, 0.1)
f_norm <- dnorm( ( y - mue )/sig )/sig


# Estimate the nonparametric density with h = 1
h <- 1
f1 <- rep(0, length(y))
for (i in seq(y)) {
  z     <- (y[i] - yt)/h
  f1[i] = mean( dnorm(z)/h )  
} 

# Estimate the nonparametric density with h = 2 
f2 <- rep(0, length(y))
h  <- 2
for (i in seq(y)) {
  z     <- (y[i] - yt)/h
  f2[i] <- mean( dnorm(z)/h ) 
}
    
# Estimate the nonparametric density with h based on rot  
f3 <- rep(0, length(y))
h  <- 1.06*sd(yt)*t^(-1/5)
for (i in seq(y)) {
    z     <- (y[i] - yt)/h    
    f3[i] <- mean( dnorm(z)/h )
}

# Estimate the parametric density assuming normality   
m <- mean(yt)
s <- sd(yt)
f_para <- dnorm( ( y - m )/s )/s


#****************************************************************************
#
#     Generate graphs
#
#****************************************************************************


figure()
par(xaxs="i", yaxs="i", mfrow=c(1,3))

#--------------------------------------------------------#
# Panel (a)
plot(y,f_para,type="l",
     main = "(a)",
     ylab = expression(f(y[t])),
     xlab = expression(y),     
     bty="l")
lines(y,f_norm, lty=2)
lines(y,f1,lty=3, col="red")

      

#--------------------------------------------------------#
# Panel (b)
plot(y,f_para,type="l",
     main = "(b) h = 2",
     ylab = expression(f(y[t])),
     xlab = expression(y),     
     bty="l")
lines(y,f_norm, lty=2)
lines(y,f2,lty=3, col="red")


#--------------------------------------------------------#
# Panel (c)
plot(y,f_para,type="l",
     main = "(c) h based on rot",
     ylab = expression(f(y[t])),
     xlab = expression(y),     
     bty="l")
lines(y,f_norm, lty=2)
lines(y,f3,lty=3, col="red")
