#============================================================================
#
#     Normal distribution example using different parametric distributions.
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
#--------------------  Parametric Approximations ----------------------------
#

t <- 200        #      Sample size 
mue <- 0
sig <- 3

yt <- mue + sig*rnorm(t)

# Generate the population (normal) distribution       
y <- seq(-15, 15, 0.1)

f_norm <- dnorm( ( y - mue )/sig )/sig

# Estimate the parametric density assuming normality   
m <- mean(yt)
s <- sd(yt)
f_npara <- dnorm( ( y - m )/s )/s

# Estimate the parametric density assuming Student t       
m <- mean(yt)
s <- sd(yt)
v <- 2*s^2/(s^2-1)

f_spara <- (1+((y-m)/s)^2/v)^(-(v+1)/2)*gamma((v+1)/2)/(sqrt(pi*v)*gamma(v/2)*s)


#****************************************************************************
#
#     Generate graphs
#
#****************************************************************************

figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

#--------------------------------------------------------#
# Panel (a)
plot(yt, type="l",
     main = "(a) Data",
     xlab = expression(t),
     ylab = expression(y[t]),     
     bty="l")

#--------------------------------------------------------#
# Panel (b)
plot(y, f_norm, type="l",
     main = "(b) Distribution",
     xlab = expression(y),
     ylab = expression(f(y)),     
     bty="l")

figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

#--------------------------------------------------------#
# Panel (a)
plot(y,f_norm,type="l",
     main = "(a)",
     xlab = expression(y),
     ylab = expression(f(y)),     
     ylim = c(0, 0.15),
     bty="l")
lines(y, f_npara, lty=2)

#--------------------------------------------------------#
# Panel (b)
plot(y,f_norm,type="l",
     main = "(b)",
     xlab = expression(t),
     ylab = expression(y[t]),  
     ylim = c(0, 0.15),
     bty="l")
lines(y,f_spara, lty=2)

