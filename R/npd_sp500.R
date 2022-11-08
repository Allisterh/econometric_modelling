#============================================================================
#
#     Empirical nonparametric examples (oznyse)
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions - figure
source("EMTSUtil.R")
# Load required library - repmat
library("matlab")

#
#-----------------------  KD S&P Share Index -------------------------------
#

# Load equities data on 2&p500 9 Feb. 1954 to 23 July 2001
yt <- as.matrix(read.table("sp500.dat"))
t  <- length(yt)

# Generate the Student t distribution for comparison   
mue <- mean(yt)
sig <- apply(yt, 2, sd)
nue <- 5

y <- seq(-0.25, 0.25, 0.001)
z <- ( y - mue )/sig

f_studt <- gamma((nue+1)/2)*(1/gamma(0.5))*(1/gamma(nue/2))*(1/sqrt(nue))*(1 + z^2/nue)^(-(1+nue)/2)*(1/sig)


# Estimate the nonparametric density
h  <- 1.06*sig*t^(-1/5)
wt <- t( dnorm((repmat(y,1,t) - repmat(yt,501,1))/h)/h )
f  <- colMeans(wt) 

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
     ylab = expression(y[t]),
     xlab = expression(t),
     xlim = c(0, 12000),
     bty="l")

#--------------------------------------------------------#
# Panel (b)
plot(y,f,type="l", lty=2,
     xlab = expression(y),
     ylab = expression(f(y)),     
     bty="l")
lines(y, f_studt)
