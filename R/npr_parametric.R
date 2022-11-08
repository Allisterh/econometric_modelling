#============================================================================
#
#    Parametric solutions of the nonlinear example 
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure, inv
source("EMTSUtil.R")

#
#---------------  Parametric Conditional Mean Estimator ---------------------
#
  

# Simulate the model
t  <- 500
ut <- 0.1*rnorm(t)                # N(0,0.1^2) 
xt <- -2 + 4*runif(t)              # U(-2,2) 
mx <- 0.3*exp( -4*(xt + 1)^2 ) + 0.7*exp( -16*(xt - 1)^2 )  
yt  <- mx + ut


# Estimate a linear model
y  <- yt
x  <- cbind(rep(1, t), xt)
b1 <- inv(t(x) %*% x) %*% t(x) %*% y

yfit1 <- x %*% b1


# Estimate a nonlinear (polynomial) model
y  <- yt
x  <- cbind(rep(1,t), xt, xt^2, xt^3, xt^4)
b2 <- inv( t(x) %*% x) %*% t(x) %*% y

yfit2 <- x %*% b2



#**********************************************************************
#***
#***     Generate graphs
#***
#**********************************************************************

figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

tmp <- cbind(xt, mx, yfit1, yfit2)
# sort by rows
tmp <- tmp[order(tmp[,1]), ]

#--------------------------------------------------------#
# Panel (a)
plot(tmp[,1],tmp[,2],type="l",
     main = "(a) Linear",
     ylab = expression(m(x[t])),
     xlab = expression(x[t]),
     ylim = c(0, 0.7),
     xlim = c(-2, 2))
lines(tmp[,1],tmp[,3], lty=2)


#--------------------------------------------------------#
# Panel (b)
plot(tmp[,1],tmp[,2],type="l",
     main = "(b) Nonlinear",
     ylab = expression(m(x[t])),
     xlab = expression(x[t]),
     ylim = c(-0.3, 0.7),
     xlim = c(-2, 2))
lines(tmp[,1],tmp[,4], lty=2)

