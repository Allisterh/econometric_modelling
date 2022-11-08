#============================================================================
#
#     Nonaparametric density of ftse returns 
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")

#
#-----------------------  KD FTSE Share Index -------------------------------
#

# Load equities returns data on ftse 20 Nov. 1973 to 23 July 2001 

yt <- as.matrix(read.table("ftse.dat"))
t  <- length(yt)
y  <- seq(-0.10, 0.10, 0.001)


# Generate the Student t distribution for comparison
mu  <- mean(yt)
sig <- apply(yt, 2, sd)
nu  <- 5

f_studt <- dt( (y - mu)/sig,nu )*(1/sig)

# Estimate the nonparametric density
f <- rep(0, length(y) )
h <- 1.06*sig*t^(-1/5)

for (i in seq(y)) {
  z     <- (y[i] - yt)/h    
  f[i] <- mean( dnorm(z)/h )  
}

#**********************************************************************
#***
#***     Generate Graphs
#***
#**********************************************************************
figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

#--------------------------------------------------------#
# Panel (a)
plot(yt,type="l",
     main = "(a) Data",
     xlab = expression(t),
     ylab = expression(y[t]),
     xlim = c(0, 7000),
     ylim = c(-0.12, 0.08),
     bty = "l")

#--------------------------------------------------------#
# Panel (b)
plot(y,f_studt, type="l",
     main = "(b) Distributions",
     xlab = expression(y),
     ylab = expression(y[t]),
     ylim = c(0, 60),
     bty="l")
lines(y, f, lty=2, col="blue")
