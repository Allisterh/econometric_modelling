#============================================================================
#
#      Program to compute and plot the News Impact Curve 
#
#============================================================================
  
rm(list = ls(all=T))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
# ------------------------ Helper Function ----------------------------------
#
# Load required functions - figure
source("EMTSUtil.R")

#
# ------------------------ Garch NIC model ----------------------------------
#
# Parameter values
a0 <- 0.1
a1 <- 0.5
b1 <- 0.4


# Compute news impact curve of an ARCH(1) model
u <- seq(from=-5, to=5, by=0.1)

a0 <- 1
a1 <- 0.2
sig21 <- a0 + a1*u^2

a0 <- 1
a1 <- 0.5
sig22 <- a0 + a1*u^2

a0 <- 1
a1 <- 0.8
sig23 <- a0 + a1*u^2


#*********************************************************************
#**     Generate graph of news impact curve
#*********************************************************************
figure()

matplot( u,cbind(sig21, sig22, sig23), type="l",
         ylab = expression(h[t]),
         xlab = expression(u),
         bty = "l",
         lty = c(1,3,5),
         col = c("green", "blue", "red"), lwd = 1)

legend("topright", 
       c("alpha=0.2","alpha=0.5", "alpha=0.8"), 
       lty = c(1,3,5),
       col = c("green", "blue", "red"), lwd = 1)
