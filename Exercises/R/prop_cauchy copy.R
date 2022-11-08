#=========================================================================
#
#     Program to demonstrate the consistency property of MLE for the
#     location parameter (theta) of the Cauchy distribution.
#
#     For this example the median is the MLE and is thus a consistent estimator
#     but the sample mean is an inconsistent estimator.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(123, kind="Mersenne-Twister")


theta  <- 1                         # Location parameter 
nu  <- 1                            # Cauchy = Student t nu =1 dof
t <- 500

# Generate sample means and medians from sample of size t=1,2,3,...,t        

zeros <- rep(0, t)
ybar <- zeros
ymed <- zeros

for (i in seq(t) ) {  
   y <- theta + rt(i, nu)
   ybar[i] <- mean(y)
   ymed[i] <- median(y)
}


#**************************************************************************
#**     
#**     Generate graph      
#**                                         
#**************************************************************************

figure()
layout (matrix(c(1,2), byrow = TRUE, ncol = 2))

# generate line graph
plot(1:t, ybar, type="l",
     col="darkblue",
     main = "(a) Mean",
     ylab = expression(hat(theta)),
     xlab = "Progressive Sample Size",
     bty="l")

plot(1:t, ymed, type="l",
     col="darkblue",
     main = "(b) Median",
     ylab = expression(hat(theta)),
     xlab = "Progressive Sample Size",
     bty="l")
lines(1:t, theta * rep(1, t))

