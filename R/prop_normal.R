#=========================================================================
#
#     Program to demonstrate the consistency property of MLE for the
#     mean of the normal distribution.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(123, kind="Mersenne-Twister")

t    <- 500    
mu   <- 1               # Population mean
sig2 <- 2               # Population variance  

# Generate sample means from sample of size t=1,2,3,...,t       

yb  <- rep(0, t)
for (t in seq(t)) {
  y <- mu + sqrt(sig2) * rnorm(t)
  yb[t] <- mean(y)  
}
    

#**************************************************************************
#**     
#**     Generate graph     
#**                                         
#**************************************************************************

figure()

plot(1:t, yb, type="l",     
     main = "Consistency of the sample mean assuming normality",
     ylab = expression(bar(y)),       
     xlab = "T",         
     bty="l")
lines(1:t, mu*rep(1, t))
