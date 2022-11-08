#=========================================================================
#
#     Program to demonstrate the Lindberg-Levy central limit theorem
#     using the uniform distribution.
#
#=========================================================================


rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(12, kind="Mersenne-Twister")

t      <- 10
Ndraws <- 5000                            # Number of draws                                     
mu     <- 1/2                             # Population mean of the uniform distribution         
sig2   <- 1/12                            # Population variance of the uniform distribution     
u      <- matrix(runif(t*Ndraws), nrow=t) # Generate uniform random numbers                 

# For the sample size of 2        
y     <- u[1:2, ]                          
ybar1 <- colMeans(y)
z1    <- sqrt(2)*(ybar1 - mu)/sqrt(sig2)    # Standardized variable 

# For the sample size of 10          
y     <- u[1:10, ]
ybar2 <- colMeans(y)
z2    <- sqrt(10)*(ybar2 - mu)/sqrt(sig2)    # Standardized variable 


#**************************************************************************
#**     
#**     Generate graph       
#**                                         
#**************************************************************************

figure()
par(mfrow=c(2,2))

hist(z1, breaks=21,     
     main = expression(bold("(a) Distribution of z (T = 2)")),
     ylab = expression(f(z)),       
     xlab = expression(z), 
     bty="l")

hist(z2, breaks=21,     
     main = expression(bold("(b) Distribution of z (T = 10)")),
     ylab = expression(f(z)),       
     xlab = expression(z), 
     bty="l")

hist(ybar1, breaks=21,     
     main = expression(bold(paste("(c) Distribution of ", bar(y), " (T = 2)") ) ),
     ylab = expression(paste("f(", bar(y), ")") ),
     xlab = expression(bar(y)),
     bty="l")

hist(ybar2, breaks=21,     
     main = expression(bold(paste("(d) Distribution of ", bar(y), " (T = 10)") ) ),
     ylab = expression(paste("f(", bar(y), ")") ),
     xlab = expression(bar(y)),
     bty="l")
