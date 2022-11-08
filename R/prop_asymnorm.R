#=========================================================================
#
#     Monte Carlo program to demonstrate the asymptotical normality
#     of MLE with R = 5000 replications based on an exponential distribution
#     with theta = 1.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(1234, kind="Mersenne-Twister")

R     <- 5000                           
theta <- 1


#     For the sample size of T = 5      
t <- 5                                        #   Sample size 
u <- matrix(runif(t*R), nrow = t)             #   Generate uniform random numbers         
y <- -theta * log(1 - u)                      #   Generate realizations of y              

ybar <- apply(y, 2, mean)
z1 <- sqrt(t) * (ybar - theta)/sqrt(theta)  #   Standardized random variable    


#     For the sample size of T = 100      
t <- 100                            #  Sample size                     
u <- matrix(runif(t*R), nrow = t)             #   Generate uniform random numbers         
y <- -theta * log(1 - u)                      #   Generate realizations of y              

ybar <- apply(y, 2, mean)                     
z2   <- sqrt(t) * (ybar - theta) / sqrt(theta)  #   Standardized random variable     


#**************************************************************************
#**
#**     Generate graph
#**
#**************************************************************************
figure()
layout(matrix(c(1,1,2,3), nrow = 2, ncol=2, byrow = TRUE))
s <- seq(0, 12, 0.001)

fy <- exp(-s / theta) / theta
plot(s, fy, type="l",     
     main = "(a) Exponential distribution",
     ylab = expression(f(y)),       
     xlab = expression(y),    
     xlim = c(0, 6),
     ylim = c(0, 1.5),
     bty="l")

hist(z1, breaks = 21,
     main = "(b) T = 5",
     xlab = expression(z),
     ylab = expression(f(z)))


hist(z2, breaks = 21,
     main = "(c) T = 100",
     xlab = expression(z),
     ylab = expression(f(z)))
