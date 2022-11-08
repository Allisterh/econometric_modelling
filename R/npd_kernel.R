#============================================================================
#
#      Gaussian kernel example
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")
# load required functions - repmat
library("matlab")
#
#-------------------- Kernel Density Estimator ------------------------------
#

xi <- c(2, 3, 5, 7, 8, 8, 9, 9, 10, 10, 11, 11, 11, 12, 12, 14, 15, 17, 18, 20)
n  <- length(xi)

# Compute kernel with a window of h = 1       

h <- 1.0
x <- c(5, 10, 15)
 
w1 <-  t( (dnorm(repmat(x,1,n)-repmat(xi,3,1))/h) )
f1 <- colSums(w1)/(n*h) 

# Compute kernel with a window of h = 2 
h <- 2.0
 
w2 <-  t( (dnorm(repmat(x,1,n)-repmat(xi,3,1))/h) )
f2 <- colSums(w2)/(n*h)
print(cbind(w1, w2))

  
figure()  
matplot(x,cbind(f1, f2),type="l",
     ylab="",
     xlab="",
    lty = c(1,2))
  

