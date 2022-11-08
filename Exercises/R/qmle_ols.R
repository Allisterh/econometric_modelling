#=============================================================================
#
#   Program to compute the White and Newey-West standard errors 
#   for the OLS estimator
#
#=============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions -----------------------------------
#

#load required functions - inv, trimr
source("EMTSUtil.R")

#
#-----------------------  Linear Regression Model -----------------------------
#

# Data to reproduce the numbers in the chapter 
y <- c(3,2,4,6,10)
x <- c(1,2,3,4,5)

t <- length(y)

# Estimate by OLS
x  <- cbind(rep(1, t), x)
model  <- lm(y ~ x - 1)
u  <- model$residuals
s2 <- mean( u^2 )

# Compute gradients at each t
g <- u * x/s2

# Compute covariance based on Hessian
h <- -t(x) %*% x/s2

cat('\nCovariance matrix (Hessian)\n')
print(unname(-h))
cat('\nStandard errors   (Hessian)\n')
print(unname(sqrt(diag(-h))))

# Compute OPG
j <- t(g) %*% g

cat('\nCovariance matrix (OPG)\n')
print(unname(j) )
cat('\nStandard errors   (OPG)\n')
print(unname(sqrt(diag(j))))

# White covariance matrix
covw <- inv(h) %*% j %*% inv(h)
 
cat('\nCovariance matrix (White)\n')
print(unname(covw))
cat('\nStandard errors   (White)\n')
print(unname(sqrt(diag(covw))))
cat('\n')

# Newey-West covariance matrix
p <- floor( 4*(t/100)^(2/9) )

for (i in seq(1:p)) {
  gmat <- t( g[((i+1):t),] ) %*% g[(1:(t-i)),]
  gmat <- t( trimr(g,i,0) ) %*% trimr(g,0,i)
  j    <- j + (1.0 - i/(p+1))*(gmat + t(gmat))  
}
covnw <- inv(h) %*% j %*% inv(h)

cat('\nCovariance matrix (Newey-West)\n')
print(unname(covnw))
cat('\nStandard errors   (Newey-West)\n')
print(unname(sqrt(diag(covnw))))
cat('\n ')

