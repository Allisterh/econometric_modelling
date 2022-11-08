#============================================================================
#
#   Program to compute gradient, Hessian, information matrix and OPG
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

set.seed(12, kind="Mersenne-Twister")


# 
#------------------------- Properties of gradient function ------------------
#
t    <- 5000         

# Normal distribution  
mu   <- 1                  
sig2 <- 1                  
u    <- rnorm(t)
y    <- mu  + sqrt(sig2 )*u     
m    <- mean(y)              # maximum likelihood estimate      
g    <- (y - m)/sig2         # gradient at t    
h   <- -1/sig2               # hessian at t    
j   <- t(g) %*% g                  # opg  

cat('\n Normal distribution results ')
cat('\nGradient vector    = ',mean(g))
cat('\nHessian            = ',mean(h))
cat('\nInformation matrix = ',-mean(h))
cat('\nOPG matrix         = ',j/t)
cat('\n ')

# Exponential distribution 
theta <- 2               

y  <- theta*rgamma(t, shape=1, scale=1)
th <- 1/mean(y)             # maximum likelihood estimate
g  <- 1/th - y              # gradient at t 
h  <- -1/th^2               # hessian at t
j <- t(g) %*% g                   # opg

cat('\n Exponential distribution results ')
cat('\nGradient vector    = ',mean(g)) 
cat('\nHessian            = ',mean(h))
cat('\nInformation matrix = ',-mean(h)) 
cat('\nOPG matrix         = ',j/t)
cat('\n ')
