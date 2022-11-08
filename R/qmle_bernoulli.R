#===============================================================================
#
#   Program to compute the QMLE standard errors where the true 
#   distribution is Bernoulli and the misspecified distribution is normal
#
#===============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")

#
#------------------------- Helper Functions -----------------------------------
#

#load required functions - inv
source("EMTSUtil.R")

#
#--------------------  QMLE Bernoulli Distribution ----------------------------
#

mu <- 0.6         
t  <- 100         

#  Generate data from the true model (Bernoulli)  
y <- runif(t) < mu                  
y <- as.numeric(y)

# Estimate the true model (Bernoulli) 
theta0 <- mean(y)

# Gradient of the true model (at each observation)  
g0 <- y/theta0 - (1-y)/(1-theta0)

# Hessian of the true model                             
h0 <- mean( -y/theta0^2 - (1-y)/(1-theta0)^2 )

# Information of the true model                       
i0 <- 1/(theta0*(1-theta0))


# OPG of the true model 
j0 <- t(g0) %*% g0/t

# Estimate the misspecified model (normal)    
theta1  <- mean(y)

# Gradient of the misspecified model (at each observation)     
g1 <- y - theta1

# Hessian of the misspecified model                       
h1 <- -1

# Information of the misspecified model   
i1 <- 1


# OPG of the misspecified model (independence)    
j1 <- t(g1) %*% g1/t


cat('\nCovariance matrix true (i)            = ', inv(i0)/t)
cat('\nCovariance matrix true (h)            = ', inv(-h0)/t)
cat('\nCovariance matrix true (j)            = ', inv(j0)/t)
cat('\n')
cat('\nCovariance matrix misspecified (i)    = ', inv(i1)/t )
cat('\nCovariance matrix misspecified (h)    = ', inv(-h1)/t)
cat('\nCovariance matrix misspecified (j)    = ', inv(j1)/t)
cat('\nCovariance matrix misspecified (qmle) = ', inv(i1)*j1*inv(i1)/t)
