#================================================================================
#
#      True model is exponential, misspecified model is normal.
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
#--------------------  QMLE Exponential Distribution -----------------------------
#
t  <- 500000
mu <- 0.5                                          

# Generate data from the true model (Exponential)   
z <- runif(t)                             
y <- -mu*log(1 - z)                       


# Estimate the misspecified model (Normal)   
m <- mean(y)


# Compute gradients of the misspecified model (Normal)
g <- (y - m)


cat('\nGradients ', sum(g))
cat('\n')

# Compute hessian of the misspecified model (Normal) 
h <- -t

cat('\nAverage Negative Hessian ', -h/t )
cat('\n' )

# Compute opg of the misspecified model (Normal)        
cat('\nAverage OPG Matrix ', t(g) %*% g/t) 

 
