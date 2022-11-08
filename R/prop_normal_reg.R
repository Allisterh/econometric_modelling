#=========================================================================
#
#    Program to demonstrate the convergence properties of the average
#    log likelihood to its expectation assuming a normal distribution.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(12345, kind="Mersenne-Twister")

t    <- 1000    
mu   <- 1               # Population mean
sig2 <- 2               # Population variance  

# Generate data from a normal distribution  
y <- mu + sqrt(sig2) * rnorm(t)

# Compute progressive sample means of the log-likelihood 
lnf <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*(y - mu)^2/sig2
tt  <- seq(t)

s_lnf <- cumsum(lnf)/tt
e_lnf <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5

cat('\nPopulation mean')
cat('\n----------------\n')
cat(e_lnf)

cat('\n\nSample mean')
cat('\n-----------\n')
cat(mean(s_lnf))


#**************************************************************************
#**     
#**     Generate graph       
#**                                         
#**************************************************************************

figure()

plot(tt, s_lnf, type="l",     
     ylab = expression(A(theta)),       
     xlab = expression(T),
     bty="l")

lines(tt, e_lnf*rep(1,t), lty=2)
