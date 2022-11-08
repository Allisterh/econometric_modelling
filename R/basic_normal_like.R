#=========================================================================
#
#   Program to estimate the parameters of a normal distribution by
#   maximum likelihood and plot the log of the likelihood function.
#
#=========================================================================

# clear all
rm(list = ls(all = TRUE))
graphics.off()

# load any compatibility requirements
source("EMTSUtil.R")

# 
# Read in data
y <- c(1,2,5,1,2)
t <- length(y)         # Define the sample size  

# Compute the MLEs
mu_mle   <- mean(y)
sig2_mle <- mean( (y - mu_mle)^2 )

cat('\nMLE(mu)\n')
cat('-------\n')
cat(mu_mle, '\n')

cat('\nMLE(sig2)\n')
cat('-------\n')
cat(sig2_mle, '\n')

#**************************************************************************
#**
#**     Generate graph
#**
#**************************************************************************

mu   <-  seq(1.5, 3.0, 0.05)
sig2 <-  seq(1.5, 3.0, 0.05)

nrows <- length(mu)
ncols <- length(sig2)
lnl <- array(0, c(nrows, ncols))

for(i in seq(nrows))
{
  for(j in seq(ncols))
  {
    lnl[i,j] <- -0.5*t*log(2*pi) - 0.5*t*log(sig2[j]) - sum( ( (y - mu[i])^2 ) / sig2[j] )
  }
}
figure()
persp(sig2, mu, lnl, 
      theta=40, 
      phi=15, 
      expand=0.5,      
      ticktype='detailed',      
      col="blue",      
      ylab = "mu",
      xlab = "sigma^2",
      zlab = "ln LT((mu, sigma^2)")
