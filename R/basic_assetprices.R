#=========================================================================
#
#   Program to compute the maximum likelihood estimates
#   of the asset price model.
#
#   The data consist of the Australian, Singapore and NASDAQ stock
#   market indexes for the period 3 January 1989 to 31 December 2009.
#
#=========================================================================

rm(list = ls(all = TRUE))
graphics.off()

# Load data
load("assetprices.RData")

pt <- as.matrix(pt_aus)     # Choose the asset price 
pt <- pt/1000;

# MLE based on the log-normal distribution of the share index 
# Transitions are from log price at t-1 to log price at t
y   <- log( pt[-1] ) - log( pt[ -length(pt) ] )
  
alpha <- mean(y)
sig2  <- mean( (y - alpha)^2 )

cat('\nMLE of alpha based on the log-normal distribution = ', alpha, '\n')
cat('\nMLE of sig2 based on the log-normal distribution  = ', sig2, '\n')
