#=========================================================================
#
#   Program to compute the maximum likelihood estimates of the portfolio
#   diversification model
#
#   Asset price data from 6 August 2010 to 2 January 2001 (note that the
#   data are in reverse order ie from recent to past)
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

# Load data
load("diversify.RData")

# Select appropriate sample  
pt_apple     <- as.matrix(pt_apple)[1:2413]
pt_ford      <- as.matrix(pt_ford)[1:2413]


# Compute percentage returns  
r_apple <- 100*diff(log(pt_apple)) 
r_ford  <- 100*diff(log(pt_ford))


# Compute statistics  
m1 <- mean(r_apple)
m2 <- mean(r_ford)

s11 <- mean((r_apple - m1)^2)
s22 <- mean((r_ford - m2)^2)
rho12 <- cor(r_apple,r_ford)


cat('\nSample mean (Apple)         = ', m1)
cat('\nSample mean (Ford)          = ', m2)
cat('\nSample variance (Apple)     = ', s11)
cat('\nSample variance (Ford)      = ', s22)
cat('\nSample correlation          = ', rho12 , '\n')


# Compute weights and risk of optimal portfolio
cov12   <- rho12*sqrt(s11*s22) 
w1      <- ( s22 - cov12) / ( s11 + s22 - 2*cov12 )
w2      <- 1 - w1
s2_port <- w1^2*s11 + w2^2*s22 + 2*w1*w2*cov12

cat('\nOptimal weight (Apple)      = ' , w1)
cat('\nOptimal weight (Ford)       = ', w2)
cat('\nRisk of optimal portfolio   = ', s2_port, '\n')


#********************************************************************
#***
#***     Generate graph
#***
#********************************************************************

figure()
plot(r_ford, r_apple, 
     ylab = "Apple",
     xlab = "Ford",
     pch=19,
     bty="l")

