#=========================================================================
#
#     Program to plot the profile likelihood for the portfolio
#     diversification model
#     
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

#
#-------------------------------- Helper Functions ---------------------------------
#

# 
# Calculates the log-likelihood function for a bivariate normal distribution
#
lnl <- function(th, y) {  
  n <- ncol(y)
  
  mu1 <- th[1]        # Mean of 1
  mu2 <- th[2]        # Mean of 2    
  s11 <- th[3]        # Variance 1
  s22 <- th[4]        # Variance 2  
  rho <- th[5]        # Correlation of 1 and 2
  
  z1 <- (y[, 1] - mu1)/sqrt(s11)
  z2 <- (y[, 2] - mu2) /sqrt(s22)
  
  lnfn <- -0.5*n*log(2*pi) - 0.5*(log(s11) + log(s22) + log(1 - rho^2) ) - (1/(2*(1 - rho^2))) * ( z1^2 - 2*rho12*z1*z2 + z2^2 )  
  return ( sum(lnfn) )
}

#
#------------------------------ Profile Log-likelihood method ------------------------------
#
# Asset price data from 6 August 2010 to 2 January 2001 
#(note that the data are in reverse order ie from recent to past  

# Load data
load("diversify.RData")

# Select appropriate sample  
pt_apple     <- as.matrix(pt_apple)[1:2413]
pt_ford      <- as.matrix(pt_ford)[1:2413]


# Compute percentage returns  
r_apple <- 100*diff(log(pt_apple)) 
r_ford  <- 100*diff(log(pt_ford))


#     Compute maximum likelihood estimates 
m1 <- mean(r_apple)
m2 <- mean(r_ford)

s11 <- mean((r_apple - mean(r_apple))^2)
s22 <- mean((r_ford - mean(r_ford))^2)

rho12 <- cor(r_apple,r_ford)
cat('\nSample correlation          = ', rho12 , '\n')

#     Generate the profile log-likelihood 
y <- cbind(r_apple, r_ford)

r <- seq( -0.9, -0.9+0.01*(186-1), 0.01 )

a <- rep(0, length(r))

mean <- c(m1, m2)

for (i in seq(r)) {
  theta <- c(m1, m2, s11, s22, r[i])
  
  # Compute average log-likelihood value for alternative values of rho
  a[i] <- lnl(theta, y)/nrow(y)
}   

#********************************************************************
#***
#***     Generate graph
#***
#********************************************************************

figure()

#     Plot the profile log-likelihood
par(xaxs="i")
plot(r,a, type="l",    
     xlab = expression(rho),
     ylab = "Average lnl",
     xlim = c(-1, 1))

