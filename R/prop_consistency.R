#=========================================================================
#
#   Program to demonstrate two aspects of consistency
#
#=========================================================================


rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(66, kind="Mersenne-Twister")

mu   <- seq(0, 20, 0.1)
sig  <- 4.0

zeros <- rep(0, length(mu))

lnl1 <- zeros
lnl2 <- zeros
lnl3 <- zeros


# Population parameters 
mu_0  <- 10.0
sig_0 <-  4.0


# Sample T=5 
t <- 5
y <- mu_0 + sig_0 * rnorm(t)

cat('\nSample mean (T=5)   =', mean(y), '\n')


for (i in seq(mu)){
  lnl1[i] <- mean( -0.5*log(2*pi*sig^2) - 0.5*(y-mu[i])^2/sig^2 )
}
  

# Sample T=20
t <- 20
y <- mu_0 + sig_0 * rnorm(t)

cat('\nSample mean (T=20)   =' , mean(y), '\n')

for (i in seq(mu)){
  lnl2[i] <- mean( -0.5*log(2*pi*sig^2) - 0.5*(y-mu[i])^2/sig^2 )
}



# Sample T=500
t <- 500
y <- mu_0 + sig_0 * rnorm(t)

cat('\nSample mean (T=500)   =', mean(y), '\n')

for (i in seq(mu)){
  lnl3[i] <- mean( -0.5*log(2*pi*sig^2) - 0.5*(y-mu[i])^2/sig^2 )
}


# Compute population log-likelihood  

e_lnl <- -0.5*log(2*pi*sig_0^2) - 0.5 - 0.5*(mu - mu_0)^2/sig_0^2

tmp <- -0.5*(log(2*pi*sig_0^2) + 1)
cat('\nMaximum value at theta = theta_0 =', tmp, '\n')



#********************************************************************
#***
#***     Generate graph
#***
#********************************************************************

figure()

plot(mu, e_lnl, type="l",     
     main = "Demonstration of consistency property",
     ylab = expression(paste("ln L"[T]*"","(", mu, ")")),       
     xlab = expression(mu),    
     xlim = c(3, 17),
     ylim = c(-3.5, -2.5),
     bty="l", 
     lty = 1)

lines(mu, lnl1, lty=2)
lines(mu, lnl2, lty=4)
lines(mu, lnl3, lty=6)

# Add a legend to the plot  
legend("topright",                       
      legend=c("Population log-likelihood", 
                "log-likelihood T = 5",
                "log-likelihood T = 20",
                "log-likelihood T = 500"),             
      lty=c(1, 2, 4, 6),                    
      lwd=c(1,1,1,1))

