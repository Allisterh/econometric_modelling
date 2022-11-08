#=========================================================================
#
#     Program to demonstrate the convergence of the sample log-likelihood
#     function
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(123456, kind="Mersenne-Twister")


t    <- c(5, 20, 100)
mu   <- seq(-3, 3, 0.1)
true <- -0.5*(1 + mu^2)
samp <- array(0, c(length(mu),length(t)))

for (i in seq(t)) {    
    yt   <- rnorm(t[i])
    
    for (j in seq(mu)) {
        samp[j,i] <- -0.5*mean((yt-mu[j])^2)
    }
}

#**************************************************************************
#**     
#**     Generate graph       
#**                                         
#**************************************************************************

figure()

plot(mu, true, type="l",          
     ylab = expression(paste("ln L"[T]*"","(", mu, ")")),       
     xlab = expression(mu),
     bty="l", 
     lty = 1)
lines(mu, samp[,1], lty=2)
lines(mu, samp[,2], lty=4)
lines(mu, samp[,3], lty=6)

# Add a legend to the plot  
legend("topright",                       
      legend=c("True", 
                "t = 5",
                "t = 20",
                "t = 100"),             
      lty=c(1, 2, 4, 6),                    
      lwd=c(1,1,1,1))

