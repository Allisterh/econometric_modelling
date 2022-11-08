#=========================================================================
#
#   Program to estimate an exponential model and plot the 
#   log-likelihood function.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# load any compatibility requirements
source("EMTSUtil.R")

# Sample data for y  
y <- c(2.1, 2.2, 3.1, 1.6, 2.5, 0.5)

t     <- length(y)           
theta <- 1/mean(y)
lnlt  <- log(theta) - theta * y
gt    <- 1/theta - y
ht    <- -1/theta^2

cat('\nSum of y  = ', sum(y), '\n')
cat('Mean of y = ', mean(y), '\n')
cat('MLE       = ', theta, '\n')
cat('Log-like  = ', mean(lnlt), '\n')
cat('Gradient  = ', mean(gt), '\n')
cat('Hessian   = ', mean(ht), '\n')
    
rnames <- seq(t)
cnames <- c("yt", "lnlt", "gt", "ht")
ones <- array(1, c(t,1))
print(matrix(c(y, lnlt, gt, ones*ht), 
             dimnames=list(rnames, cnames), nrow=t))

# *************************************************************************
# ***
# ***     Generate graphs
# ***
# *************************************************************************

theta <- seq(0.001, 1.5, 0.01)
lnl <- t*log(theta) - theta * sum(y)

# create subplots
figure()
layout(matrix(c(1,2), byrow=TRUE, ncol=2))
plot(theta, lnl, type="l",
     col="darkblue",
     main = "(a) Log-likelihood function",
        ylab = expression(paste("ln L"[T]*"","(",theta, ")")),
        xlab = expression(theta),
        bty="l")

plot(theta, exp(lnl)* 1e5, type="l",
     col="darkblue",
     main = "(b) Likelihood function",
     ylab = expression(paste("ln L"[T]*"", "(",theta, ")", x10^{5})),
     xlab = expression(theta),
        bty="l")

