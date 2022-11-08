#============================================================================
#
#     Cauchy nonparametric example.
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")

#
#------------------------  Cauchy Distribution ------------------------------
#

t <- 100         #      Sample size            
r <- 1000        #      Number of replications 
nue <- 1
yt  <- 1/(0.5 - atan(matrix(runif(t*r), nrow=t, ncol=r) )/pi)  # Inverse of cauchy cdf. 

sdyt <- apply(yt, 2, sd)
tstat <- sqrt(t)*colMeans(yt)/sdyt

# Generate the asymptotic (normal) distribution
y <- seq(-50, 49.90, 0.1)
f_norm <- dnorm(y)

# Estimate the nonparametric density 
h  <- 1.06*sd(tstat)*t^(-1/5)
wt <- dnorm(((y - tstat)/h))/h
f  <- mean(wt) 



#****************************************************************************
#
#     Generate graphs
#
#****************************************************************************

figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

#--------------------------------------------------------#
# Panel (a)
plot(tstat, type="l",
     main = "(a) T-stats",
     ylab = expression(lt),
     xlab = expression(lr),     
     bty="l")

#--------------------------------------------------------#
# Panel (b)
plot(y, f_norm, type="l",
     main = "(b) Distribution",
     ylab = expression(f(t)),
     xlab = expression(t),
     xlim = c(-6,6),
     bty = "l")
lines(y,rep(f, length(y)), lty=2, col="red")
