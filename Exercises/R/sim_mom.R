# ===========================================================================
#
#     Method of moments estimation of a first order MA model. 
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, trimr
source("EMTSUtil.R")

#
#------------------------- Moments of MA model ------------------------------
#

sim_mom <- function() {
  
t     <- 250        # Sample size             
theta <- 0.5        # MA1 Population parameter

# Generate the data for a MA(1) model      
u  <- rnorm(t)
y  <- trimr(u, 1,0) - theta* trimr(u, 0, 1)

# Estimate the first order autocorrelation coefficient
y   <- y - mean(y)
rho <- lm(trimr(y, 1, 0) ~ trimr(y, 0, 1) - 1)$coef

# Estimate theta using the method of moments estimator   
b_mom <- ( -1 + sqrt(1 - 4*rho^2) ) / (2*rho)


cat('\n')
cat('\nSample size                        = ', t)
cat('\nTrue population parameter (theta)  = ', theta)
cat('\nMethod of moment estimate          = ', b_mom)
cat('\n')
cat('\nTrue population parameter (AR(1))  = ', -theta/(1+theta^2))
cat('\nFirst order AR(1)                  = ', rho)


#**************************************************************************
#**
#**     Generate graph
#**
#**************************************************************************

figure()
par(xaxs="i", yaxs="i")
plot(seq(t-1),y,type="l",
     main = "Simulated data from an AR(1) model",
     xlab = expression(t),
     ylab = expression(y[t]),
     xlim = c(0, t),
     ylim = c(-4, 4),
     bty = "l")
}


