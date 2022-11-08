#=========================================================================
#**
#**     Program to demonstrate multiple roots of the bivariate normal model
#**
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

# This commented out code generates the data from scratch
# set.seed(1234, kind="Mersenne-Twister")
# t   <- 4                              # Sample size                                 
# rho <- 0.6
# u1  <- rnorm(t)                     # Generate independent N(0,1)                 
# u2  <- rnorm(t)
# x   <- u1                             # Generate dependent normal random numbers   
# y   <- rho*u1 + sqrt(1 - rho^2)*u2

# Take the GAUSS data for consistency
x <-   c(-0.60303846, -0.098331502, -0.15897445, -0.65344600)

y <-   c(0.15367001, -0.22971064, 0.66821992, -0.44328369)
   
sxy <- mean(x*y)
sxx <- mean(x^2)
syy <- mean(y^2)
zeros <- rep(0, 199)
A   <- zeros                                # Average log-likelihood                    
G   <- zeros                                # Gradients                                 
rho <- seq(-0.99, 0.99, 0.01)               # Grid values of rho      

tg  <- length( rho )

for (i in seq(tg)) {

    A[i] <- -log(2*pi) - 0.5*log(1 - rho[i]^2) - 0.5*(sxx - 2*rho[i]*sxy + syy)/(1 - rho[i]^2)
    G[i] <- rho[i]*(1 - rho[i]^2) + (1 + rho[i]^2)*sxy - rho[i]*(sxx + syy)
}


#**************************************************************************
#**
#**     Generate graph
#**
#**************************************************************************

figure()
par(mfrow=c(1,2), xaxs="i", yaxs="i")


plot(rho, G, type="l",     
     main = "(a) Gradient",
     ylab = expression(G(rho)),       
     xlab = expression(rho),    
     xlim = c(-1, 1),
     ylim = c(-0.6, 0.6),
     bty="l")
lines(rho, rep(0, length(rho)), lty=1)

plot(rho[5:length(rho)-5], A[5:length(A)-5], type="l",     
     main = "(b) Average log-likelihood",
     ylab = expression(A(rho)),       
     xlab = expression(rho),    
     xlim = c(-1, 1),
     ylim = c(-3, -1.5),
     bty="l")

