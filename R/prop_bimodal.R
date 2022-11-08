#============================================================================
#
#     Program to demonstrate multiple roots of the bivariate normal model
#
#============================================================================
rm(list = ls(all=TRUE))
graphics.off()

#set.seed(1234, kind="Mersenne-Twister")


# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - figure, seqa
source("EMTSUtil.R")

# 
#------------------------- Bimodal likelihood ------------------------------
#
# t   <- 4                              # Sample size                                 
# rho <- 0.6
# u1  <- rnorm(t)                     # Generate independent N(0,1)                 
# u2  <- rnorm(t)
# x   <- u1                             # Generate dependent normal random numbers   
# y   <- rho*u1 + sqrt(1 - rho^2)*u2

# Take the GAUSS data for consistency
x <-   c(-0.60303846,
       -0.098331502,
       -0.15897445, 
       -0.65344600 ) 

y <-   c(0.15367001,
        -0.22971064, 
        0.66821992, 
        -0.44328369)



sxy <- mean(x*y)
sxx <- mean(x^2)
syy <- mean(y^2)
lnl <- rep(0, 199)                  # log-likelihood         
g   <- rep(0, 199)                  # gradients                                 
rho <- seq(-0.99, 0.99, 0.01)                # grid values of rho      

tg <- length( rho )

for (i in seq(tg)) {  
  lnl[i] = -log(2*pi) - 0.5*log(1 - rho[i]^2) - 0.5*(sxx - 2*rho[i]*sxy + syy)/(1 - rho[i]^2)
  g[i]   = rho[i]*(1 - rho[i]^2) + (1 + rho[i]^2)*sxy - rho[i]*(sxx + syy)
  
}
 

#**************************************************************************
#**
#**     Generate graph
#**
#**************************************************************************
figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2)) 
 
matplot(rho,cbind(g, rep(0, tg)), type="l",
        main = '(a) Gradient',
        ylab = expression(G[T](rho)),
        xlab = expression(rho),
        xlim = c(-1, 1),
        ylim = c(-0.6, 0.6),
             bty = "l")


plot(rho[5:length(rho)-5],lnl[5:length(lnl)-5], type="l",
     main = '(b) Log-likelihood function',
     ylab = expression(L[T](rho)),
     xlab = expression(rho),
     xlim = c(-1, 1),
     ylim = c(-3, -1.5),
     bty = "l")
