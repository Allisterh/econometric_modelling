#============================================================================
#
#    Gaussian bivariate nonparametric kernel regression example 
#    
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

# load required library - ones, zeros, meshgrid, reshape
library("matlab")

#
#------------------------  Multivariate KD Estimation ------------------------
#
# Simulate the nonlinear model
n <- 2          # Number of explanatory variables
t <- 20000      # Number of observations 

xt <- matrix(rnorm(t*n), nrow=t, ncol=n)
ut <- 0.1*rnorm(t)
yt <- xt[,1]*xt[,2]^2 + ut

# Set up a grid between -3 and 3
xmesh <- seq(-3, 3, 0.2)

# Need to work out all possible points at which to compute bivariate
# density on this mesh and express results as vectors
mesh <- meshgrid(xmesh)
x1m <- mesh$x
x2m <- mesh$y

x <- cbind( c(x1m), c(x2m) )


# Estimate density using product kernel
fac <- 1/(4.0+n)
sdxt <- apply(xt, 2, sd)
h   <- sdxt/(t^fac)
ph  <- prod( h )

ker  <- zeros(n,1)
fx   <- zeros(nrow(x),1)
fxy  <- zeros(nrow(x),1)
pker <- zeros(t,n)

pb <- txtProgressBar(min=0, max=nrow(x), style=3)

for (j in seq(nrow(x))) {
  for (i in seq(t)) {
    for (p in seq(n)) {
       ker[p] <- dnorm( (x[j,p] - xt[i,p])/h[p])      
    }
    pker[i,1]  <- prod(ker)
    pker[i,2] <- prod(ker)*yt[i]
  } 
  fx[j]  <- mean( pker[,1])/ph
  fxy[j] <- mean( pker[,2] )/ph
  setTxtProgressBar(pb, j)
} 
close(pb)
mx  <- fxy/fx 


#**************************************************************************
#**
#**     Generate graphs
#**
#**************************************************************************

figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

#--------------------------------------------------------#
# Panel (a)
f <- reshape( cbind(x[,1]*(x[,2]^2)),nrow(x1m),nrow(x2m) )

persp(x2m[,1], x1m[1,], f, theta = 45, phi = 35, 
      main = "(a) True Surface",
      ticktype="detailed", nticks = 5,
      xlab = "x2t",
      ylab = "x1t",
      zlab = "\nf(x1t,x2t)",     
      bty="l")


#--------------------------------------------------------#
# Panel (b)
f <- reshape( mx,nrow(x1m),nrow(x2m) )
persp(x2m[,1], x1m[1,], f, theta = 45, phi = 35, 
      main = "(b) Estimated Surface",
      ticktype="detailed", nticks = 5,
      xlab = "x2t",
      ylab = "x1t",
      zlab = "\nf(x1t,x2t)",     
      bty="l")

