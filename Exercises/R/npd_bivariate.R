#============================================================================
#
#    Bivariate kernel example
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

# Generate data
n  <- 2                      # Dimension                          
t  <- 20000                 # Sample size
yt <- matrix(rnorm(t*n), nrow=t, ncol=n)


# Set up a grid between -4 and 4
ymesh <- seq(-4, 4, 0.2)

# Need to work out all possible points at which to compute bivariate
# density on this mesh and express results as vectors
mesh <- meshgrid(ymesh)

y1m <- mesh$x
y2m <- mesh$y

y <- cbind(c(y1m), c(y2m) )

# Estimate density using product kernel
fac <- 1/(4.0+n)
sdyt <- apply(yt, 2, sd)
h   <- sdyt/(t^fac)
ph  <- prod(h)

ker  <- zeros(n,1)
pdf  <- zeros(nrow(y),1)
pker <- zeros(t,1)

pb <- txtProgressBar(min=0, max=nrow(y), style=3)
for (j in seq(nrow(y))) {
  for (i in seq(t)) {
    for (p in seq(n)) {      
      ker[p] <- dnorm( (y[j,p] - yt[i,p])/h[p] )      
    }
    pker[i] <- prod(ker)
  }
  pdf[j] <- mean( pker )/ph  
  setTxtProgressBar(pb, j)
}  
close(pb)

#**************************************************************************
#**
#**     Generate graph
#**
#**************************************************************************

figure()
f <- reshape(pdf, nrow(y1m), nrow(y2m))
persp(y2m[,1], y1m[1,], f, theta = 45, phi = 35, 
      ticktype="detailed", nticks = 5,
      xlab = "y2t",
      ylab = "y1t",
      zlab = "\nf(y1t,y2t)",
      zlim = c(0,0.2),
      bty="l")
