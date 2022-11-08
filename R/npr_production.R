#============================================================================
#
#    Gaussian nonparametric kernel regression estimator of a
#    production function using a product kernel. 
#
#   Example taken from Pagan and Ullah (1999) Non-Parametric
#    Econometrics.
#    The variables are yt (log of output), x1t and x2t (logs of inputs)
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure
source("EMTSUtil.R")

# load required library - ones, zeros, meshgrid, reshape
library("matlab")

#----------------------------------------------------------------------------
#
#    Procedure to compute the Nadayara-Watson bivariate product kernel 
#    production function using a product kernel. 
#
#    Output is:
#
#     mx is a (n1 x n2) matrix of the conditional mean
#     h1 is the bandwidth for x1
#     h2 is the bandwidth for x2
#
#
#    Inputs are:
#
#     yt is a (t x 1) vector of the dependent variable
#     xt is a (t x 2) vector of the two independent variables
#     x1 is a (n1 x 1) vector of grid points associated with x1
#     x2 is a (n2 x 1) vector of grid points associated with x2
#
#----------------------------------------------------------------------------


kernbiv <- function(yt,xt,x1,x2,s) {
  t <- nrow(xt)
  n <- ncol(xt)
  
  x <- cbind(x1, x2)
  
  tx <- nrow(x)
  nx <- ncol(x)
  
  fac <- -1/(4 + n + 2*s)
  sdxt <- apply(xt, 2, sd)
  h   <- sdxt*(t^fac)
  ph  <- prod(h)
  
  # Estimate density using product kernel
  ker  <- zeros(n,1)
  pker <- zeros(t,n)
  fx   <- zeros(tx,1)
  fxy  <- zeros(tx,1)  
  
  #pb <- txtProgressBar(min=0, max=tx, style=3)
  
  for (j in seq(tx)) {
    for (i in seq(t)) {
      for (p in seq(n)) {
         ker[p] <- dnorm( (x[j,p] - xt[i,p])/h[p] )        
      }
      pker[i,1]  <- prod(ker )
      pker[i,2] <- prod(ker)*yt[i]      
    }
    fx[j]  <- mean( pker[,1] )/ph
    fxy[j] <- mean( pker[,2] )/ph    
   # setTxtProgressBar(pb, j)
  }
  mx  <- fxy/fx 
  #close(pb)
  return(list(mx=mx, h=h))  
}


#
#-----------------  Nonparametric Estimates of Derivatives ------------------
#

# Simulate the nonlinear production function
n <- 2                  # Number of explanatory variables
t <- 200                # Number of observations

xt <- matrix(runif(t*n), nrow=t)        # Uniform random numbers  
ut <- 0.1*rnorm(t)
yt <- -0.2*log( exp(-5*xt[,1]) + 2*exp(-5*xt[,2]) ) + ut

cat('\nStandard deviation of x1t and x2t\n')
sdxt <- apply(xt, 2, sd)
cat(sdxt, '\n') 


# Compute the Nadaraya-Watson bivariate product kernel estimator
xmesh   <- seq(0.0,1,0.02)
mesh <- meshgrid(xmesh)
x1 <- mesh$x
x2 <- mesh$y

# True conditional mean
mx_pop <- -0.2*log( exp(-5*x1) + 2*exp(-5*x2) )    

# Estimated conditional mean
res <- kernbiv(yt,xt,c(x1) ,c(x2),0)  
mx <- res$mx
h <- res$h

# Compute the derivatives of the kernel conditional mean at x1 = x2 = 0.5 
    
# Compute correct bandwidths with s=1
    
s<-1
fac <- -1/(4 + n + 2*s)
sdxt <- apply(xt, 2, sd)
h   <- sdxt*(t^fac)

res <- kernbiv(yt,xt,0.5-h[1],0.5-h[2],1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5,0.5-h[2],1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5+h[1],0.5-h[2],1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5-h[1],0.5,1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5,0.5,1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5+h[1],0.5,1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5-h[1],0.5+h[2],1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5,0.5+h[2],1) 
mx1 <- res$mx
cat(mx1, '\n')

res <- kernbiv(yt,xt,0.5+h[1],0.5+h[2],1)
mx1 <- res$mx
cat(mx1, '\n')


# x1 derivative
res <- kernbiv(yt,xt,0.5-h[1],0.5,1)
mxa <- res$mx

res <- kernbiv(yt,xt,0.5+h[1],0.5,1)
mxb <- res$mx

d1 <- (mxb-mxa)/(2*h[1])                          

# x2 derivative
res <- kernbiv(yt,xt,0.5,0.5-h[2],1)
mxa <- res$mx

res <- kernbiv(yt,xt,0.5,0.5+h[2],1)
mxb <- res$mx

d2 <- (mxb-mxa)/(2*h[2])                          

# Compute the derivatives of the population conditional mean
xm1 <- 0.5      
xm2 <- 0.5

d1_pop <- exp(-5*xm1)/(exp(-5*xm1) + 2*exp(-5*xm2))
d2_pop <- 2*exp(-5*xm2)/(exp(-5*xm1) + 2*exp(-5*xm2))


cat('\nPopulation derivative of x1 (evaluated at x1=0.5)    = ',d1_pop) 
cat('\nNonparametric derivative of x1 (evaluated at x1=0.5) = ',d1) 

cat('\n\nPopulation derivative of x2 (evaluated at x2=0.5)    = ', d2_pop) 
cat('\nNonparametric derivative of x2 (evaluated at x2=0.5) = ', d2) 
    

#**************************************************************************
#**
#**     Generate graphs
#**
#**************************************************************************
figure()

#--------------------------------------------------------#
f <- reshape( mx,nrow(x1),nrow(x2) )
persp(x1[1,], x2[,1], f, theta = -50, phi = 35,
      ticktype="detailed", nticks = 5,
      xlab = "x1t",
      ylab = "x2t",
      zlab = "\nm(x1t,x2t)",     
      bty="l")

