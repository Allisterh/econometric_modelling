#============================================================================
#
#     Cross validation approach to compute optimal bandwidth.
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(11, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - figure
source("EMTSUtil.R")
# load required library - ones, zeros, repmat
library("matlab")

#
#------------------  Cross Validation Bandwidth Selection -------------------
#

# Simulate the model 
t <- 500
ut  <- 0.1*rnorm(t)               #     N(0,0.1^2)                      
xt  <- sort(-2 + 4*runif(t))      #     U(-2,2)  

mx  <- 0.3*exp( -4*(xt + 1)^2 ) + 0.7*exp( -16*(xt - 1)^2 )  
yt  <- mx + ut

# Compute weight function for cross validation (exclude lowest and highest 5# from sample)
xts <- sort(xt)
xl  <- xts[round(0.05*t)]
xu  <- xts[round(0.95*t)]
    
wt  <- zeros(t,1)
for (i in seq(t)) {
  
  if(xt[i]-xl > 0) {
    
    if(xt[i]-xu < 0) {
       wt[i]<- 1      
    }
  }  
}

# Perform a grid search for h over the grid [0.01, 1.0] 
hopt <- seq(0.001, 0.3, 0.001)
n    <- length(hopt)
ace  <- zeros(n,1)

for (i in seq(n)) {
  h   <- hopt[i]
  tmp <- (repmat(xt,1,t) - repmat(xt,t,1))/h  
  fx  <- colMeans( diagrv( t( dnorm(tmp) )/h, zeros(t,1) ) )                  # Leave-one-out   
  fyx <- colMeans( diagrv( t( dnorm(tmp) ) * repmat(yt,1,t)/h, zeros(t,1)) )  # Leave-one-out   
  mx  <- fyx/fx
  
  ace[i] <- sum( ((yt - mx)^2 )*wt )  
}

# Repeat calculations without cross-validation
nce  <- zeros(n,1)

for (i in seq(n)) {
  h   <- hopt[i]
  tmp <- (repmat(xt,1,t) - repmat(xt,t,1))/h
  
  fx  <- colMeans( t( dnorm(tmp) )/h  )                
  fyx <- colMeans( t( dnorm(tmp) ) * repmat(yt,1,t)/h ) 
  mx  <- fyx/fx
  
  nce[i] <- sum( ((yt - mx)^2 )*wt )
}

value <- min(ace)
index <- which.min(ace)
cat('\nBandwidth (with cross-validation)             = ', index )
cat('\nObjective function (with cross-validation)    = ', value)

value <- min(nce)
index <- which.min(nce)
cat('\nBandwidth (without cross-validation)          = ', index )
cat('\nObjective function (without cross-validation) = ', value)

 
#**********************************************************************
#***
#***     Generate graphs
#***
#**********************************************************************

figure()
par(xaxs="i", yaxs="i")
plot(hopt,ace, type="l",
     xlab = expression(h),
     ylab = expression(S(h)),
     ylim = c(0,9),
     bty = "l")
lines(hopt,nce, lty=2, col="blue")


