#=========================================================================
#
#   Maximum likelihood estimation of the stationary distribution of the 
#   CIR model of interest rates using Ait Sahalia's (1996) data.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
# 

#-------------------------------------------------------------------------
# Wrapper function to calculate inverse of a given matrix
#-------------------------------------------------------------------------
inv <- function (M) {
  return(solve(M))
}

#-------------------------------------------------------------------------
# Likelihood wrapper function
#-------------------------------------------------------------------------
neglog <- function( p,data ) {
   f <- -sum( lnlt( p,data) )
   return(f)
}

#-------------------------------------------------------------------------
# Likelihood function for stationary distribution of CIR model
#-------------------------------------------------------------------------
lnlt <- function(p,data) {
  v <- abs(p[1])
  w <- abs(p[2])    
  f <- v*log( w ) - lgamma( v ) + (v-1)*log(data) - w*data  
  return(f)
}

# Load any compatibility requirements
source("EMTSUtil.R")

#
#--------------------------- Algorithms -----------------------------------
#

max_stationary <- function() {
  
  # Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
  #   1. year
  #   2. day
  #   3. date stamp
  #   4. interest rates
  load ("eurodollar.RData")
  
  rt <- sort( eurodata[,4]*100 )    
  t  <- length( rt )  
  pstart <- c(5.6556, 0.6763)
  
  results <- optim(pstart, neglog, data = rt, method="BFGS", hessian=TRUE)
  
  phat <- results$par
  hessian <- results$hessian
  
  hessinv <- inv(hessian)
  cat('\nParameter estimates\n')
  cat('\t', phat)
  cat('\nStandard errors based on inverse Hessian\n')
  cat('\t', c(sqrt( hessinv[1,1] ), sqrt( hessinv[2,2] )) )
  
  #********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************

  figure()
  par(xaxs="i", yaxs="i", yaxt='n')
  hist(rt, freq=FALSE, breaks = 50,
       main = paste("Estimated stationary gamma distribution for\n",
                    "Eurodollar interest rates\n",
                    "(1/06/1973-25/02/1995)"),
       ylab = expression(f(r)),
       xlab = expression(r))
  gpdf <- dgamma(rt, shape=phat[1], scale= 1/phat[2])
  lines(rt, gpdf, type="l", lwd=2)
  box(bty = "l")
}
