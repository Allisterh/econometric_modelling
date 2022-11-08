#============================================================================
#
#   Program to demonstrate the consistency of GMM
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# GMM objective function  
#----------------------------------------------------------------------------
gmmcrit <- function ( alpha,y ) {
  d1 <- y - alpha
  d2 <- y^2 - alpha*(alpha+1)
  d3 <- 1/y - 1/(alpha-1)
  
  d <- cbind(d1, d2, d3)
  g <- colMeans(d)
  w <- t(d) %*% d/nrow(d)  
  Q <- g %*% inv(w) %*% cbind(g)
  return(Q)
}

#
#------------------ Graphical Demonstration of Consistency-------------------
#
gmm_consistency <- function () {
   # Population distribution: T=100000  
  t      <- 100000
  alpha0 <- 10.0
  alpha  <- seq(6,14,0.1)

  y      <- rgamma(t, shape=alpha0,rate=1)   

  Q0     <- rep(0, length(alpha))

  for (k in seq(alpha)) {
    Q0[k] <- gmmcrit( alpha[k],y )   
  }
  
  # Finite sample distributions
  Q <- array(0, c(length(alpha),4))

  t <- 50
  for (k in seq(alpha)) {
    tmp    <- y[1:t]
    Q[k,1] <- gmmcrit( alpha[k],tmp )    
  }

  t <- 100
  for (k in seq(alpha)) {
    tmp  <- y[1:t]
    Q[k,2] <- gmmcrit( alpha[k],tmp )
  }

  t <- 200
  for (k in seq(alpha)) {
    tmp  <- y[1:t]
    Q[k,3] <- gmmcrit( alpha[k],tmp )
  }
  t <- 400
  for (k in seq(alpha)) {
    tmp  <- y[1:t]
    Q[k,4] <- gmmcrit( alpha[k],tmp )
  }
  #**********************************************************************
  #***
  #***     Generate graphs
  #***
  #**********************************************************************

  figure()
  matplot(alpha, cbind(Q0, Q[,1], Q[,2], Q[,3], Q[,4]), type="l",
          xlab = expression(alpha),
          ylab = expression(Q(alpha)),
          bty = "l")
}
