#=========================================================================
#
#   Effect of the initial condition on mztOLS and mztGLS 
#   Union of rejections included
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# load required functions -  figure, recserar, trimr, inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Detrending function: 
#       cbar = -7 constant 
#       cbar = -13.5 linear trend  
#       cbar = -T for OLS detrending
#----------------------------------------------------------------------------

glsdetrend <- function(y,x,cbar) {
  t <- length(y)
  yc <- cbind(c(y[1], (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1))) )
  xc <- rbind(x[1,], cbind(trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)))    
  b <- lm(yc ~ xc - 1)$coef
  u <- y - x %*% b  
  return(u)
}


#----------------------------------------------------------------------------
#  ar1rvar
#----------------------------------------------------------------------------

ar1rvar <- function(u,k) {
  du <- trimr(u,1,0)-trimr(u,0,1) 
  x  <- cbind(trimr(u,0,1))
  
  b <- lm(trimr(du,k,0) ~ x - 1)$coef
  e <- trimr(du,k,0) - x %*% b 
  s2 <- (t(e) %*% e)/length(e)
  if (k > 0) 
    s2 <- s2/(1-sum(trimr(b,1,0)))^2    
  return (s2)
}

#----------------------------------------------------------------------------
#  M tests
#----------------------------------------------------------------------------
mtests <- function( u,k ){
  
  t  = length(u)
  s2 = ar1rvar( u,k )
  u2 = sum( u[1:(t-1)]^2/t^2 )
  
  mza = (u[t]^2/t-s2)/(2*u2)
  msb = sqrt(u2/s2)
  mzt = msb*mza

  mz = c(mza,  msb, mzt)
  return(mz)
}

# 
#----------------------Role of Union of Rejections---------------------------
#
unit_poweru0UR <- function() {
  # Alternative initial conditions
  u0v <- c(0, 2.5, 5)
  
  # Critical values to compute the power envelope with cv=0 -> size of 0.05           
  cv   <- seq(from=-30, by=1, length.out=31)  
  n    <- length( cv )
  t    <- 200                       
  nreps <- 50000

  # Test parameters
  x    <- cbind(rep(1,t), seq(from=1, by=1, length.out=t) )
  cbar <--13.5 
  cvgls <- -2.867
  cvols <- -3.130 
  tau   <- 1.038

  mztols <- array(0, c( nreps,n,3 ) )
  mztgls <- array(0, c( nreps,n,3 ) )
  rejols <- array(0, c( 3,n ))
  rejgls <- array(0, c( 3,n ))
  rejur  <- array(0, c( 3,n ))
  
  for (i in seq(u0v)) {
    pb <- txtProgressBar(min=0, max=n, style=3)
    for (k in seq(n) ) {
      set.seed(1234, kind="Mersenne-Twister")
      c <- cv[k]      
      for (j in seq(nreps)) {
        u <- cbind(c(u0v[i], rnorm(t)))
        y <- trimr( recserar(cbind(u),cbind(u[1]),cbind(1+c/t)),1,0 )
        
        mztols[j,k,i] <- trimr( mtests(glsdetrend(y,x,-t),0),2,0 )        
        mztgls[j,k,i] <- trimr( mtests(glsdetrend(y,x,cbar),0),2,0 )
        
      }
      setTxtProgressBar(pb, k)      
    }
    close(pb)
    
    # Compute rejection frequencies 
    rejols[i,] <- colMeans(mztols[,,i] < cvols)
    rejgls[i,] <- colMeans(mztgls[,,i] < cvgls)
    rejur[i,]  <- colMeans(mztols[,,i] < (cvols*tau) | mztgls[,,i] < (tau*cvgls))
  }  

  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************

  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,3))

  #--------------------------------------------------------#
  matplot(cv, cbind(rejols[1,],rejgls[1,], rejur[1,]), type="l",
          main = expression(paste("(a) v", ""[0]*"", " = 0.0")),
          ylab = "Power",
          xlab = "Critical Values",
          xlim = c(-30, 0),
          ylim = c(0, 1),
          bty = "l")

  #--------------------------------------------------------#
   matplot(cv, cbind(rejols[2,],rejgls[2,], rejur[2,]), type="l",
          main = expression(paste("(b) v", ""[0]*"", " = 2.5")),
          ylab = "Power",
          xlab = "Critical Values",
          xlim = c(-30, 0),
          ylim = c(0, 1),
          bty = "l")

    
  #--------------------------------------------------------#
   matplot(cv, cbind(rejols[3,],rejgls[3,], rejur[3,]), type="l",
          main = expression(paste("(c) v", ""[0]*"", " = 5.0")),
          ylab = "Power",
          xlab = "Critical Values",
          xlim = c(-30, 0),
          ylim = c(0, 1),
          bty = "l") 
}
