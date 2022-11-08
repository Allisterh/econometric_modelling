#=========================================================================
#
#   Approximate asymptotic power envelope and power curves of ADF test
#   using alternative detrending methods
#
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, recserar, trimr, inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Detrending function: 
#       cbar = -7 constant 
#       cbar = -13.5 linear trend  
#       cbar = -T for OLS detrending
#----------------------------------------------------------------------------

glsdetrend <- function( y,x,cbar ){
  t <- length(y)
  yc <- cbind(c(y[1], (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1))) )
  xc <- rbind(x[1,], cbind(trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)))    
  b <- lm(yc ~ xc - 1)$coef
  u <- y - x %*% b  
  return(u)
}


#-------------------------------------------------------------------------
#  ADF coefficient and t tests: u must be already de-trended
#-------------------------------------------------------------------------
adftests <- function(u,k) {
  du <- trimr(u,1,0) - trimr(u,0,1)
  x  <- trimr(u,0,1)
  
  xxi <- inv(t(x) %*% x) 
  b   <- xxi %*% t(x) %*% trimr(du,k,0) 
  e   <- trimr(du,k,0)- x %*% b   
  s2  <- t(e) %*% e/length(e)
  
  adfstats <- cbind(c(length(u)*b[1], b[1]/sqrt(s2 %*% xxi[1,1])) )  
  return(adfstats)
}

#
#----------------------- Asymptotic Local Power -----------------------------
#
unit_asypower <- function () {
  # Critical values to compute the power envelope with cv=0 -> size of 0.05           
  cv   <- seq(from=-30, by=1, length.out=31)  
  n    <- length( cv )
  t    <- 1000                       
  nreps <- 50000                 
 
  # Detrending parameters
  x    <- cbind(rep(1, t), seq(from=1, by=1, length.out=t))     
  cbar <- -13.5                     
  dfols  <- array(0, c(nreps,n ) )
  dfgls  <- array(0, c( nreps,n ) )
  
  pb <- txtProgressBar(min=0, max=n, style=3)
  for (k in seq(n)) {
    set.seed(1234, kind="Mersenne-Twister")
    c <- cv[k]    
    for (j in seq(nreps)) {
      u  <- rnorm(t)
      yc <- recserar(cbind(u),cbind(u[1]),cbind(1+c/t))
      
      dfols[j,k] <- trimr(adftests(glsdetrend(yc,x,-t),0),1,0)               
      dfgls[j,k] <- trimr(adftests(glsdetrend(yc,x,cbar),0),1,0)      
    }    
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  # Compute rejection frequencies for alternative detrending methods
  
  rejdfols <- colMeans(dfols < (quantile(dfols[, ncol(dfols)], probs=0.05)) )
  rejdfgls <- colMeans(dfgls < (quantile(dfgls[, ncol(dfgls)], probs=0.05)) )
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************

  figure()
  par(xaxs="i", yaxs="i")
  
  #--------------------------------------------------------#
  matplot(cv,cbind(rejdfols, rejdfgls),type="l",
          ylab = 'Power',
          xlab = "Critical Values",
          ylim = c(0.0, 1.0),
          xlim = c(-30, 0),
          bty = "l")  
}
