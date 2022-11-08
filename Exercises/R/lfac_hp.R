#============================================================================
#
#    Kalman Filter implementation of the Hodrick Prescott Filter 
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()


#
#--------------------------- Helper Functions -------------------------------
#

# Load required functions - inv
source("EMTSUtil.R")

library(matlab) # flipud

#-------------------------------------------------------------------------
#       Univariate Hodrick Prescott Filter in Kalman Filter Form
#       Unrestricted so parameters can be estimated.
#-------------------------------------------------------------------------
neglog <- function(b,y) {
  # Set up filter 
  
  Lam  <- c(1,0)
  Phi  <- matrix(c(1, 1,
                   0, 1), 2,2, byrow=T)  
  R <- b[1]^2
  Q <- matrix(c(0, 0,
                0, b[2]^2), 2,2, byrow=T)
  
  lf <- -mean( kalman(y,Phi,Lam,R,Q) )
  return (lf)
}

#--------------------------------------------------------------------------
# Kalman filter
#--------------------------------------------------------------------------
kalman <- function(y,Phi,Lam,R,Q) {
  # Allocate arrays  
  t <- nrow(y)
  n <- ncol(y)
  k         <- nrow(Q)
  lnl       <- rep(0, t)
  
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- diag(k)*1000
  
  
  mt <- Lam %*% st  
  vt <- Lam %*% pt %*% cbind(Lam) + R
  ut <- y[1,] - mt
  
  
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5* ut %*% inv(vt) %*% cbind(ut)
  Kal <- pt %*% Lam %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt
  
  
  # Main loop over observations
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0     
    pt <- Phi %*% p0 %*% t(Phi) + Q          
    
    # Observation
    mt <- Lam %*% st
    vt <- Lam %*% pt %*% cbind(Lam) + R
    ut <- y[i,] - mt
    
    # Construct log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5* ut %*% inv(vt) %*% cbind(ut)
    
    # Update        	
    Kal <- pt %*% Lam %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt 
  }    
  return(lnl)
}


#-------------------------------------------------------------------------
#       Univariate Hodrick Prescott Filter in Kalman Filter Form
#-------------------------------------------------------------------------
hpsmooth <- function( y ) {
  # Set up filter
  sigu <- 1
  sigz <- sqrt(1/1600)
  
  Lam <- c(1,0)
  Phi <- matrix(c(1, 1,
                  0, 1), 2, 2, byrow=T)
  # R
  R <- sigu
  
  # Q
  Q <- matrix(c(0, 0,
                0, sigz), 2,2, byrow=T)
  
  # Allocate arrays
  t <- nrow(y)
  n <- ncol(y)
  lnl     <- rep(0, t)
  k       <- nrow(Q)
  s10     <- array(0, c(t,k))    # st|t-1
  s11     <- array(0, c(t,k))    # st|t  
  p10     <- array(0, c(k,k,t))
  p11     <- array(0, c(k,k,t))
  ss      <- array(0, c(t,k))
  
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- diag(k)*1000  
  s10[1,] <- st 
  p10[,,1] <- pt
  
  mt <- Lam %*% st  
  vt <- Lam %*% pt %*% cbind(Lam) + R
  ut <- y[1,] - mt
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5* ut %*% inv(vt) %*% cbind(ut)
  
  Kal <- pt %*% Lam %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt
  s11[1,] <- s0
  p11[,,1] <- p0
  
  
  # Main loop over observations  
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0     
    pt <- Phi %*% p0 %*% t(Phi) + Q  
    s10[i,] <- st
    p10[,,i] <- pt
    
    # Observation
    mt <- Lam %*% st  
    vt <- Lam %*% pt %*% cbind(Lam) + R
    ut <- y[i,] - mt
    
    # Construct log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut %*% inv(vt) %*% cbind(ut)
    
    # Update          
    Kal <- pt %*% Lam %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt
    s11[i,] <- s0
    p11[,,i] <- p0
  }
  
  # Now smooth the factor    
  ss <- s11
  for (j in 1:(t-1)) {    
    Jt      <- p11[,,t-j] %*% t(Phi) %*% inv(p10[,,t-j])
    ss[t-j,] <- s11[(t-j),] + Jt %*% t(( t(ss[(t-j+1),]) - t(s10[(t-j+1),]) ) )
  }  
  lfac <- ss   
  lf   <- -mean(lnl)
  return(list(lfac=lfac, lf=lf))
}

#-------------------------------------------------------------------------
#  Univariate Hodrick Prescott Filter in Original Form
#-------------------------------------------------------------------------
HPfilter <- function ( y,lambda ) {
  t <- nrow(y)
  m <- ncol(y)
  
  if (t < m) {
    y <- t(y)
    t <- m
  }
  # Setting up
  td   <- t - 4
  tmp  <- rep(0, t)
  tmpd <- rep(0, td)
  
  # Set up the first and two rows of the filter matrix F
  r1 <- tmp
  r1[1] <- 1.0  
  r1[2] <- -2.0 
  r1[3] <- 1.0
  
  r2 <- tmp
  r2[1] <- -2.0 
  r2[2] <- 5.0  
  r2[3] <- -4.0 
  r2[4] <- 1.0
  
  
  # Create the diagonals for the filter matrix F
  tmp  <- rep(1, t)
  tmpd <- rep(1, td)
  
  D <- cbind(tmpd, -4.0*tmpd, 6.0*tmpd, -4.0*tmpd, tmpd)  
  
  # Construct the filter matrix F
  F <- cbind(r1, r2, spdiags(diag(D), t,td ), flipud(r2), flipud(r1) )
  
  
  HPtrend <- inv(lambda * F +diag(t)) %*% as.matrix(y)
  return(HPtrend)
}

#----------------------------------------------------------------------------
# Creates a sparse matrix or row and col, given the vector v
#----------------------------------------------------------------------------
spdiags <- function(vec, rows, cols) {
  m <- matrix(0, rows, cols) 
  for(j in seq(cols)) {
    i <- j
    for(n in seq(vec)) {     
      m[i, j] <- vec[n]    
      i <- i + 1
    }   
  }
  return(m)
}


#
#--------------------------- Hodrick-Prescott Filter -------------------------
#
lfac_hp <- function( ) {
  # Load quarterly US data for the period 1940:1 to 2000:4 (T = 244)
  RGDP <- read.table("lfac_usgdp.dat")
  
  y <- log(RGDP)-colMeans(log(RGDP))
  y <- y*100
  t <- nrow(y)
  y <- as.matrix(y)
  
  # Estimate unconstrained model
  start <- c(1,1)
  estResults <- optim(start, neglog, y=y, method="BFGS")
  bhat <- estResults$par
  lf1 <- estResults$val
  
  lf1 <- -lf1
  
  cat('\nParameter Estimates\n')
  print( bhat )
  
  # Compute likelihood from the Kalman filter with HP restrictions
  hp.smooth <- hpsmooth(y)
  fac <- hp.smooth$fac
  lf0 <- hp.smooth$lf
  
  
  trend_kf <- fac[,1]
  
  lf0 <- -lf0
  
  cat('\nLog-likelihood function (unconstrained)   = ',lf1) 
  cat('\nLog-likelihood function (constrained)     = ',lf0) 
  
  lr <- -2*t*(lf0 - lf1)
  
  cat('\nLR statistic    = ',lr)
  cat('\np-value         = ',1-pchisq(lr,1))
  
  
  
  # Get trend component from HP filter and plot results
  
  trend_hp <- HPfilter( y,1600 )
  
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  dt <- 1:54
  figure()
  
  matplot( dt,cbind(y[1:54],trend_hp[1:54]), type='l',
           main = '',
           ylab = 'Log real US GDP (times 100)',
           xlab = 't',
           bty = 'l')
}

