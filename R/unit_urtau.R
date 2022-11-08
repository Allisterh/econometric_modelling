#============================================================================
#
#   Compute tau for Union of Rejections MZt unit root test
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

set.seed(123) # mt19937

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr, lag.matrix
source("EMTSUtil.R")

#-------------------------------------------------------------------------
#  Detrending function: 
#       cbar = -7 constant 
#       cbar = -13.5 linear trend  
#       cbar = -T for OLS detrending
#-------------------------------------------------------------------------
glsdetrend <- function(y,x,cbar)  {  
  t <- length(y) 
  yc <- matrix(c(y[1],(trimr(y,1,0)-(1 + cbar/ t) * trimr(y,0,1))), nrow=t)
  xc <- matrix(c(x[1, ],(trimr(x,1,0)-(1 + cbar/ t) * trimr(x,0,1))), nrow=t )
  b <- lm(yc ~ xc - 1)$coef
  u <- y- x %*% b 
  
  return(u)
} 

#---------------------------------------------------------------------------- 
# M tests 
#---------------------------------------------------------------------------- 
mtests <- function(u,k)  {  
  s2 <- ar1rvar(u,k) 
  n <- length(u) 
  tmp <- sum(u[1:(n-1)]^2)   
  
  u2 <- tmp/n^2 
  mza <- (u[n]^2/n-s2)/ (2*u2)
  msb <- sqrt(u2/ (s2))
  mzt <- msb %*% mza 
  
  tests <- cbind(mza,msb,mzt) 
  
  return(tests)
}
#---------------------------------------------------------------------------- 
# Autoregressive long run variance estimator 
#---------------------------------------------------------------------------- 
ar1rvar <- function(u,k)  {
  
  du <- cbind(trimr(u,1,0)-trimr(u,0,1) )
  x <- trimr(u,0,1) 
  
  if (k > 0)  { 
    x <- cbind(x,lag.matrix(du,seqa(1,1,k))) 
    x <- cbind(x[-(1:k), ])
  }
  x <- cbind(x)  
  b <- lm(trimr(du,k,0) ~ x - 1)$coef
  
  e <- trimr(du,k,0)- x %*% b 
  s2 <- t(e) %*% e/ length(e)  
  if (k > 0)  { 
    s2 <- s2/ (1-sum(trimr(b,1,0)))^2
  }   
  return(s2)
}


# 
#------------------------------ Union of rejection tests --------------------
#
unit_urtau <- function( ) {
  t    <- 1000 
  reps <- 100000
  
  # Constant case
  x     <- cbind(rep(1, t))
  cbar  <- -7 
  cvols <- -2.467 
  cvgls <- -1.944
  
  # Obtain replications of MZt-OLS and MZt-GLS tests
  MZols <- rep(0, reps)
  MZgls <- rep(0, reps)
  
  pb <- txtProgressBar(min=0, max=reps, style=3)
  for (rep in seq(reps)) { 
    y          <- cumsum(rnorm(t))
    MZols[rep] <- trimr(t(mtests(glsdetrend(y,x,-t),0)),2,0)
    MZgls[rep] <- trimr(t(mtests(glsdetrend(y,x,cbar),0)),2,0)
    setTxtProgressBar(pb, rep)
  }
  close(pb)
  
  
  
  # Search for tau such that UR test has size of 0.05
  rejtau <- 1 
  tau    <- 1
  
  # Second decimal place:
  while (rejtau >= 0.05) {
    rejtau <- mean(MZols < (tau*cvols) | MZgls < (tau*cvgls))
    tau    <- tau + 0.01  
  }
  
  rejtau <- 1 
  tau <- tau - 0.02
  
  # Third decimal place:
  while (rejtau >= 0.05) {
    rejtau <- mean(MZols < (tau*cvols) | MZgls < (tau*cvgls))
    tau    <- tau + 0.001
  }
  cat('\ntau     = ',tau-0.001)  
}


