#============================================================================
#
#   Program to find the value of cbar such that the power of the point
#   optimal unit root test against phi = 1+cbar/T has asymptotic power of 0.5
#
#============================================================================
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, inv, recserar
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Search for cbar for a given x
#----------------------------------------------------------------------------
cbar <- function(x,cbar0,reps,t) {
  cbarv <- cbar0 
  powerv <- power(x,cbarv,reps,t)
  
  while (powerv[1] < 0.5) {
    powerv <- rbind(power(x,cbarv[1],reps,t), powerv)
    cbarv  <- rbind(cbarv[1]-0.5,  cbarv)    
  }  
  cbarc <- cbarv[2] - 0.5*(0.5-powerv[2]/(powerv[1]-powerv[2]))  
  return(cbarc)
}

#----------------------------------------------------------------------------
#  Power envelope
#----------------------------------------------------------------------------
power <- function(x,c,reps,t) {
  set.seed(45)      
  
  p0 <- rep(0, reps)
  p1 <- rep(0, reps)
  
  for (rep in seq(reps)) {
    u  <- rnorm(t)
    y0 <- cumsum(u)
    y1 <- recserar(cbind(u),cbind(u[1]),cbind(1+c/t))
    p0[rep] <- Ptest(y0,x,c,0)
    p1[rep] <- Ptest(y1,x,c,0)    
  }
  pow <- mean(p1 < quantile(p0,0.05))  
  return(pow)
}


#----------------------------------------------------------------------------
#  P-test
#----------------------------------------------------------------------------
Ptest <- function(y,x,cbar,k) {
  n  <- length(y)
  uc <- glsdetrend(y,x,cbar) 
  u0 <- glsdetrend(y,x,0)
  s2 <- ar1rvar(uc,k)
  uc <- matrix(c(uc[1],trimr(uc,1,0)-(1+cbar/n)*trimr(uc,0,1)), nrow=n)
  u0 <- matrix(c(u0[1], trimr(u0,1,0)-trimr(u0,0,1)), nrow=n)
  
  pt <- (t(uc) %*% uc-(1+cbar/n) %*% t(u0) %*% u0)/s2
  return(pt)
}

#----------------------------------------------------------------------------
#  Detrending function: 
#       cbar = -7 constant 
#       cbar = -13.5 linear trend  
#       cbar = -T for OLS detrending
#----------------------------------------------------------------------------
glsdetrend <- function( y,x,cbar ) {  
  t <- length(y) 
  yc <- matrix(c(y[1],(trimr(y,1,0)-(1 + cbar/ t) * trimr(y,0,1))), nrow=t)
  xc <- matrix(c(x[1, ],(trimr(x,1,0)-(1 + cbar/ t) * trimr(x,0,1))), nrow=t )
  b <- lm(yc ~ xc - 1)$coef
  u <- y- x %*% b
  return(u)
}

#----------------------------------------------------------------------------
#  ar1rvar
#----------------------------------------------------------------------------
ar1rvar <- function(u,k) {
  du <- cbind(trimr(u,1,0)-trimr(u,0,1) )
  x  <- trimr(u,0,1)
  
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


unit_cbar <- function( ) {
  reps <- 1000
  t    <- 100
  
  cbarc <- cbar(matrix(1, nrow=t, ncol=1),-6,reps,t)
  cbart <- cbar(cbind(rep(1, t), seqa(1,1,t)),-12,reps,t)
  
  tauv    <- seqa(0.15,0.05,15)
  cbartau <- rep(0, length(tauv))
  for (tc in seq(tauv)) {
    TB <- floor(tauv[tc]*t)
    cbartau[tc] <- cbar(cbind(rep(1,t), seqa(1,1,t), cbind(c(rep(0, TB), seqa(1,1,t-TB)))),-13,reps,t)
  }
  
  
  
  cat('\nConstant               = ',cbarc)
  cat('\nLinear trend           = ',cbart)
  cat('\n')
  cat('\nBreak in trend\n')
  print(cbind(tau=tauv, cbar=cbartau))
  
}

