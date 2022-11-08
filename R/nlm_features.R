#=========================================================================
#
#      Program to demonstrate a few nonlinear features  
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  trimr, figure
source("EMTSUtil.R")

#
#--------------------- Features of Nonlinear Models -------------------------
#
nlm_features <- function() 
{
  # Stationary time series
  #---------------------------------------------------------------------
  t <- 40
  y <- 0.5 + rep(0, t)
  for (i in 3:t)
  {
    y[i] <- 0.9 + 0.7 * y[i-1] - 0.6 * y[i-2]
  }
    
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
  plot(seq(t), y, type="l",
       main='(a)',
       xlab = 't',
       ylab = expression(y[t]),
       bty="l")
  
  plot(trimr(y, 0, 1), trimr(y, 1, 0), type="l",
       main="(b)",
       xlab = expression(y[t-1]),
       ylab = expression(y[t]),
       bty = "l")
  
  # Limit Cycles (taken from Tong (1983, p.85)  
  # ---------------------------------------------------------------------
  t <- 100
  y <- 2 + rep(0,t)
  
  threshold <- 3.05
  for(i in 7:t)
  {
    if(y[i-2] <= threshold)
    {
      y[i] <- 0.8023 + 1.0676*y[i-1] - 0.2099*y[i-2] + 0.1712*y[i-3] - 0.4528*y[i-4] + 0.2237*y[i-5] - 0.0331*y[i-6]
    }else{
      y[i] <- 2.2964 + 1.4246*y[i-1] - 1.0795*y[i-2] - 0.090*y[i-3]
    }    
  }
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
  plot(seq(t), y, type="l",
       main='(a)',
       xlab = 't',
       ylab = expression(y[t]),
       bty="l")
  
  plot(trimr(y, 0, 1), trimr(y, 1, 0), type="l",
       main="(b)",
       xlab = expression(y[t-1]),
       ylab = expression(y[t]),
       bty = "l")
  
  
  # Strange attactor and chaos based on the Kaldor model 
  # (Lorenz (1989, p.130)) 
  #--------------------------------------------------------------------
  
  
  nobs <- 10000
  
  a <- 5.0
  c <- 20.0
  d <- 0.01
  e <- 0.05
  f <- 280.0
  g <- 4.5
  s <- 0.21
  
  alfa <- 20.0
  delta <- 0.05
  eta   <- 0.00001
  
  
  t <- nobs + 5000
  y <- rep(0, nobs)
  k <- rep(0, nobs)
  
  y[1] <- 65.0
  k[1] <- 265.0
  
  for (i in seq(t-1)) {
    invest <- c*2^(-1/((d*y[i] + eta)^2)) + e*y[i] + a*(f/k[i])^g
    saving <- s*y[i]
    y[i+1] <- y[i] + alfa*(invest - saving)
    k[i+1] <- k[i] + invest - delta*k[i]
  }
  y <- trimr(y,5000,0)
  k <- trimr(k,5000,0)
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(2,2))
  plot(y, k,
       main='(a)',
       ylab = expression(k[t]),
       xlab = expression(y[t]),
       bty="l")
  
  plot(trimr(y,0,1),trimr(y,1,0),
       main='(b)',
       ylab = expression(y[t+1]),
       xlab = expression(y[t]),
       bty="l")
  
  plot(trimr(k,0,1),trimr(k,1,0),
       main='(c)',
       ylab = expression(k[t+1]),
       xlab = expression(k[t]),
       bty="l")
  
  # Show dependence on initial conditions
  nobs <- 100
  
  y0 <- rep(0, nobs)
  k0 <- rep(0, nobs)
  y0[1] <- 65.0
  k0[1] <- 265.0
  
  for (i in seq(nobs-1)) {
    invest <- c*2^(-1/((d*y0[i] + eta)^2)) + e*y0[i] + a*(f/k0[i])^g
    saving <- s*y0[i]
    y0[i+1] <- y0[i] + alfa*(invest - saving)
    k0[i+1] <- k0[i] + invest - delta*k0[i]
  }
  
  y1 <- rep(0, nobs)
  k1 <- rep(0, nobs)
  y1[1] <- 65.1
  k1[1] <- 265.1
  
  for (i in seq(nobs-1)) {
    invest <- c*2^(-1/((d*y1[i] + eta)^2)) + e*y1[i] + a*(f/k1[i])^g
    saving <- s*y1[i]
    y1[i+1] <- y1[i] + alfa*(invest - saving)
    k1[i+1] <- k1[i] + invest - delta*k1[i]
  }
  
  
  plot(seq(nobs),y0, type="l",
       lty=1,
       main = "d",
       xlab = "Time",
       ylab = expression(y[t]),
       bty = "l")
  lines(seq(nobs),y1, lty=2)
  
  # Multiple equilibria (taken from Ericsson, Hendry and Prestwich 
  # Scandinavian Journal of Economics (1998, p.296)    
  #--------------------------------------------------------------------
  nsim <- 50000
  rm   <- rep(0, nsim) + 0.1
  u    <- rep(0, nsim)
  v    <- 0.005*rnorm( nsim)
  
  for (i in 5:nsim) {
    rm[i] <- rm[i-1] + 0.45*(rm[i-1] - rm[i-2]) - 0.10*(rm[i-2] - 2*rm[i-3] + rm[i-4])- 2.54*(rm[i-1] - 0.1)*rm[i-1]^2 + v[i]
  }
  drm <- trimr( rm,1,0 ) - trimr( rm,0,1 )
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i")
  hist(rm, breaks=51,
       xlab = "Midpoint",
       ylab = "Frequency",
       bty="l")
}


