#=========================================================================
#
#   Maximum likelihood estimation of the transitional distribution of the 
#   CKLS model of interest rates using Ait Sahalia's (1996) data.
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions ----------------------------
# 
#load required functions - inv, figure
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Quasi-likelihood function for transitional distribution of CKLS model  
#-------------------------------------------------------------------------
ckls <- function(p,data) {
   t <- length(data)-1
   
   alpha <- abs(p[1])
   mu    <- abs(p[2])
   sigma <- abs(p[3])
   gamma <- abs(p[4])
    
   rnow  <- data[-1]
   rlag  <- data[-length(data)]
   dt    <- 1/250
   drift <- alpha*(mu-rlag)*dt
   stdev <- sigma*rlag^(gamma)*sqrt(dt)

   ut  <-  (rnow - rlag - drift)/stdev

   tmp <- -0.5*ut^2 - 0.5*log(stdev^2) - 0.5*log(2*pi)
   lf  <- -sum(tmp)
   return(lf)
}

#
#--------------------  QMLE CKLS Model -----------------------------
#

qmle_ckls <- function() {
  # Load data (5505x4 array called eurodata, 1 Jun 1973 - 25 Feb 1995)
  #   1. year
  #   2. day
  #   3. date stamp
  #   4. interest rates
  load('eurod.RData')
 
  dt <- 1/250
  rt <- eurodata[,4]

  # Starting values
  x  <- rt[1:(length(rt)-1)]
  dx <- diff(rt)           
  dx <- dx/x^0.5
  regressors <- cbind(dt/x^0.5, dt*x^0.5)
  drift <- lm(dx ~ regressors - 1)$coef
  res   <- regressors %*% drift - dx
  alpha <- -drift[2]
  mu    <- -drift[1]/drift[2]
  sigma <- sqrt(var(res)/dt)
    
  p0 <- c(abs(alpha), abs(mu), abs(sigma), 0.5)
    
  # Estimation based on scaled Bessel function
  estResults <- optim(p0, ckls, data=rt, method="BFGS", hessian=T)
  phat <- estResults$par
  hessian <- estResults$hessian

  
  hessinv <- inv(hessian)
  cat('\nParameter estimates\n')
  cat(phat, '\n')
  cat('\nStandard errors based on inverse Hessian\n')
  cat(sqrt( diag(hessinv) ), '\n')
  
  #********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************  
  figure()
  par(xaxs="i", yaxs="i", yaxt='n', mar=c(5, 5, 4, 1))
  rnow  <- rt[-1]
  rlag  <- rt[-length(rt)]

  tmp0 <- seq(0.03, 0.25, 0.01)
  tmp1 <- phat[3]*tmp0^phat[4]

  plot(rlag, (rnow-0.029)^2, pch=19, col="darkgray",
       xlab = expression(r[t-1]),
       ylab = expression(r[t]^2))
  lines(tmp0,tmp1^2, lwd=2)  
}


