#========================================================================
#
#      Foreign exchange market efficient
#
#========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions------------------------------
#

#load required functions - inv, trimr
source("EMTSUtil.R")

#
#--------------------  Model estimation --------------------------------
#

qmle_sf <- function() {
  # Load data
  data <- as.matrix(read.table("sf.dat"))
  f <- data[,1]
  s <- data[,2]
  
  
  # Spread
  y <- trimr(s,3,0)-trimr(f,0,3)
  
  #********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************
  dates <- seq(from=1979+(4/12), to=2012-(4/12), by=1/12)
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
  plot(dates,trimr(s,3,0), type="l",
       main = "Spot and Forward Rates",
       xlab = "Date",
       ylab = "Rate",
       xlim = c(1978, 2012),
       ylim = c(1.0, 2.8),
       bty = "l")
  lines(dates, trimr(f,0,3), lty=2)
  
  plot(dates,y, type="l",
       main = "Spot and Forward Spread",
       xlab = "Date",
       ylab = "Spread",
       xlim = c(1978, 2012),
       ylim = c(-0.4, 0.3),
       bty = "l")
  
  
  t <- length(y)-3
  x <- cbind(rep(1,t), trimr(y,0,3))
  y <- trimr(y,3,0)
  
  xxi     <- inv(t(x) %*% x)
  bhat    <- lm(y ~ x - 1)$coef
  uhat    <- y-x %*% bhat
  
 
  figure()
  acf(uhat, 6, main='Autocorrelations of uhat')  
  
  figure() 
  acf(x[,2]*uhat, 6, main="Autocorrelations y(t-3)*uhat")
  
  cols <- ncol(x)
  k        <- 1        
  xuhat    <- array(0, c(t,cols))
  
  for (j in seq(cols)) {  
    xuhat[,k] <- x[,j]*uhat
    k<-k+1  
  }
  
  s2      <- t(uhat) %*% uhat/t
  seols   <- sqrt(diag(c(s2) * xxi))
  white   <- xxi %*% (t(xuhat) %*% xuhat) %*% xxi
  sewhite <- sqrt(diag(white))
  
  # Initial estimate of optimal lag length
  P <- floor(4*(t/100)^(2/9)) 
  
  cat('\nInitial estimate of lag length = ', P)
  J0 <- t(xuhat) %*% xuhat
  cat("\nJ0 = \n")
  print(J0)
  
  J1 <- 0
  for (j in seq(P)){
    Gj <- t( trimr(xuhat,j,0) ) %*% trimr(xuhat,0,j)
    J0 <- J0 + Gj + t(Gj)
    J1 <- J1 + 2*j*Gj  
  }
  cat("\nJ1 = \n")  
  print(J1)
  
  # Revised estimate of optimal lag length
  i  <- rbind(1, rep(1, cols-1))
  v0 <- t(i) %*% J0 %*% i
  v1 <- t(i) %*% J1 %*% i
  P  <- floor(1.1447*((v1/v0)^2*t)^(1/3))
  cat('\nRevised estimate of lag length = ', P)
  
  # Compute Newey-West estimate of variance
  JT <- t(xuhat) %*% xuhat 
  for (j in seq(P)) {
    Gj <- t( trimr(xuhat,j,0) ) %*% trimr(xuhat,0,j)
  	JT <- JT + c(1-j/(P+1))*(Gj + t(Gj))  
  }
  
  varNW <- xxi %*% JT %*% xxi
  seNW <- sqrt(diag(varNW))
  
  cat('\n')
  print(cbind(bhat, seols, sewhite, seNW))
  
  # Wald test that both beta1 and beta2 = 0
  wd <- t(bhat) %*% inv(varNW) %*% bhat
  cat('\nWald statistic          = ',wd) 
  cat('\np-value                 = ', 1- pchisq(wd, 2))
}
