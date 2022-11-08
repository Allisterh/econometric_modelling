#============================================================================
#
#   Estimate a MIDAS regression model for the hedge funds 
#
#============================================================================
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, figure, inv, reshapeg
source("EMTSUtil.R")


#----------------------------------------------------------------------------
#  Log-likelihood function
#----------------------------------------------------------------------------
neglog <- function(b,rh,rm,rhd2,d_gfc,tm,n,lag) {
  theta1 <- b[5]
  theta2 <- b[6]    
  
  # Precompute weights
  maxlag <- n*lag
  w      <- getweights(theta1,theta2,maxlag)
  
  # Now need 3 counters, k,i,j to compute midas variance
  # k is a simple placeholder for the calculated monthly variance
  # i loops over the daily data and is incremented by the number of days in month
  # j loops from 1 to maxlag over the daily data to compute the required sum  
  
  vm <- rep(0, tm)
  
  # If maxlag is greater than number of days in a month then pad the returns
  if (maxlag > n ) {
    rhd2 <- cbind(c(rep(0, maxlag-n),  rhd2))
  }  
  k <- 1

  for (j in seq(from=(maxlag+1), to=(length(rhd2)+1), by=n)) {
    for (i in seq(maxlag)) {
      vm[k] <- vm[k] + rhd2[j-i]*w[i]
    }    
    k <- k+1    
  }
      
  # Trim extraneous start up months and lag the result  
  if (lag > 1) {
    vtmp <- trimr(vm,lag-1,1)
  } else {
    vtmp <- trimr(vm,0,1)
  }    
  v  <- trimr(rh,lag,0) - b[1] - b[2]*trimr(rm,lag,0) - b[3]*vtmp - b[4]*trimr(d_gfc,lag,0)   
  s2 <- t(v) %*% v/length(v)                                              
  z  <- v/sqrt(s2)
  
  f  <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*z^2   
  lf <- -mean( f)
  return(lf)
}

#----------------------------------------------------------------------------
#  Return weights normalised to sum to unity
#----------------------------------------------------------------------------
getweights <- function(theta1,theta2,maxlag) {
  i <- c(0, seq(maxlag-1))  
  w <- exp( theta1*i/1.0e2 + theta2*i^2/1.0e4  )
  w <- w/sum(w)
  return(w)
}

#
# ------------------------ GARCH MIDAS Model --------------------------------
# 
garch_midas <- function( ) {
  # Load daily data: March 31st 2003 to Aug. 4th 2011 
  data <- read.table("daily_hedge.dat", header=T)
  
  # Compute daily returns
  rhd <- trimr(log(data[,'p_hedge']),1,0) - trimr(log(data[,'p_hedge']),0,1)                  
  rmd <- trimr(log(data[,'p_market']),1,0) - trimr(log(data[,'p_market']),0,1)                 
  
  # Compute monthly returns from daily returns 
  # assuming a month has 22 days (starting with April 2003)
  rh <- rowSums( reshapeg( rhd,length(rhd)/22,22))
  rm <- rowSums( reshapeg( rmd,length(rmd)/22,22))
  
  #  Lag length to compute the MIDAS 
  # lag=1 -> conditional variance based on daily returns in the previous month
  # lag=2 -> conidtional variance based on daily returns in the previous two months         
  lag <- 18           
  tm <- length(rh)           #  Number of months in the monthly data set
  td <- length(rhd)          #  Number of days in the daily data set
  n <- round(td/tm)          #  Number of days in a month
  
  # Squared daily returns scaled by 22 to convert to monthly returns
  # already lagged given the way that the daily data are defined 
  rhd2 <- 22*rhd^2                  
  
  
  # Construct global financial crisis dummy variable
  d_gfc          <- rep(0, tm)
  k <- which.min(rh)
  
  d_gfc[k] <- 1
  
  
  # Estimate the MIDAS model
  start <- 0.1*rep(1, 6)
  estResults <- optim(start, neglog, rh=rh, rm=rm, rhd2=rhd2, d_gfc=d_gfc, tm=tm, n=n, lag=lag, method="BFGS", hessian=T)
  theta <- estResults$par
  lf1 <- estResults$val
  hess <- estResults$hessian
  
  lf1 <- -lf1
  vc  <- (1/(tm-lag))*inv(hess)
  
  # Plot estimated weight function 
  w <- getweights(theta[5],theta[6],n*lag)  
  figure()
  plot(w, type="l", 
       main = "",
       ylab = "MIDAS weights",
       xlab = "i",
       bty = "l")
  
  # Wald test of theta(5)=theta(6)=0
  r <- matrix(c(0,  0,  0,  0,  1,  0,
                0,  0,  0,  0,  0,  1), nrow=2, ncol=6, byrow=T)
  q <- rbind(0,0)
  w <- t( (r %*% theta - q) ) %*% inv(r %*% vc %*% t(r)) %*% (r %*%theta - q)
  
  cat('\nWald statistic     = ',w)
  cat('\np-value            = ',1-pchisq(w,2))
  
}

