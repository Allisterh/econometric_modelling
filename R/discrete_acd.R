#============================================================================
#
#   Estimate the ACD model of Engle-Russell for AMR on Aug01, 2006 
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  inv, trimr, figure, recserar
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  The log-likelihood of the ACD model 
#----------------------------------------------------------------------------
neglog <- function(b,y) {
  b <- abs(b)
  mu <- recserar( cbind((b[1]^2) + (b[2]^2)*trimr( c(0.0,y),0,1)),cbind(mean(y)),cbind(b[3]^2))     
  lf <- -mean ( -log(mu) - y/mu )
  return(lf)
}


#----------------------------------------------------------------------------
#  The log-likelihood of the GARCH model assuming conditional normality 
#----------------------------------------------------------------------------
negloggarch <- function(b,y){
  u  <- sqrt(y)
  s2 <- recserar( cbind((b[1]^2) + (b[2]^2)*trimr( c(0.0, u^2),0,1)),cbind(var(y)),cbind((b[3]^2)))     
  z  <- u/sqrt(s2)                                                     
  
  lf <- -mean (- 0.5*log(2*pi) - 0.5*log(s2) - 0.5*z^2 )
  return(lf)
}


#
#---------------- ACD and GRACH models  -------------------------------------
#
discrete_acd <- function( ) { 
  
  # Load data for AMR on Aug 1, 2006.
  # Order of the variables:
  #   1.  Hour
  #   2.  Minute
  #   3.  Second
  #   4.	Time (time effectively starts at 9.30am and ends at 4.00pm with the time interval being one second)
  #   5.	y (a binary variable giving a one if a trade has occured, and 0 otherwise
  #   6.	N (a counter variable which increases by one when a trade occurs)
  #   7.	u (gives the duration between last trade)
  #   8.	Firm News dummy (AMR)
  #   9.	Macro News dummy
  
  data <- read.table("amr_aug1.dat")
  
  # Trim the first 32 observations and exclude the last observation
  # Ensure postive durations
  data <- trimr(data,32,1)
  
  # Define durations
  y    <- data[,7]
  t    <- length(y)                       
  
  # Descriptive statistics
  cat('\nTotal number of observations             = ',t)
  cat('\nAverage duration time (in seconds)       = ',mean(y))
  cat('\nMinimum duration time (in seconds)       = ',min(y))
  cat('\nMaximum duration time (in seconds)       = ',max(y)) 
  cat('\nUnconditional Hazard rate                = ',1/mean(y))
  cat('\n')
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************** 
  
  
  #--------------------------------------------------------#
  # Panel (a)  
  figure()
  par(mfrow=c(1,2))
  hist(y,21, col="darkblue",
       main = expression(paste("(a) Histogram of duration times y"[t]*"")),
       xlab = "",
       ylab = "",
       bty="l")  
  
  #--------------------------------------------------------#
  # Panel (b)  
  plot(seqa(1,1,51), acf(y, 50, plot=F)$acf,col="darkblue", type="h", lwd=3,
       main = expression(paste("(b) ACF of duration times y"[t]*"")),
       xlab = "",
       ylab="",
       bty = "l")
  
  # Estimate the ACD model
  theta0 <- c(0.1, 0.1, 0.9)
  estResults <- optim(theta0, neglog, y=y, method="BFGS")
  thetac <- estResults$par
  lf <- estResults$val  
  
  lf <- -lf  
  cat('\nLog-likelihood (ACD) model    = ',lf)  
  cat( '\n' )
  
  # Estimate the GARCH model and show the equivalence with the ACD model    
  
  theta0 <- c(0.9670198086627571,
              0.9323466594671925, 
              0.0000091497133080)
  estResults <- optim(theta0, negloggarch, y=y, method="BFGS")
  thetag <- estResults$par
  lf <- estResults$val
  lf <- -lf
  
  cat('\nLog-likelihood (GARCH) model  = ',lf)
  cat( '\n' )
  
  cat('\n Parameter Estimates \n')
  print(cbind(ACD=abs(thetac), GARCH=abs(thetag)))
}

