#============================================================================
#
#  Program to compute the statistical properties of daily equity returns 
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr, figure
source("EMTSUtil.R")

# load required library - bkde
library("KernSmooth")

#
#------------------- Statistical Properties of Equity Returns ---------------
#

garch_statistics <- function( )
{
  # ftse, dow, nk, dummy vars for day of week (tue, wed, thurs, fri, holiday_ftse, holiday_dow, holiday_nk)
  data <- as.matrix(read.table("equity.dat", 
                               col.names=c('FSTSE', 'DOW', 'NIKKEI', 
                                           'tdum',  'wdum',  'thdum',	'fdum',	'hdum_ftse',	'hdum_dow',	'hdum_nk')))
  
  # Choose equity index and compute percentage return
  equity <- data[,1] # FTSE
  y <- 100*(trimr( log( equity ),1,0 ) - trimr( log( equity ),0,1 ))
  
  
  # Compute the moments of the unconditional (empirical) distribution
  y_z <- (y - mean(y))/sd(y)
  cat('\n Skewness\n', mean(y_z^3))
  cat()
  cat('\n Kurtosis\n', mean(y^4))
  
  #*********************************************************************
  #**     Generate graphs
  #*********************************************************************
  
  figure()
  par(mfrow=c(2,2), mar=c(5,5,5,5))
  
  # max lags for autocorrelation of returns and returns squared
  mlag  <- 20
  
  t <- seq(length(y))
  plot(t, y, type="l",
       main="(a) Returns",
       xlab="t",
       ylab=expression(y[t]),
       bty="l")
  
  plot(t, y^2, type="l",
       main="(b) Squared Returns",
       xlab="t",
       ylab=expression(y[t]^2),
       bty="l")
  
  acf(y, lag.max=mlag, plot=TRUE,
      main="(c) ACF Returns",
      xlab="t",
      ylab=expression(acf(y[t])),
      bty="l")
  
  acf(y^2, lag.max=mlag, plot=TRUE,
      main="(d) ACF Squared Returns",
      xlab="t",
      ylab=expression(acf(y[t]^2)),
      bty="l")
  
  
  figure()
  xi <- seq(-6, 6, 0.1)
  f <- bkde(y, range.x=c(-6, 6))
  ft <- dnorm(xi, 0, 1)
  
  plot(xi, ft, type="l", col="blue",
       xlab="y",
       lty=3,
       ylab=expression(f(y)),
       xlim = c(-6, 6),
       ylim=c(0, 0.6),      
       bty="l")
  lines(f, col="red")
  
  detach("package:KernSmooth")
}


