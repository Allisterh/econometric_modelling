#============================================================================
#
#   Program to identify properties of US macro variables
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure, trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Computes the acf and pacf for the given variables
#----------------------------------------------------------------------------

compute_vars <- function(y, lags, doPlot=F) {
  
  acf1  <<- acf( y[,1],lags, plot=doPlot )$acf
  pacf1 <<- pacf( y[,1],lags, plot=doPlot )$acf
  acf2  <<- acf( y[,2],lags, , plot=doPlot )$acf
  pacf2 <<- pacf( y[,2],lags, plot=doPlot )$acf
  acf3  <<- acf( y[,3],lags, , plot=doPlot )$acf
  pacf3 <<- pacf( y[,3],lags, , plot=doPlot )$acf
  acf4  <<- acf( y[,4],lags, , plot=doPlot )$acf
  pacf4 <<- pacf( y[,4],lags, , plot=doPlot )$acf
}

#
#------  Time Series Properties of US Macrieconomic Variables --------------
#
stsm_properties <- function() {
    # Read the data: quarterly US data from Jan-1959 to Dec-1998
  ytdata <- as.matrix(read.table("sims_data.dat"))
 
  # Define varaibles
  r    <- ytdata[,1]        
  #lex  <- log( ytdata[,2] )
  #lcp  <- log( ytdata[,3] )
  lm   <- log( ytdata[,4] )
  lp   <- log( ytdata[,5] )
  lo   <- log( ytdata[,6] )
  #sdum <- ytdata[,7:17]

  # Number of lags used to compute the ACF and PACF
  lags <- 6 
    
  # Compute the ACF and PACF on levels of variables
  y     <- cbind(r, lm, lp, lo)
  compute_vars(y, lags) 
     
  cat('\n ACF of series in levels ')
  cat('\n------------------------------------------------------------------\n')
  rnames = rep("Lag", lags)
  cnames = c("", "y1", "y2", "y3", "y4")
  tblResults <- matrix(c(seq(lags), acf1[-1], acf2[-1], acf3[-1], acf4[-1]),
                       ncol=length(cnames), dimnames = list(rnames, cnames))
  print(tblResults)
    
  cat('\n\n PACF of series in levels')
  cat('\n------------------------------------------------------------------\n')
  tblResults <- matrix(c(seq(lags), pacf1, pacf2, pacf3, pacf4),
                       ncol=length(cnames), dimnames = list(rnames, cnames))
  print(tblResults)

  # Compute the ACF and PACF on first difference of variables
  dy   <- trimr(y,1,0) - trimr(y,0,1)         
  compute_vars(dy, lags) 
  
  cat('\n ACF of first difference of the series ')
  cat('\n------------------------------------------------------------------\n')
  tblResults <- matrix(c(seq(lags), acf1[-1], acf2[-1], acf3[-1], acf4[-1]),
                       ncol=length(cnames), dimnames = list(rnames, cnames))
  print(tblResults)
    
  cat('\n PACF of first difference of the series')
  cat('\n------------------------------------------------------------------\n')
  tblResults <- matrix(c(seq(lags), pacf1, pacf2, pacf3, pacf4),
                       ncol=length(cnames), dimnames = list(rnames, cnames))
  print(tblResults)
    
  # Compute the ACF and PACF on second difference of variables
  ddy <- trimr(y,12,0) - trimr(y,0,12)       
  compute_vars(ddy, lags)

  cat('\n ACF of 12th difference of the series ')
  cat('\n------------------------------------------------------------------\n')
  tblResults <- matrix(c(seq(lags), acf1[-1], acf2[-1], acf3[-1], acf4[-1]),
                       ncol=length(cnames), dimnames = list(rnames, cnames))
  print(tblResults)
    
  cat('\n PACF of 12th difference of the series ')
  cat('\n------------------------------------------------------------------\n')
  tblResults <- matrix(c(seq(lags), pacf1, pacf2, pacf3, pacf4),
                       ncol=length(cnames), dimnames = list(rnames, cnames))
  print(tblResults)

  #**************************************************************************
  #   
  # Generate graphs
  #
  #**************************************************************************
  figure()
  par(yaxt='n')
  t <- t <- seq(1959+1/12, 1999, 1/12)
  
  matplot(t, y, type="l",
          main = expression(y[t]),
          ylab="",
          xlab="t",
          lty=c(2,5), col=2:5)
  legend("topleft", 
         legend=c("Interest rate", "log money stock", "log of prices", "log of output"),
         lty=c(2,5), col=2:5)
    
}


