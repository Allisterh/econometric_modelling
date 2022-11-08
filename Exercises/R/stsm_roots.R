#=========================================================================
#
#   Program to identify properties of US macro variables
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - trimr
source("EMTSUtil.R")
#
#------  Stationarity Propperties of US Macroeconomic variables -------------
#
stsm_roots <- function() {
  
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
  
  
  # Construct variables for use in VAR: from 1960:1 to 1998:12    
  # interest rate and the annual # growth rates in money, price and output
  yvar <- cbind(r,  lm,  lp,  lo)             
  
  # Compute the roots of the AR(4) model for the interest rate
  y     <- r
  X     <- cbind( rep(1, length(y)-4), trimr(y,3,1), trimr(y,2,2), trimr(y,1,3), trimr(y,0,4) )    
  theta <- lm(trimr(y,4,0) ~ X - 1)$coef
  
  c  <- c(1, -theta[2], -theta[3], -theta[4], -theta[5])
  rt <- polyroot( c )
  
  od <- options(digits=6)
  cat('\nAR(4) model for interest rates')
  cat('\n------------------------------\n')
  print(cbind(Roots=rt, Abs.Roots= abs(rt)))
  cat('\n  ')
  
  # Compute the roots of the AR(2) model for real output growth
  y     <- 100*(trimr(lo,12,0) - trimr(lo,0,12))
  X     <- cbind(rep(1, length(y)-2), trimr(y,1,1), trimr(y,0,2))
  theta <- lm(trimr(y,2,0) ~ X - 1)$coef
  
  c <- c(1, -theta[2], -theta[3])
  rt <- polyroot( c )
  
  cat('\nAR(2) model for real output growth')
  cat('\n------------------------------\n')
  print(cbind(Roots=rt, Abs.Roots= abs(rt)))
  cat('\n  ')
  
  # Compute the roots of the AR(2) model for log output 
  y     <- log( ytdata[,6] )
  X     <- cbind(rep(1, length(y)-2), trimr(y,1,1), trimr(y,0,2))
  theta <- lm(trimr(y,2,0) ~ X - 1)$coef
  
  c <- c(1, -theta[2], -theta[3])
  rt <- polyroot( c )
  
  cat('\nAR(2) model for log of real output')
  cat('\n------------------------------\n')
  print(cbind(Roots=rt, Abs.Roots= abs(rt)))
  cat('\n  ')
  
  options(od) # reset the options
}

