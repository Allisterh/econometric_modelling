# ===========================================================================
#
#   Program to estimate heteroskedastic models of the business cycle 
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#

# load utility functions - trimr 
source("EMTSUtil.R")

 
#----------------------------------------------------------------------------
# Unrestricted likelihood function 
#----------------------------------------------------------------------------
neglog <- function(p,y,x,w) {
  m  <- p[1] + p[2]*x[,1] + p[3]*x[,2]
  s2 <- exp(p[4] + p[5]*w[,1] + p[6]*w[,2])           
  lt <- -0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m)^2)/s2        
  
  lf <- - mean( lt )
  return(lf)
}

#----------------------------------------------------------------------------
# Restricted likelihood function 
#----------------------------------------------------------------------------
neglog0 <- function(p,y,x) {
  m  <- p[1] + p[2]*x[,1] + p[3]*x[,2]
  s2 <- exp(p[4])           
  lt <- -0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m)^2)/s2        
  
  lf <- - mean( lt )
  return(lf)
}


#
#---------------- heteroskedastic of the business cycle----------------------
#
hetero_g7 <- function() {
  load('G7Data.RData')
  # Set country
  data <- canada
  # france
  # germany
  # italy
  # japan
  # uk
  # us
  
  # Get the data 
  y  <- data[,1]                # Growth rate of GDP
  x  <- data[,2]                 # Spread
  
  
  # Mean variables: lagged gdp growth and lagged interest spread
  x <- trimr(cbind(y, x), 0,1)       
  
  # Variance variables: lagged gdp growth and lagged interest spread
  w <- x
  
  y <- trimr(y,1,0)
  t <- length(y)
  
  # Estimate unrestriced model
  start <- 0.1*rep(1, 6)   
  estResults <- optim(start, neglog, y=y, x=x,w=w, method="BFGS")
  bhat1 <- estResults$par
  lf1 <- estResults$val
  lf1 <- -lf1
  
  # Estimate restriced model
  start <- 0.1*rep(1, 4)
  estResults <- optim(start, neglog0, y=y, x=x, method="BFGS")
  bhat0 <- estResults$par
  lf0 <- estResults$val
  lf0 <- -lf0
  
  # LR test
  lr  <- -2*t*(lf0 - lf1)
  dof <- length(bhat1) - length(bhat0)
  
  cat('\n')
  cat('\nLog-likelihood function (unrestricted) = ',lf1)
  cat('\nLog-likelihood function (restricted)   = ',lf0)
  cat('\nLR statistic                           = ',lr)
  cat('\np-value                                = ',1-pchisq(lr,dof))      
}

# call main function
hetero_g7()

