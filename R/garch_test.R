#=========================================================================
#
#      Program to test for and estimate an ARCH(p) model
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr, lag.matrix
source("EMTSUtil.R")

#-------------------------------------------------------------------------
#  Test for ARCH(maxlag)
#-------------------------------------------------------------------------
testarch <- function(data,maxlag) {  
  # Create lags of y
  tmp <- lag.matrix(cbind(data),1:maxlag)
  tmp <- tmp[-(1:maxlag),]
    
  y <- trimr(data,maxlag,0)^2
  x <- cbind(rep(1,nrow(tmp)), tmp^2)
  b <- lm(y ~ x - 1)$coef
  e <- y - x %*% b
  
  t   <- length(y)    
  yc  <- y-mean(y)
  rss <- t(e) %*% e
  tss <- t(yc) %*% yc
  r2  <- 1- rss/tss
  lm  <- t*r2
  pv  <- 1-pchisq(lm,maxlag)
  
  return(list(lm=lm, pv=pv))  
}


#-------------------------------------------------------------------------
# Likelihood function for an ARCH(p) model
#-------------------------------------------------------------------------
neglog <- function( b,data,maxlag ) 
{
  b <- abs(b)
  
  # Compute u and lags of u
  u   <- data 
  tmp <- lag.matrix(cbind(data),1:maxlag)
  tmp <- tmp[-(1:maxlag),]
  
  # Compute h
  h <- cbind(rep(1, nrow(tmp)), tmp^2 ) %*% b
  h <- cbind(c(sd(u)^2*rep(1, maxlag), h))
  
  # Likelihood function
  f  <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*(u/sqrt( h ))^2
  lf <- -mean( f )
  return(lf)  
}

#
#---------------------------- ARCH Tests ------------------------------------
#

garch_test <- function() 
{
  # ftse, dow, nk, dummy vars for day of week (tue, wed, thurs, fri, holiday_ftse, holiday_dow, holiday_nk)
  data <- as.matrix(read.table("equity.dat", 
                               col.names=c('FSTSE', 'DOW', 'NIKKEI', 
                                           'tdum',  'wdum',  'thdum',  'fdum',	'hdum_ftse',	'hdum_dow',	'hdum_nk')))
  
  # Choose equity index and compute percentage return
  equity <- data[,1] # FTSE
  y <- 100*(trimr( log( equity ),1,0 ) - trimr( log( equity ),0,1 ))
  
  
  
  y <- y-mean(y)
  
  # Set ARCH order
  p <- 2
  
  # Test for ARCH(p)
  results <- testarch(y,p)
  lm <- results$lm
  pv <- results$pv
  
  # Estimate ARCH(p)
  start      <- runif(p+1)
  estResults <- optim(start, neglog, data=y, maxlag=p, method="BFGS")
  theta <- estResults$par
  lf <- estResults$val
  
  theta      <- abs(theta)    
  
  cat('\nARCH Estimation and Testing')
  cat('\n---------------------------')
  cat('\nARCH order          = ',p)
  cat('\nalpha_0             = ',theta[1])
  cat('\nalpha_1             = ',theta[2])
  cat('\nLikelihood function = ',-lf)
  cat('\n ')
  cat('\nLM test             = ',lm) 
  cat('\np-value             = ',pv)
}


