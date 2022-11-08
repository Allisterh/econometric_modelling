#============================================================================
#
#   Program to test a liquidity trap for the United States
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv, numhess
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#   Negative log-likelihood function
#----------------------------------------------------------------------------
neglog <- function(b,y,x1,x2) {
  t   <- length(y)
  m   <- b[1]*x1 +  b[2]/(x2 -  b[3])
  u   <- y - m
  s2  <- t(u) %*% u/t
  lnl <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m)^2)/s2
  
  lf <- -mean(lnl)
  return(lf)
}
# 
#-------------------- Testing for a Liquidity Trap --------------------------
#
nls_liquiditytrap <- function( ) {
  # Load data for the United States: January 1959 to December 2011   
  # Variables are
  #       m2
  #       gdp (real)
  #       cpi
  #       interest
  load('us_liquiditytrap.Rdata')
  
  cpi <- data[,1]
  gdp <- data[,2]
  interest <- data[,3]
  m2 <- data[,4]
  
  # Construct variables 
  y  <- log(m2/cpi)        
  x1 <- log(gdp)           
  x2 <- interest/100       
  
  t <- length(y)
  
  # Estimate linear model by OLS         
  x    <- cbind(rep(1, t),   x1,   x2)
  bols <- lm(y ~ x - 1)$coef
  u    <- y - x %*% bols
  s2 <- t(u) %*% u/t
  
  # LM test (2-step) 
  # Intercept is used in the construction of the test 
  # so the restricted model has an intercept. **/
  
  # Stage 1 regression
  x <- cbind(rep(1, t), x1, 1/x2 )
  b1 <- lm(y ~ x - 1)$coef
  u <- y - x %*% b1      
  
  # Stage 2 regression
  z  <- cbind(x,  1/(x2^2))
  b2 <- lm(u ~ z -1)$coef
  v  <- u - z %*% b2                      
  r2 <- 1 - (t(v) %*% v)/( t(u) %*% u)                   
  lm <- t*r2              
  
  cat('\n')
  cat('\nLM test (2-step) with intercept  = ',lm)
  cat('\np-value                          = ',1-pchisq(lm,1)) 
  
  # Estimate model by MLE
  start <- bols[c(2, 1, 3)]    # Note change of order of OLS estimates
  estResults <- optim(start, neglog, y=y, x1=x1, x2=x2, method="BFGS", hessian=T)
  bhat <- estResults$par
  hess <- estResults$hess
  
  vc <- (1/t)*inv(hess)
  
  # Wald test
  wd <- (bhat[3] - 0)^2/vc[3,3]
  
  cat('\n')
  cat('\nWald test        = ', wd)
  cat('\np-value          = ', 1-pchisq(wd,1))  
}


