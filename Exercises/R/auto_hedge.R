#=========================================================================
#
#       Program to compute the MIC model using the hedge data
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
#

#load required functions - inv, recserar
source("EMTSUtil.R")
library("matlab")

# -----------------------------------------------------------------------
#   Negative log-likelihood function
# -----------------------------------------------------------------------
neglog <- function(b,y,m) {
  beta0 <- b[1]                             
  beta1 <- b[2]
  rho1  <- b[3]
  u     <- y - beta0 - beta1*m
  v     <- trimr(u,1,0) - rho1*trimr(u,0,1)
  sig2  <- mean(v^2)
  lnl   <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2
  lf    <- -mean(lnl)
  return(lf)  
}

#
#---------------------------  Hedge Fund model -----------------------------------
#

auto_hedge <- function( ) {
  # Load daily hedge fund data (1 April 2003 - 28 May 2010)
  load("hedge.RData")
  
  # 1 - 7 daiy hedge fund returns (not excess returns)
  # 8.  Market excess return (this IS adjusted for risk free)
  # 9.  Risk free rate (expressed on a daily basis)
  # 10. Market excess rets with gaps filled in (this IS adjusted for risk free)
  # 11. Risk free rate with gaps filled in (expressed on a daily basis)
  # 12. Dow daily returns in percentage (not excess)
  # 13. NASDAQ daily returns in percentage (not excess)
  # 14. SP500 daily returns in percentage (not excess)
  # 15. Tuesday dummy
  # 16. Wednesday dummy
  # 17. Thursday dummy
  # 18. Friday dummy
  # 19. Holiday dummy

  # Hedge fund excess ret        
  hedge     <- hedgedata[,1:7] - repmat(hedgedata[,11],1,7) 
  # Market excess returns
  m_dow     <- hedgedata[,12] - hedgedata[,11]                          
  m_nasdaq  <- hedgedata[,13] - hedgedata[,11]                       
  m_sp500   <- hedgedata[,14] - hedgedata[,11] 
  # Dummy variables 
  d_season  <- hedgedata[,15:18]                           
  d_hol     <- hedgedata[,19]

  # Choose a hedge fund
  y <- hedge[,3]
  m <- m_sp500
  t <- length(y)

  # LM test applied to the CAPM without AR(1) disturbances         
  x <- cbind(ones(t,1), m )
  b1 <- lm(y ~ x - 1)$coef

  v <- y - x %*% b1
  z <- cbind(ones(t-1,1), trimr(m,1,0), trimr(v,0,1))
  v <- trimr(v,1,0)       
  b2 <- lm(v ~ z - 1)$coef 

  e <- v - z %*% b2            
  r2 <- 1 - t(e) %*% e/( t( (v-mean(v)) ) %*% (v-mean(v)))                      
  LM  <- (t-1)*r2
      
  cat('\nLM test applied to the CAPM without AR(1) disturbances\n')
  cat('\nLM test                   = ', LM)
  cat('\np-value                   = ', 1-pchisq(LM,1), '\n')

  # Estimate a CAPM model with AR(1) disturbances
  theta0 <- c(0.1 , 0.1 , 0.0)
  estResults <- optim(theta0, neglog, y=y, m=m, method = "BFGS", hessian=TRUE)
  theta <- estResults$par
  H <- estResults$hessian
    
  # Compute residuals
  u  <- y - theta[1] - theta[2]*m
  v  <- trimr(u,1,0) - theta[3]*trimr(u,0,1)
  z  <- cbind(ones(t-2,1), trimr(m,2,0), trimr(v,0,1))
  v  <- trimr(v,1,0)       
  b2 <- lm(v ~ z - 1)$coef                
  e  <- v - z %*% b2            
  r2 <- 1 - t(e) %*% e/(t( (v-mean(v))) %*% (v-mean(v)))                      
  LM <- (t-1)*r2
 
  cat('\nLM and Wald tests applied to the CAPM with AR(1) disturbances')
  cat('\nLM test                   = ', LM)
  cat('\np-value                   = ', 1-pchisq(LM,1), '\n')

  # Wald test (applied to the model with AR(1) disturbances)     
  r <- rbind(c(0, 0, 1))
  q <- 0
  W <- t* t( (r %*% theta - q) ) %*% inv(r %*% (H) %*% t(r)) %*% (r %*% theta - q)

  cat('\nWald test                 = ', W)
  cat('\np-value                   = ', 1-pchisq(W, 1))  
}



