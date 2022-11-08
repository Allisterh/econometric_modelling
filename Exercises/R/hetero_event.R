#============================================================================
#
#   A heteroskedastic model of money shocks in US asset markets
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv, trimr
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Unrestricted log-likelihood function
#-------------------------------------------------------------------------
neglog1 <- function(b,y,d) {  
  m  <- b[1] + b[2]*d         # Mean
  s2 <- exp(b[3] + b[4]*d)    # Variance  
  lt <- -0.5*log(2*pi) -0.5*log(s2) - 0.5*((y - m)^2)/s2        
  lf <- -mean(lt)
  return(lf)
}

#-------------------------------------------------------------------------
# Restricted log-likelihood function
#-------------------------------------------------------------------------
neglog0 <- function(b,y,d) {
  m  <- b[1] + b[2]*d         # Mean
  s2 <- exp(b[3])             # Variance
  
  lt <- -0.5*log(2*pi) -0.5*log(s2) - 0.5*((y - m)^2)/s2        
  lf <- -mean(lt)
  return(lf)
}


hetero_event <- function( ) {
  # Load data
  ydata <- read.table("usasset.dat")
  
  #   1.  tcm3m           (from FRB webpage)
  #  2.	tcm6m  
  #	3.	tcm1y  
  #	4.	tcm2y 
  #	5.	tcm3y
  #	6.	tcm5y
  #	7.	tcm7y
  #   8.  tcm10y
  #   9.  edm1m
  #   10. xrate ($US/$AU)
  #	11.	frb target in US
  #   12. change in target variable US
  #   13. dummy variable corresponding to event dates US (taken from Poole and Kuttner)
  #   14. value-weighted returns in US    (from CRSP)
  #   15. value-weighted index US         (from CRSP)
  
  # Different to Poole dating
  ydata[1274,13] <- 1.0
  ydata[1275,13] <- 0.0
  
  # Define variables
  yield    <- ydata[,1:8]  # 3m, 6m 1y, 2y, 3y, 5y, 7y, 10y as percentages      
  euro     <- ydata[,9]                 # 1m  euro rate as  percentage                        
  event    <- ydata[,13]                 # dummy variable on meeting dates                                 
  nyse     <- ydata[,15]                 # NYSE equity price index                                         
  
  # Dependent variables
  dyield <- (trimr(yield,1,0) - trimr(yield,0,1))            
  deuro  <- trimr(euro,1,0) - trimr(euro,0,1)                  
  rnyse  <- 100*(trimr(log(nyse),1,0) - trimr(log(nyse),0,1))   
  
  y <- dyield[,1]
  d <- trimr(event,1,0)   # Event day
  x <- cbind(rep(1, length(y)),  d)
  t <- length(y)
  
  # Compute descriptive statistics on non-event and event days    
  y_nonevent <- y[d == 0]
  y_event    <- y[d == 1]
  
  cat('\nSample mean (event)         = ',mean(y_event))
  cat('\nSample mean (non-event)     = ',mean(y_nonevent))
  cat('\nSample variance (event)     = ',mean((y_event-mean(y_event))^2)) 
  cat('\nSample variance (non-event) = ',mean((y_nonevent-mean(y_nonevent))^2))
  
  #  Estimate the unrestricted model by MLE with starting values based on OLS
  bols <- lm(y ~ x - 1)$coef
  u    <- y - x %*% bols
  s2   <- mean(u^2) 
  
  start <- c(bols,  s2,  0 )
  estResults <- optim(start, neglog1, y=y, d=d, method="BFGS", hessian=T)
  bhat1 <- estResults$par
  lf1 <- estResults$val
  hess <- estResults$hess
  
  lf1 <- -lf1
  vc1 <- (1/t)*inv(hess)
  cat('\n')
  cat('\nUnrestricted model')
  cat('\nMLE estimate of the mean       (event days)     = ',(bhat1[1] + bhat1[2]))
  cat('\nMLE estimate of the mean       (non-event days) = ',bhat1[1])
  cat('\nMLE estimate of the volatility (event days)     = ',exp(bhat1[3] + bhat1[4])) 
  cat('\nMLE estimate of the volatility (non-event days) = ',exp(bhat1[3]))
  cat('\n')
  
  #  Estimate the restricted model by MLE with starting values based on OLS
  start <- c(bols,  s2)
  estResults <- optim(start, neglog0, y=y, d=d, method="BFGS", hessian=T)
  bhat0 <- estResults$par
  lf0 <- estResults$val
  hess <- estResults$hess
  
  lf0 <- -lf0
  vc0 <- (1/t)*inv(hess)
  
  cat('\n')
  cat('\nRestricted model')
  cat('\nMLE estimate of the mean       (event days)     = ',(bhat0[1] + bhat0[2])) 
  cat('\nMLE estimate of the mean       (non-event days) = ',bhat0[1])
  cat('\nMLE estimate of the volatility (event days)     = ',exp(bhat0[3]))
  cat('\nMLE estimate of the volatility (non-event days) = ',exp(bhat0[3]))
  cat('\n')
  
  # LR test   
  lr <- -2*t*(lf0 - lf1)
  cat('\nLR statistic            = ',lr)
  cat('\np-value                 = ',1-pchisq(lr,1))
  
  # Wald test 
  r <- cbind(0 , 0 , 0 , 1)
  q <- 0
  wd <- t(r %*% bhat1 - q) %*% inv(r %*% vc1 %*% t(r)) %*% (r %*% bhat1 - q)
  cat('\nWald statistic          = ',wd)
  cat('\np-value                 = ',1-pchisq(wd,1))
  
  # LM test (regression form)                                   
  # Stage 1 regression
  b <- lm(y ~ x - 1)$coef
  u <- y - x %*% b    
  w <- x                               
  v <- u^2
  # Stage 2 regression
  b  <- lm(v ~ w - 1)$coef
  e  <- v - w %*% b
  r2 <- 1 - sum(e^2)/sum( (v-mean(v))^2 )
  lm <- t*r2
  cat('\nLM statistic (regression) = ',lm)
  cat('\np-value                   = ',1-pchisq(lm,1))
}



