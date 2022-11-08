#============================================================================
#
#   Program to estimate a probit model of US monetary policy
#   The data file is based on the Hamilton and Jorda (2002) data set.
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  inv, figure, seqa
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Unrestricted Probit negative log-likelihood function 
#----------------------------------------------------------------------------
lprobit <- function (b,y,x) {
  f  <- pnorm(x %*% b)
  lf <- -mean( y*log(f) + (1 - y)*log(1 - f) )
  return (lf)
}

#----------------------------------------------------------------------------
#  Restricted Probit negative log-likelihood function 
#----------------------------------------------------------------------------
l0probit <- function(b,y) {  
  f  <- pnorm(b) 
  lf <- -mean( y*log(f) + (1 - y)*log(1 - f) )
  return(lf)
}


#
#---------------- Probit Model for Monetory Policy  -------------------------
#

discrete_probit <- function() {
  # Read the data: US weekly yields   
  # from 1st week of February 1984 to first week of June 1997
  
  usmoney <- as.matrix(read.table("usmoney.dat"))
  #event   <- usmoney[,1]
  target  <- usmoney[,2]
  #change  <- usmoney[,3]
  bin     <- usmoney[,4]
  fomc    <- usmoney[,8]
  spread6 <- usmoney[,20]
  inf     <- usmoney[,23]
  #unemp   <- usmoney[,25]
  gdp     <- usmoney[,28]
  
  # Reverse the spread so it is the Federal funds rate 
  # less 6-month Treasury bill rate       
  spread <- -spread6                                                  
  
  #  Redefine the target rate based on the consolidated series 
  # constructed in Hamilton and Jorda (2002)       
  target_adj <- cumsum( c(target[1], bin[2:length(bin)] ))
  
  figure()
  plot(seqa(1984+5/52,1/52,length(target_adj)),target_adj, type="l",
       main="Federal Funds target rate(%) from 1984 to 1997", 
       ylab = "Federal Funds Rate (%)",
       xlab = "Time")
  
  # Choose data based on fomc days
  ind       <- fomc == 1
  data      <- cbind(bin, spread, inf, gdp)
  data_fomc <- data[ind,]
  
  # Dependent and independent variables 
  y <- as.numeric(data_fomc[,1] > 0.0)
  t <- length(y)
  x <- cbind(rep(1, t), data_fomc[,2:4])
  
  # Estimate model by OLS (ie ignore that the data are binary
  reg <- lm(y ~ x - 1)
  b <- reg$coef
  u <- reg$residuals
  s <- sqrt( mean(u^2) )
  
  theta0 <- b/s
  
  estResults <- optim(theta0, lprobit, x=x, y=y, method="BFGS", hessian=T)
  
  theta1 <- estResults$par
  l1 <- estResults$val
  h <- estResults$hessian
  
  cov <- (1/t)*inv(h)
  l1 <- -l1
  
  cat('\nUnrestricted log-likelihood function =     ',l1)
  cat('\nT x unrestricted log-likelihood function = ',t*l1)
  
  cat('\n\nUnrestricted parameter estimates')
  cat('\n',theta1)
  
  # Estimate the restricted probit regression model by MLE     
  theta0 <- qnorm(mean(y),0,1)
  estResults <- optim(theta0, l0probit, y=y, method="BFGS", hessian=T)
  
  theta <- estResults$par
  l0 = estResults$val
  h <- estResults$hessian
  
  l0 <- -l0
  cat('\n\nRestricted log-likelihood function =     ',-l0)
  cat('\nT x restricted log-likelihood function = ',-t*l0)
  
  cat('\n\nRestricted parameter estimates')
  cat('\n', theta)
  #  Likelihood ratio test  
  lr <- -2*t*(l0 - l1)
  
  cat('\n\nLR Statistic         = ',lr)
  cat('\np-value              = ',1-pchisq(lr,ncol(x)-1))
  
  cat('\n')
  
  # Wald test  
  r <- matrix(c(0,   1,   0,   0,
                0,   0,   1,   0,
                0,   0,   0,   1), byrow=T, nrow=3)
  
  q <- rbind(0,0,0)
  
  wd <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov %*% t(r)) %*% (r %*% theta1 - q)
  
  cat('\nWald Statistic       = ',wd)
  cat('\np-value              = ',(1-pchisq(wd,ncol(x)-1)))
  
  cat('\n')
  
  # LM test of the joint restrictions  
  u  <- y - mean(y)
  b  <- lm(u ~ x - 1)$coef
  e  <- u - x %*% b
  r2 <- 1 - (t(e) %*% e)/(t(u) %*% u)
  lm <- t*r2
  
  cat('\nRegression estimates = ', b)
  cat('\nSample size (t)      = ', t )
  cat('\nR2                   = ', r2 )
  cat('\nLM Statistic         = ',lm)
  cat('\np-value              = ',1-pchisq(lm,ncol(x)-1))  
}
