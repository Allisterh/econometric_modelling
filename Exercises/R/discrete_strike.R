#============================================================================
  #
#   Program to estimate a duration model of strikes
#   using the Kennan (1985) data set.
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
# 
#----------------------------------------------------------------------------
#  Unrestricted log-likelihood function of duration model
#----------------------------------------------------------------------------
neglog <- function(b,y,x) { 
  m  <- x %*% b[1:2]   
  s2 <- abs(b[3])   
  z  <- (log(y) - m)/sqrt(s2)
  lf <- -mean( 0.5*log(s2) + z - exp(z) )    
  return(lf)
}


#----------------------------------------------------------------------------
#  Restricted log-likelihood function of duration model
#----------------------------------------------------------------------------
neglog0 <- function(b,y,x) {
  m  <- x %*% b[1:2]   
  s2 <- 1   
  z  <- (log(y) - m)/sqrt(s2)
  lf <- -mean( 0.5*log(s2) + z - exp(z) )
  return(lf)
}

#
#------------------- A Duration Model of Strikes  ---------------------------
#
#discrete_strike <- function() {
  # Load data: US strike data from 1968 to 1976
  data <- as.matrix(read.table("strike.dat"))
  
  y   <- data[,1]    # Duration of strike in days    
  out <- data[,2]    # Unanticipated output     
  t   <- length(y)
  x   <- cbind(rep(1,t), out)
  
  # Estimate the Weibull model by MLE 
  reg <- lm(log(y) ~ x - 1)$coef
  theta_0       <- c(reg,  1)
  
  estResults <- optim(theta_0, neglog, y=y, x=x)
  theta1 <- estResults$par
  l1 <- estResults$val
  
  l1 <- -l1   
  cat('\nLog-likelihood (unrestricted) = ',l1)
  cat('\n(T-1)xLog-likelihood        = ',(t-1)*l1)
  cat('\n')
  
  # Estimate the Exponential model by MLE
  theta_0       <- c(reg)
  
  estResults <- optim(theta_0, neglog0, y=y, x=x)
  theta0 <- estResults$par
  l0 <- estResults$val
  
  l0 <- -l0   
  cat('\nLog-likelihood (restricted) = ',l0)
  cat('\n(T-1)xLog-likelihood        = ',(t-1)*l0)
  cat('\n')
  
  # Likelihood ratio test 
  lr <- -2*t*(l0 - l1)
  cat('\nLR Statistic               = ',lr)
  cat('\np-value                    = ',1-pchisq(lr,ncol(x)-1))
#}


