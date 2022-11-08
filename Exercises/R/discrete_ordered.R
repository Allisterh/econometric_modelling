#=========================================================================
#
#   Program to estimate an ordered probit model of US monetary policy
#   The data file is based on the Hamilton and Jorda (2002) data set.
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Unrestricted Probit negative log-likelihood function 
#----------------------------------------------------------------------------
lprobit <- function(b,x,d) {
  # Cut off points
  c <- b[1:4]

  # Regression part excluding the intercepts
  xb <- x[,2:4] %*% b[5:7]
  
  # Cut off points
  f1 <- pnorm( c[1] - xb,0,1)
  f2 <- pnorm( c[2] - xb,0,1) - f1
  f3 <- pnorm( c[3] - xb,0,1) - f1 - f2
  f4 <- pnorm( c[4] - xb,0,1) - f1 - f2 - f3
  f5 <- 1 - f1 - f2 - f3 - f4
  f  <- cbind(f1,  f2,  f3,  f4,  f5)   
  
  # Log-likelihood function  
  tp <- d * log(f)  
  lf <- -mean( rowSums(tp) )  
  return(lf)
}


#----------------------------------------------------------------------------
#  Restricted Probit negative log-likelihood function 
#----------------------------------------------------------------------------
l0probit <- function(b,d) {
  # Cut off points
  f1 <- pnorm( b[1],0,1)
  f2 <- pnorm( b[2],0,1) - f1
  f3 <- pnorm( b[3],0,1) - f1 - f2
  f4 <- pnorm( b[4],0,1) - f1 - f2 - f3
  f5 <- 1 - f1 - f2 - f3 - f4
  f  <- cbind(f1, f2,  f3,  f4,  f5)
    
  # Log-likelihood function 
  tp <- t(apply(d, 1, '*', log(f)))
  lf <- -mean( rowSums(tp) )
  return(lf)
}


#
#---------------- Ordered Probit Model for Monetory Policy  -----------------
#

discrete_ordered <- function() {
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
  
  # Choose data based on fomc days
  ind       <- fomc == 1
  data      <- cbind(bin, spread, inf, gdp)
  data_fomc <- data[ind,]
  
  # Dependent and independent variables 
  y <- data_fomc[,1]
  t <- length(y)
  x <- cbind(rep(1, t), data_fomc[,2:4])   
  
  # Create dummy variables for each interest rate change
  d1 <- as.numeric(y == -0.50)
  d2 <- as.numeric(y == -0.25)
  d3 <- as.numeric(y ==  0.00)
  d4 <- as.numeric(y ==  0.25)
  d5 <- as.numeric(y ==  0.50)
  
  d  <- cbind(d1, d2, d3, d4, d5)
  
  # Estimate model by OLS (ie ignore that the data are ordered
  reg <- lm(y ~ x - 1)
  b <- reg$coef
  u <- reg$residuals
  s <- sqrt( mean(u^2) )
  
  # Compute the unconditional probabilities of each regime    
  p <-  cumsum( colMeans(d) ) 
  
  # Estimate the unrestricted ordered probit regression model by MLE   
  theta0 <- c(qnorm(p[1:4],0,1), b[2:4]/s )
  estResults <- optim(theta0, lprobit, x=x, d=d, method="BFGS", hessian=T)
  theta1 <- estResults$par
  l1 <- estResults$val
  h <- estResults$hessian
  
  cat('\nUnrestricted parameter estimates')
  cat('\n', theta1 )
  
  cov <- (1/t)*inv(h)
  l1 <- -l1
  
  cat('\n\nUnrestricted log-likelihood function =     ',l1)
  cat('\nT x unrestricted log-likelihood function = ',t*l1)
  
  # Estimate the restricted probit regression model by MLE     
  theta0 <- qnorm(p[1:4],0,1)
  estResults <- optim(theta0, l0probit, d=d, method="BFGS", hessian=T)
  theta <- estResults$par
  l0 <- estResults$val
  h <- estResults$hessian
  
  cat('\n\nRestricted parameter estimates')
  cat('\n', theta )
  
  l0 <- -l0
  cat('\n')
  cat('\nRestricted log-likelihood function =     ',-l0)
  cat('\nT x restricted log-likelihood function = ',-t*l0)
  
  #  Likelihood ratio test  
  lr <- -2*t*(l0 - l1)
  
  cat('\n\nLR Statistic         = ',lr)
  cat('\np-value              = ',1-pchisq(lr,ncol(x)-1))
  
  cat('\n')
  
  # Wald test  
  r <- matrix(c(0,   0,   0,   0,   1,   0,   0,
                0,   0,   0,   0,   0,   1,   0,
                0,   0,   0,   0,   0,   0,   1), byrow=T, nrow=3)
  
  q <- rbind(0, 0, 0)
  
  wd <- t( r %*% theta1 - q) %*% inv(r %*% cov %*% t(r)) %*% (r %*% theta1 - q)
  
  cat('\nWald Statistic       = ',wd)
  cat('\np-value              = ',1-pchisq(wd,ncol(x)-1))
  
}

