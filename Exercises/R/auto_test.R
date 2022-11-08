#=========================================================================
#
#     Program to estimate and test a dynamic investment model
#     U.S. quarterly data for period 1957 to 2007
#     
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions --------------------------
#
#load required functions - inv, trimr, numgrad
source("EMTSUtil.R")
library("matlab")

#-----------------------------------------------------------------------
# Negative unconstrained log-likelihood  
#-----------------------------------------------------------------------
neglog1 <- function(b,ri,ry,rint) {
  lf <- -mean( lnlt1(b,ri,ry,rint) )
  return(lf)
}

#-----------------------------------------------------------------------
# Unconstrained log-likelihood function at each observation
#-----------------------------------------------------------------------
lnlt1 <- function(b,ri,ry,rint) {
  beta0 <- b[1]                                          
  beta1 <- b[2]
  beta2 <- b[3]    
  rho1  <- b[4]
  sig2  <- b[5]

  u     <- ri - beta0 - beta1*ry - beta2*rint
  v     <- trimr(u,1,0) - rho1*trimr(u,0,1)
     
  #  Log-likelihood for t=2,3,...   
  lnl <- - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2
  return(lnl)  
}

#-----------------------------------------------------------------------
# Negative constrained log-likelihood function
#-----------------------------------------------------------------------
neglog0 <- function(b,ri,ry,rint) {
  beta0 <- b[1]                                                    
  beta1 <- b[2]
  beta2 <- b[3]
  rho1  <- 0.0
  sig2  <- b[4]
  u     <- ri - beta0 - beta1*ry - beta2*rint
  v     <- trimr(u,1,0) - rho1*trimr(u,0,1)
     
  # Log-likelihood for t=2,3,...
  lnl <- - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2     
  lf  <- - mean(lnl)  
  return(lf)
}

#
#--------------------------- Model Tests -----------------------------
#

auto_test <- function() {
    # Load data
  data <- read.table("usinvest.txt")
    
  cpi    <- data[,1] 
  gdp    <- data[,2]
  invest <- data[,3]
  r10yr  <- data[,4]
  r3yr   <- data[,5]
  tbill  <- data[,6]
     
  # Generate variables: data start in 1958Q1    
  gfc  <- cbind( c(zeros(nrow(data)-13,1),  ones(13,1)) )
  dri  <- 100*(log( trimr(invest/cpi,1,0) ) - log( trimr(invest/cpi,0,1) ))   
  inf  <- 100*log(trimr(cpi,1,0)/trimr(cpi,0,1))                           
  rint <- trimr(r10yr/4,1,0) - inf                                         
  dry  <- 100*( log(trimr(gdp/cpi,1,0)) - log(trimr(gdp/cpi,0,1)) )       
  gfc  <- trimr(gfc,1,0)
    
  # OLS regression
  t <- length(dry)
  y <- dri
  x <- cbind(ones(t,1),dry,rint)
  
  b  <- lm(y ~ x - 1)$coef
  e  <- y - x %*% b
  s2 <- mean(e^2)

  # Unconstrained model
  theta  <- c(b, 0.02, sqrt(s2))
  estResults <- optim(theta, neglog1, ri=dri, ry=dry, rint=rint, method="BFGS", hessian=T)
  theta1 <- estResults$par
  f1 <- estResults$value
  H1 <- estResults$hessian

  lnl1  <- -f1
  invH1 <- inv(H1)
    
  cat('\nResults for conditional MLE')
  cat('\nLog-likelihood function = ', lnl1)
  cat('\nParameter estimates and std. errors\n')
  sterr <- diag(invH1)/t
  print( unname(cbind(theta, sterr)) )

  # Constrained model 
  theta <- c(b, sqrt(s2))
  estResults <- optim(theta, neglog0, ri=dri, ry=dry, rint=rint, method="BFGS", hessian=T)
  theta0 <- estResults$par
  f0 <- estResults$value
  H0 <- estResults$hessian
  
  invH0 <- inv(H0)
  lnl0  <- -f0
    
  cat('\nResults for conditional MLE')
  cat('\nLog-likelihood function = ', lnl0)
  cat('\nParameter estimates and std. errors\n')
  sterr <- diag(invH0)/t
  print( unname(cbind(theta0, sterr)) )

  # LR test      
  lr <- -2*t*(lnl0 - lnl1)
  cat('\nLR statistic            = ',lr) 
  cat('\np-value                 = ',1-pchisq(lr, 1))

  # Wald test   
  r <- rbind(c( 0 , 0 , 0 , 1 , 0))
  q <- 0
  wd <- t* t( (r %*% theta1 - q) ) %*% inv(r %*% inv(H1) %*% t(r)) %*% (r %*% theta1 - q)
  cat('\nWald statistic          = ', wd) 
  cat('\np-value                 = ', 1 - pchisq(wd, 1))

  # LM test      
  theta <- c(theta0[1:3],0.0, theta0[4])
  gmat  <- numgrad(lnlt1,theta,dri,dry,rint)  
  g      <- cbind(colMeans(gmat))
  j      <- t(gmat) %*% gmat/t
  lm     <- t* t( g ) %*% inv(j) %*% g

  cat('\nGradient evaluated at contrained estimates\n')
  print(g)
  cat('\nOuter product of gradients matrix\n')

  print(j)     
  cat('\nLM statistic               = ', lm) 
  cat('\np-value                    = ', 1-pchisq(lm, 1))

  # LM test (regression form)                                
  # Stage 1 regression
  #b <- lm(y ~ x - 1)$coef
  #u <- y - x %*% b  
  u <- lm(y ~ x - 1)$res
  
  # Stage 2 regression
  y <- trimr(u,1,0)
  z <- cbind(trimr(x,1,0), trimr(e,0,1))

  #b <- lm(y ~ z - 1)$coef
  #e  <- y - z %*% b
  e <- lm(y ~ z - 1)$res
  
  ete <- t(e) %*% e
  ymy <- (t(y-mean(y)) %*% (y-mean(y)))
  r2  <- 1 - lm(ete ~ ymy - 1)$coef
  lm <- (t-1)*r2
  cat('\nLM statistic (regression)  = ', lm) 
  cat('\np-value                    = ', 1-pchisq(lm, 1))
  
  #  LM test (first-order autocorrelation of residuals)           
  y <- dri
  u <- lm(y ~ x - 1)$res
  r1 <- acf(u,1, plot=FALSE)$acf     
  lm <- (t-1)*r1[2]^2

  cat('\nLM statistic (alternative) = ', lm) 
  cat('\np-value                    = ', 1-pchisq(lm, 1))  
}
