#========================================================================================
#
#     Program to estimate a dynamic investment model
#     U.S. quarterly data for period 1957 to 2007
#
#=========================================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions ---------------------------------------------
#
#load required functions - inv, trimr
source("EMTSUtil.R")

#------------------------------------------------------------------------------------------
#   Exact log-likelihood function 
#------------------------------------------------------------------------------------------
neglog <- function(b,ri,ry,rint) {
  beta0 <- b[1]                                          
  beta1 <- b[2]
  beta2 <- b[3]
  rho1  <- b[4]
  sig2  <- b[5]

  u     <- ri - beta0 - beta1*ry - beta2*rint
  v     <- trimr(u,1,0) - rho1*trimr(u,0,1)
     
  #  Log-likelihood for t=1              
  lnl_0 <- -0.5*log(2*pi) - 0.5*log(sig2) + 0.5*log(1 - rho1^2) - 0.5*(u[1] - 0)^2/(sig2/(1 - rho1^2))   
  #  Log-likelihood for t=2,3,...   
  lnl_1 <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2                                                  

  lf <- -mean( c(lnl_0, lnl_1))
  return(lf)
}

#------------------------------------------------------------------------------------------
#  Conditional log-likelihood function 
#------------------------------------------------------------------------------------------
neglogc <- function(b,ri,ry,rint) {
  beta0 <- b[1]                                          
  beta1 <- b[2]
  beta2 <- b[3] 
  rho1  <- b[4]     
  sig2  <- b[5]

  u     <- ri - beta0 - beta1*ry - beta2*rint
  v     <- trimr(u,1,0) - rho1*trimr(u,0,1) 
     
  #  Log-likelihood for t=2,3,...   
  lnl_1 <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2                                                  

  lf <- -mean( lnl_1 )
  return (lf)
}

  
     
#
#--------------------- Dynamic Model of U.S. Investment -------------------------------
#

auto_invest <- function() {
  # Load data
  data <- read.table("usinvest.txt")
    
  cpi    <- data[,1] 
  gdp    <- data[,2]
  invest <- data[,3]
  r10yr  <- data[,4]
  r3yr   <- data[,5]
  tbill  <- data[,6]
     
  # Generate variables: data start in 1958Q1    
  gfc  <- cbind( c(rep(0, nrow(data)-13),  rep(1, 13)) )
  dri  <- 100*(log( trimr(invest/cpi,1,0) ) - log( trimr(invest/cpi,0,1) ))   
  inf  <- 100*log(trimr(cpi,1,0)/trimr(cpi,0,1))                           
  rint <- trimr(r10yr/4,1,0) - inf                                         
  dry  <- 100*( log(trimr(gdp/cpi,1,0)) - log(trimr(gdp/cpi,0,1)) )       
  gfc  <- trimr(gfc,1,0)
    
  # OLS regression
  t <- length(dry)
  y <- dri
  x <- cbind(rep(1, t),dry,rint)

  b  <- lm(y ~ x - 1)$coef
  e  <- y - x*b
  s2 <- mean(e^2)
     
  # Exact MLE
  theta0  <- c(b, 0.02, sqrt(s2))   
  estResults <- optim(theta0, neglog, ri=dri, ry=dry, rint=rint, method="BFGS", hessian=TRUE)
  theta <- estResults$par
  fe    <- estResults$value
  H     <- estResults$hessian

  invH <- inv(H)
     
  cat('\nResults for exact MLE\n')
  cat('\nLog-likelihood function = ', -fe)
  cat('\nParameter estimates and std. errors\n')    
  print( unname(cbind(theta, diag(invH)/t) ))

  estResults <- optim(theta, neglogc, ri=dri, ry=dry, rint=rint, method="BFGS", hessian=TRUE)
  theta <- estResults$par
  fc    <- estResults$value
  H     <- estResults$hessian

  invH <- inv(H)

  cat('\nResults for conditional MLE\n')
  cat('\nLog-likelihood function = ', -fc)
  cat('\nParameter estimates and std. errors\n')
  print( unname(cbind(theta, diag(invH)/t) ))  
}

