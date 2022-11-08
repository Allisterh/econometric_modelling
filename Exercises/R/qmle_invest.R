#=========================================================================
#
#      Estimate investment model and compute QMLE standard errors
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions-------------------------------
#

#load required functions - trimr, inv
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Negative unconstrained log-likelihood function
#-------------------------------------------------------------------------
neglog <- function( b,ri,ry,rint ) {
  logl <- -mean( lnlt( b,ri,ry,rint ) )
  return(logl)
}

#-------------------------------------------------------------------------
#  Conditional log-likelihood function 
#-------------------------------------------------------------------------
lnlt <- function(b,ri,ry,rint) {
  beta0 <- b[1]                                          
  beta1 <- b[2]
  beta2 <- b[3]    
  sig2  <- b[4]
  u     <- ri - beta0 - beta1*ry - beta2*rint
     
  #  Log-likelihood for t=2,3,...   
  lf <- - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*u^2/sig2
  return(lf) 
}

#
#---------------------- U.S. Investment Model ---------------------------
#

qmle_invest <- function () {
  # Load data
  load('investment.Rdata')
    
  cpi    <- usinvest[,1]
  gdp    <- usinvest[,2]
  invest <- usinvest[,3]
  r10yr  <- usinvest[,4]
  r3yr   <- usinvest[,5]
  tbill  <- usinvest[,6]
     
  # Generate variables: data start in 1958Q1
  gfc  <- cbind(c( rep(0, nrow(usinvest)-13), rep(1, 13) ))
  dri  <- 100*(log(trimr(invest/cpi,1,0)) - log(trimr(invest/cpi,0,1)))   
  inf  <- 100*log(trimr(cpi,1,0)/trimr(cpi,0,1))                           
  rint <- trimr(r10yr/4,1,0) - inf                                          
  dry  <- 100*( log(trimr(gdp/cpi,1,0)) - log(trimr(gdp/cpi,0,1)) )       
  gfc  <- trimr(gfc,1,0)

  # OLS regression
  t <- length(dry)
  y <- dri
  x <- cbind(rep(1, t), dry, rint)
  reg <- lm(y ~ x - 1)
  
  b  <- reg$coef
  e  <- reg$residuals  
  s2 <- mean(e^2)
    
  # Conditional MLE 
  theta0  <- c(b, s2)
  estResults <- optim(theta0, neglog, ri=dri, ry=dry, rint=rint, method="BFGS", hessian=T)
  theta <- estResults$par
  fc <- estResults$value
  H <- estResults$hessian

  iH  <- inv(H)
  seH <- sqrt((1/t)*diag(iH))

  cat('\nResults for conditional MLE')
  cat('\nLog-likelihood function = ', -fc)
  cat('\nParameter estimates and std. errors\n')
  print( unname( cbind(theta, seH) ) )

  # Estimate the gradient and Hessian using the unconstrained parameters 
  g    <- numgrad(lnlt,theta,dri,dry,rint )

  # Compute standard errors based on the OPG matrix
  j   <- t(g) %*% g/t
  seG <- sqrt( (1/t)*diag( inv( j ) ) )
 
  # Compute QMLE standard errors  
  seQ0 <- sqrt( (1/t)*diag( iH %*% j %*% iH) )

  #  Compute qmle standard errors with lag length based on p 
  p <- floor( 4*(t/100)^(2/9) )
 
  tmax <- nrow(g) 
  for (i in seq(p)) {
    gmat <- t( g[(i+1):tmax,] ) %*% g[1:(tmax-i),]/t
    j    <- j + (1.0 - i/(p+1))*(gmat + t(gmat))    
  }

  seQp <- sqrt( (1/t)*diag( iH %*% j %*% iH) )

  cat('\n\nParameter Estimates')
  cat('\n-------------------')
  cat('\n', theta)

  cat('\n\nStandard errors (Hessian)')
  cat('\n-------------------------')
  cat('\n', seH )

  cat('\n\nStandard errors (OPG)')
  cat('\n-------------------------')
  cat('\n', seG )

  cat('\n\nStandard errors (QMLE - p=0)')
  cat('\n-------------------------')
  cat('\n', seQ0 )


  cat('\n\nStandard errors (QMLE - p=4)')
  cat('\n-------------------------')
  cat('\n', seQp )  
}
 
    

