#============================================================================
#
#      Simulation of data from an autoregressive model with heteroskedasticity
#      Estimate by OLS and generate various standard errors
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
# 

# load required functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Negative unconstrained log-likelihood
#----------------------------------------------------------------------------
neglog <- function( b,y,x ) {  
    logl <- -mean( lnlt( b,y,x ) )
    return(logl)
}

#----------------------------------------------------------------------------
#   Log-likelihood at every observation
#----------------------------------------------------------------------------
lnlt <- function( b,y,x ){
  u     <- y - b[1] - b[2]*x
  sig2  <- exp( b[3] ) 
  z     <- u/sqrt( sig2 )
  loglt <- - 0.5*log( 2*pi ) - 0.5*log( sig2 ) - 0.5*z^2
  return(loglt)  
}

#
#--------------------  QMLE Simulation Experiment ---------------------------
#

    
qmle_simul_ols <- function() {
  t     <- 500                
  beta1 <- 1.0   
  beta2 <- 0.5 
  rho   <- 0.5
  gam1  <- 1.0  
  gam2  <- 0.5

  # Construct exogenous variables 
  tmp <- seq(t)
  x <- 0.01*tmp + 10*runif(t)*cos( 2*pi*t/40 )

  # Simulate the AR(1) regression model with heteroskedasticity    
  y <- rnorm(t)
  u <- rnorm(t)
  v <- rnorm(t)

  for (i in 2:t) {
    sig2 <- exp( gam1 + gam2*x[i] )
    v[i] <- rnorm(1)*sqrt(sig2)
    u[i] <- rho*u[i-1] + v[i]
    y[i] <- beta1 + beta2*x[i] + u[i]
  }

  # Estimate the model by OLS
  X     <- cbind(rep(1, t), x)
  bols  <- lm(y ~ X - 1)$coef
  e     <- y - X %*% bols
  s2    <- mean(e^2)
  omols <- s2*inv(t(X) %*% X)
  seols <- sqrt(diag(omols)) 

  # Estimate the model by ML
  bstart <- c(bols,  log( s2 ))

  # Optimization settings
  estResults <- optim(bstart, neglog, y=y, x=x, method="BFGS", hessian=T)
  bhat <- estResults$par
  H <- estResults$hessian

  # Computing Standard Errors
  # Based on the Hessian
  iH  <- inv( H )
  seH <- sqrt( (1/t)*diag( iH ) )

  # Based on the OPG matrix
  g   <- numgrad(lnlt,bhat,y,x )
  j   <- t(g) %*% g/t
  seG <- sqrt( (1/t)*diag(inv( j )) )

  # Based on QMLE
  seQ <- sqrt( (1/t)*diag( iH %*% j %*% iH ) )

  cat('\nStandard Errors \n')
  cat('\n         Hessian     OPG         QMLE')
  cat('\n-------------------------------------------\n')
  print(unname(cbind(seH, seG, seQ)))

  #  Compute qmle standard errors with lag length based on p 
  p <- floor( 4*(t/100)^(2/9) )
 
  for (i in seq(p)) {
    gmat <- t( g[(i+1):t,] ) %*% g[1:(t-i),]/t
    j    <- j + (1.0 - i/(p+1))*(gmat + t(gmat))    
  }

  seQp <- sqrt( (1/t)*diag( iH %*% j %*% iH) )
  
  cat('\nStandard errors (QMLE - p)')
  cat('\n-------------------------------------------\n')
  cat(seQp, '\n')

  cat('\nStandard errors OLS')
  cat('\n-------------------------------------------\n')
  cat(seols, '\n')  
}
