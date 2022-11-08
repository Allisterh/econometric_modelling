#=================================================================================
#     Program to estimate a nonlinear consumption function
#     using the GAUSS-NEWTON algorithm
#
#     U.S. data on real consumption and real disposable income(2005 $)
#     1960:Q1 to 2009:Q4 200 observations
#================================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
#

# load required functions - inv
source("EMTSUtil.R")

#
#--------------------------- Nonlinear Consumption ------------------------------
#

nls_conest <- function() {  
  # Load data
  # [cons, inc]
  load('USdata.Rdata')
  
  yt <- USdata[,2]
  ct <- USdata[,1]
  t <- length(yt)
  
  # Estimate the linear model to use as initial starting values  
  b <- lm(ct ~ cbind(rep(1, t), yt) - 1)$coef
  
  # GAUSS-NEWTON algorithm   
  alpha <- b[1]
  beta  <- b[2]
  gam   <- 1.00
  
  crit  <- 0.00001                        # Convergence criterion         
  maxit <- 20                             # Maximum number of iterations  
  for (i in seq(maxit)) {
    et  <- ct - alpha - beta*yt^gam          
    z1t <- rep(1,t)                         
    z2t <- yt^gam
    z3t <- beta*log(yt)*yt^gam

    cat('\n' )
    cat('\nIteration        = ', i)
    cat('\nParameters       = ', alpha, beta, gam)
    cat('\n---------------------------------------------------------------\n')
    y       <- et
    x       <- cbind(z1t, z2t, z3t)
    bchange <- lm(y ~ x - 1)$coef

    if ( (t(bchange) %*% bchange) < crit) {
      break
    } else {
      alpha <- alpha + bchange[1]     # Update new parameters   
      beta  <- beta  + bchange[2]
      gam   <- gam   + bchange[3]

      if (i == maxit)
        cat('\nFailed to converge after iteration: ', maxit, '\n' )        
    }    
  }

  ssr  <- t(y) %*% y
  sig2 <- ssr/t
  cat('\nSum squared residuals = ', ssr)
  cat('\nResidual variance     = ', sig2)          
  cat( '\n\nInformation matrix\n' )
  z <- cbind(z1t, z2t, z3t)
  im <- t(z) %*% z/c(sig2)
  print(im)
  
  cat('\n\nEstimated asymptotic covariance matrix\n')
  print( inv(im) )
  
  bhat <- c(alpha, beta, gam)
  se   <- sqrt(diag(inv(im)))
  cat('\nParameter estimates  = ', bhat )
  cat('\nStandard errors      = ', se )
  cat('\nt-statistics         = ', bhat/se )
}
