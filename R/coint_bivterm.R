#============================================================================
#
#   Program to estimate bivariate term structure model
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr, rref
source("EMTSUtil.R")
# Load required library - flipud
library(matlab)

#-------------------------------------------------------------------------
#  Johansen procedure 
#----------------------------------------------------------------------------
johansen <- function(y,p,r,model) {
  ty <- nrow(y)
  dy <- trimr(y,1,0)-trimr(y,0,1) 
  z0 <- trimr(dy,p-1,0)
  z1 <- trimr(y,p-1,1)  
  
  z2 <- c()
  
#   for (j in 1:p-1) {
#     z2 <- cbind(z2, trimr(dy,p-1-j,j))
#   }
  if (model == 1) {
    z1 <- z1
    z2 <- z2
  } else if (model == 2) {
    z1 <- cbind(trimr(y,p-1,1),  rep(1, ty-p))
    z2 <- z2
  } else if (model == 3 ) {
    z1 <- z1
    z2 <- cbind(rep(1, ty-p), z2)    
  } else if (model == 4) {
    z1 <- cbind(z1, seqa(1,1,ty-p))
    z2 <- cbind(rep(1, ty-p),  z2) 
  }else if (model == 5) {
    z1 <- z1
    z2 <- cbind(rep(1, ty-p), seqa(1,1,ty-p),  z2)
  }  
  if (p == 1) {
    if (model <= 2) {
      r0 <- z0 
      r1 <- z1    
    }    
  }else {
    r0 <- lm(z0 ~ z2 - 1)$residuals 
    r1 <- lm(z1 ~ z2 - 1)$residuals
  }
  tr <- nrow(r0)
  tc <- ncol(r0)
  
  # Construct sums of squares matrices
  s00 <- t(r0) %*% r0/tr                                        
  s11 <- t(r1) %*% r1/tr
  s01 <- t(r0) %*% r1/tr
  s10 <- t(s01)
    
  # Sort eigenvectors on eigenvalue index and normalize
  l       <- t(chol(s11))
  eig <- eigen( inv(l) %*% s10 %*% inv(s00) %*% s01 %*% inv(t(l)) )
  tmp <- eig$values
  e <- eig$vectors  

  gam <- (inv(t(l)) %*% e)
  
  # Estimates
  beta  <- t(rref(t(gam[,1:r])))
  alpha <- t(lm(r0 ~ (r1 %*% beta) - 1)$coef)
  
    
  # Statistics
  logl   <- -0.5*tc*(log(2*pi)+1) - 0.5*log(det(s00)) - 0.5*sum( log(1-tmp[1:r]))
  tracet <- -tr*flipud(cumsum(log(1-flipud(tmp))))                     
  maxt   <- -tr*log(1-tmp) 
  return(list(alpha=alpha,beta=beta,logl=logl,maxt=maxt,tracet=tracet))
}


#-------------------------------------------------------------------------
#  Reduced rank log-likelihood function: Model 2 and p=1
#-------------------------------------------------------------------------
neglog <- function(b,z0,z1) {
  nobs <- nrow(z0)
  f   <- rep(0, nobs)
  
  v1 <- z0[,1] - b[3]*(z1[,1] - b[1] - b[2]*z1[,2])
  v2 <- z0[,2] - b[4]*(z1[,1] - b[1] - b[2]*z1[,2])
  v  <- cbind(v1,  v2)
  
  k <- nrow(v)
  n <- ncol(v)
  
  omegav <- t(v) %*% v/k

  for (t in seq(nobs)) {
    f[t] <- - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[t,] %*% inv(omegav) %*% cbind(v[t,])
  }  
  lf <- -mean( f )  
  return(lf)
}

# 
#-------------------- Term Structure of Interest Rates ----------------------
#
coint_bivterm <- function( ) {
  # Load data and choose 10-year and 1-year yields
  # data = r1yr, r5yr, r10yr
  data <- as.matrix(read.table("usmacro.dat"))
  yield <- cbind(data[,3], data[,1])
  t     <- nrow(yield)
  
  z0   <- trimr(yield,1,0) - trimr(yield,0,1)
  z1   <- trimr(yield,0,1)
  
  # Full rank case: Estimate the var in levels 
  x <- cbind(rep(1, t-1), trimr(yield,0,1))
  y <- trimr(yield,1,0)
  k <- nrow(y)
  n <- ncol(y)
  
  b      <- lm(y ~ x - 1)$coef
  v      <- y - x %*% b
  omegav <- t(v) %*% v/t
  
  lnl <- -0.5*k*n*(log(2*pi) + 1) - 0.5*k*log(det(omegav))
  cat('\nEstimated coefficients\n')
  print( b )
  cat('\n ' )
  cat('\nEstimated covariance matrix\n')
  print( omegav )
  
  cat('\nLog-likelihood of full rank model    = ', lnl/t )
  cat('\nDeterminant of omegav                = ', det(omegav))
  
  # Reduced rank case parameters
  p     <- 1      
  r     <- 1    
  model <- 2
  
  # Reduced rank case: iterative procedure (Model 2)
  # Starting values
  b <- lm(yield[,1] ~ cbind(rep(1, t), yield[,2]) - 1)$coef
  u <- yield[,1] - cbind(rep(1, t), yield[,2]) %*% b
  a <- lm(z0 ~ trimr(u,0,1) - 1)$coef
  
  pstart <- c(b, a)
  estResults <- optim(pstart, neglog, z0=z0, z1=z1, method="BFGS", hessian=T)   
  theta <- estResults$par
  logl <- estResults$val
  hess <- estResults$hess
  
  
  cat('\nIterative procedure estimates')        
  cat('\n-----------------------------')
  cat('\nLog likelihood function = ',-logl)
  cat('\nbeta_c         = ',theta[1])
  cat('\nbeta_r         = ',theta[2])
  cat('\nalpha_1        = ',theta[3]) 
  cat('\nalpha_2        = ',theta[4])
  vc <- (1/t)*inv(hess)
  cat('\nCovariance matrix (iterative estimator)\n')
  print(unname(vc))
  cat('\n ')
  
  # Perform Wald test of (1,-1) cointegrating vector)        
  wd <- (theta[2]-1)^2/vc[2,2]      
  
  cat('\nWald test of (1,-1) cointegrating vector')
  cat('\nWald statistic           = ',wd)
  cat('\np-value                  = ',1-pchisq(wd,1))
  cat('\n ')
  
  # Perform Wald test of y1 weakly exogenous    
  wd <- (theta[3]-0)^2/vc[3,3]                         
  
  cat('\nWald test of y1 weakly exogenous')
  cat('\nWald statistic           = ',wd)
  cat('\np-value                  = ',1-pchisq(wd,1))
  cat('\n ')
  
  # Perform Wald test of y2 weakly exogenous     
  wd <- (theta[4]-0)^2/vc[4,4] 
  
  cat('\nWald test of y2 weakly exogenous')
  cat('\nWald statistic           = ',wd) 
  cat('\np-value                  = ',1-pchisq(wd,1))
  cat('\n ')
  
  # Johansen estimator of reduced rank model
  johansen.model <- johansen(y,p,r,model)
  alpha <- johansen.model$alpha
  beta <- johansen.model$beta
  logl <- johansen.model$logl
  maxt <- johansen.model$maxt
  tracet <- johansen.model$tracet
  
  cat('\nJohansen procedure estimates')
  cat('\n----------------------------')
  cat('\nLog likelihood function = ',logl)
  cat('\nbeta_c         = ',-beta[3])
  cat('\nbeta_r         = ',-beta[2])
  cat('\nalpha_1        = ',alpha[1])
  cat('\nalpha_2        = ',alpha[2])
  
  
  # Zero rank case: VAR in first differences         
  # x  <- rep(1, nrow(yield)-1)
  y  <- trimr(yield,1,0) - trimr(yield,0,1)
  # b  <- lm(y ~ x - 1)$coef
  v  <- y             
  vc <- t(v) %*% v/nrow(v)
  n  <- ncol(y)
  lf <- -0.5*n*(log(2*pi) + 1) - 0.5*log(det(vc))
  
  cat('\n\nLog-likelihood of zero rank model    = ',lf)
  cat('\nDeterminant of covariance matrix     = ',det(vc))
  cat('\nCovariance matrix of residuals\n') 
  print(vc)
}




