#============================================================================
#
#   Program to estimate a bivariate vecm of the 
#    term structure of interest rates with level effects on the variance
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Reduced rank log-likelihood function: unrestricted
#----------------------------------------------------------------------------
neglog <- function(b,z0,z1) {
  nobs <- nrow(z0)
  f   <- rep(0, nobs)
  
  v1 <- z0[,1] - b[3]*(z1[,1] - b[1] - b[2]*z1[,2])
  v2 <- z0[,2] - b[4]*(z1[,1] - b[1] - b[2]*z1[,2])
  v  <- cbind(v1,  v2)  
  
  sd  <-  matrix(c(b[5],   0,
                   b[6],  b[7]), nrow=2, byrow=T)
  
  # Level effects parameters     
  kappa  <-  c(b[8], b[9])       
  
  n <- ncol(v)
  
  for (t in seq(nobs)) {
    #  Specify level-effects covariance structure using Choleski 
    l <- sd * (z1[t,]^kappa)
    omegav <- l %*% t(l)    
    f[t] <- - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[t,] %*% inv(omegav) %*% cbind(v[t,])
  }  
  lf <- -mean(f)    
  return(lf)
}


#
#--------------- Term structure of Interest Rates with Level effects --------
#

coint_level <- function() {
  # Load data and choose 10-year and 1-year yields
  data <- as.matrix(read.table("usmacro.dat"))  
  yield <- cbind(data[,3], data[,1])
  
  z0   <- trimr(yield,1,0) - trimr(yield,0,1)
  z1   <- trimr(yield,0,1)
  
  t    <- nrow(z0)
  nobs <- nrow(yield)
  
  # Starting values
  b <- lm(yield[,1] ~ cbind(rep(1, nobs), yield[,2]) - 1)$coef
  u <- yield[,1] - cbind(rep(1, nobs), yield[,2]) %*% b
  
  a <- lm(z0 ~ trimr(u,0,1) - 1)$coef
  v <- z0 - trimr(u,0,1) %*% a
  
  
  # Estimate levels effect model
  start <- c(b, a, vech(t(chol(cov(v)))),0.0,0.0)
  estResults <- optim(start, neglog, z0=z0, z1=z1, method="BFGS", hessian=T)
  
  theta <- estResults$par
  logl <- estResults$val
  hess <- estResults$hessian
  
  cat('\nLevels effect model estimates')        
  cat('\n-----------------------------')
  cat('\nLog likelihood function = ',-logl)
  cat('\nbeta_c         = ',theta[1])
  cat('\nbeta_r         = ',theta[2])
  cat('\nalpha_1        = ',theta[3])
  cat('\nalpha_2        = ',theta[4])
  
  
  vc <- (1/t)*inv(hess)
  
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
  
  
  # Wald test of levels effect
  r  <- matrix(c(0,  0,  0,  0,  0,  0,  0,  1,  0, 
                 0,  0,  0,  0,  0,  0,  0,  0,  1), nrow=2, byrow=T)
  q  <- cbind(c(0,0))
  wd <- t((r %*% theta - q)) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta - q)
  
  cat('\nWald test of levels effect')
  cat('\nWald statistic           = ',wd)
  cat('\np-value                  = ',1-pchisq(wd,2))
  cat('\n ')  
}


