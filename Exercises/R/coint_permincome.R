
#============================================================================
#
#   Program to estimate a bivariate vecm of the permanent income hypothesis
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - rref, trimr
source("EMTSUtil.R")
# Load required library - flipud, Null
library(matlab)
library(MASS)

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
  
  tmpx <- cbind(z1 %*% beta, z2)
  param <- lm(z0 ~ tmpx - 1)$coef
  
  # Statistics
  logl   <- -0.5*tc*(log(2*pi)+1) - 0.5*log(det(s00)) - 0.5*sum( log(1-tmp[1:r]))
  tracet <- -tr*flipud(cumsum(log(1-flipud(tmp))))                     
  maxt   <- -tr*log(1-tmp) 
  return(list(alpha=alpha,beta=beta,param=param, logl=logl,maxt=maxt,tracet=tracet))
}


# 
#----------------------- Permanent Income hypothesis ------------------------
#
coint_permincome <- function( ) {
  # Load data 
  data <- as.matrix(read.table("permincome.dat"))
  
  
  # Select desired sample (1984Q1 to 2005Q4)
  rc <- data[149:236,1]
  ry <- data[149:236, 2]
  y  <- cbind(log(rc),  log(ry))
  
  # Reduced rank case: Johansen estimator       
  p     <- 4     
  r     <- 1     
  model <- 3  
  
  johansen.model <- johansen(y,p,r,model)
  alpha <- johansen.model$alpha
  beta <- johansen.model$beta
  param <- johansen.model$param
  
  cat('\nEstimates of cointegrating vector   = ', beta)
  
  mu0 <- t(param[r+1,])
  cat('\nEstimates of alpha                  = ', alpha)
  cat('\nEstimates of mu0                    = ', mu0)
  
  # Compute orthonormal complement of alpha
  alpha0 <- Null(alpha)      
  cat('\nOrthogonal complement matrix alpha0 = ', alpha0)
  
  # Compute constant term in cointegrating regression
  beta0 <- lm(mu0 %*% alpha ~ t(alpha) %*% alpha - 1)$coef
  
  cat('\nEstimate of beta0                   = ', beta0)
  
}
    