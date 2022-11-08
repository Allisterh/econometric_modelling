#============================================================================
#
#   Likelihood ratio tests of deterministic trends in a vecm
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234)

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr
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
  
  tmpx <- cbind(z1 %*% beta, z2)
  param <- lm(z0 ~ tmpx - 1)$coef
  
  # Statistics
  logl   <- -0.5*tc*(log(2*pi)+1) - 0.5*log(det(s00)) - 0.5*sum( log(1-tmp[1:r]))
  tracet <- -tr*flipud(cumsum(log(1-flipud(tmp))))                     
  maxt   <- -tr*log(1-tmp) 
  return(list(alpha=alpha,beta=beta,param=param, logl=logl,maxt=maxt,tracet=tracet))
}

# 
#---------------------- Likelihood ratio test time trends -------------------
#
coint_trend <- function() {
  t <- 200 
  n <- 2
  
  # Parameters: True DGP based on Model 3
  beta   <- c(1,-1)
  beta0  <- 4 
  beta1  <- 0
  alpha  <- c(-0.2,0.2)
  alpha0 <- c(1,1)
  delta0 <- alpha0*0.2 
  delta1 <- alpha0*0
  psi1   <- matrix(c(-.2, 0, 
                     0, .4), nrow=2, byrow=T)
  
  # Simulate the model  
  v <- matrix(rnorm(t*n), nrow=t, ncol=n)
  y <- array(0, c(t+2,2))
  
  for (j in (3:t+2)) {
    u <- beta0 + beta1*j + y[j-1,]*beta        
    y[j,] <- y[j-1,] + t(delta0 + delta1*j + alpha*u) + (y[j-1,] - y[j-2,]) %*% psi1 + v[j-2,]
  }
  
  y   <- trimr(y,2,0)
  mu0 <- alpha*beta0 + delta0
  mu1 <- alpha*beta1 + delta1
  
  # Estimate the vecm by maximum likelihood using Johansen estimator 
  p <- 2      # Number of lags in VAR      
  r <- 1      # Number of cointegrating equations    
  
  # Estimate Model 5    
  lnl5 <- johansen(y,p,r,5)$logl  
  
  # Estimate Model 4    
  lnl4 <- johansen(y,p,r,4)$logl 
  
  # Estimate Model 3     
  lnl3 <- johansen(y,p,r,3)$logl     
  
  
  # LR test of Model 5 (unrestricted) against Model 4 (restricted)    
  lr <- -2*t*(lnl4 - lnl5)
  
  cat('\nLog-likelihood (unrestricted model 5) = ',lnl5)
  cat('\nLog-likelihood (restricted model 4)   = ',lnl4)
  cat('\nValue of likelihood ratio statistic   = ',lr)
  cat('\np-value                               = ',1-pchisq(lr,n-r))
  cat('\n')
  
  # LR test of Model 4 (unrestricted) against Model 3 (restricted)    
  lr <- -2*t*(lnl3 - lnl4)
  
  cat('\nLog-likelihood (unrestricted model 4) = ',lnl4)
  cat('\nLog-likelihood (restricted model 3)   = ',lnl3)
  cat('\nValue of likelihood ratio statistic   = ',lr)
  cat('\np-value                               = ',1-pchisq(lr,1))
  cat('\n')
  
}
