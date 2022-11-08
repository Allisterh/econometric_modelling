#===========================================================================
#
#   Program to estimate trivariate term structure model
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - rref, trimr
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
  return(list(alpha=alpha,beta=beta,logl=logl,maxt=maxt,tracet=tracet, lam=tmp))
}

# 
#-------------------- Term Structure of Interest Rates ----------------------
#
coint_triterm <- function() {
  # Load data and choose 10-year and 1-year yields
  data <- as.matrix(read.table('usmacro.dat')) 
  y     <- cbind(data[,3], data[,2], data[,1])
  t     <- nrow(y)
  
  # Reduced rank case parameters
  p     <- 1      
  r     <- 2    
  model <- 2
  
  
  # Johansen estimator of reduced rank model
  johansen.model <- johansen(y,p,r,model)
  alpha <- johansen.model$alpha
  beta <- johansen.model$beta
  logl <- johansen.model$logl
  maxt <- johansen.model$maxt
  tracet <- johansen.model$tracet
  lam <- johansen.model$lam
  
  
  cat('\nJohansen procedure estimates')
  cat('\n----------------------------')
  cat('\nLog likelihood function = ',logl)
  
  cat('\nEstimates of beta in row eschelon form')
  cat('\n--------------------------------------')
  cat('\nbeta\n')
  print(t(beta))
  cat('\n ')
  
  cv = rbind(34.938, 20.205,  9.142)
  
  cat('\n')
  print(cbind(Eigs=trimr(lam,0,1), "LR Test"=trimr(tracet,0,1), cv=cv))
  
}
