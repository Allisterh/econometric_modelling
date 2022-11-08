#============================================================================
#
#   Program to test for weak exogeneity using Wald and LR tests
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1)

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Reduced rank log-likelihood function (unrestricted)
#----------------------------------------------------------------------------
neglog1 <- function(b,z0,z1,z2) {
  nobs <- nrow(z0)
  f   <- rep(0, nobs)
  
  m1  <- b[1]*(z1[,1] - b[3]*z1[,2]) + z2 %*% b[c(4, 5)]
  m2  <- b[2]*(z1[,1] - b[3]*z1[,2]) + z2 %*% b[c(6, 7)]
  
  v1 <- z0[,1] - m1
  v2 <- z0[,2] - m2
  v  <- cbind(v1,  v2)
  
  k <- nrow(v)
  n <- ncol(v)
  
  omegav <- t(v) %*% v/k
  for (t in seq(nobs)) {
    f[t] <- - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[t,] %*% inv(omegav) %*% cbind(v[t,])
  }  
  lf <- - mean( f )  
  return(lf)
}


#----------------------------------------------------------------------------
#  Reduced rank log-likelihood function (restricted: y2 weakly exogenous)
#----------------------------------------------------------------------------
neglog0 <- function(b,z0,z1,z2) {
  nobs <- nrow(z0)
  f   <- rep(0, nobs)
  
  m1  <- b[1]*(z1[,1] - b[2]*z1[,2]) + z2 %*% b[c(3, 4)]
  m2  <- z2 %*% b[c(5, 6)]
  
  v1 <- z0[,1] - m1
  v2 <- z0[,2] - m2
  v  <- cbind(v1,  v2)
  
  k <- nrow(v)
  n <- ncol(v)
  
  omegav <- t(v) %*% v/k
  
  for (t in seq(nobs)) {
    f[t] <- - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[t,] %*% inv(omegav) %*% v[t,]
  }  
  lf <- - mean( f )  
  return(lf)
}


#----------------------------------------------------------------------------
#  Reduced rank log-likelihood function (restricted: y2 weakly exogenous)
#----------------------------------------------------------------------------
neglog00 <- function(b,z0,z1,z2) {
  nobs <- nrow(z0)
  f   <- rep(0, nobs)
  
  m1  <- b[1]*(z1[,1] - b[2]*z1[,2]) + z2 %*% b[c(3, 4)]
  m2  <- z2 %*% rbind(b[5], 0)
  
  v1 <- z0[,1] - m1
  v2 <- z0[,2] - m2
  v  <- cbind(v1,  v2)
  
  k <- nrow(v)
  n <- ncol(v)  
  
  omegav <- t(v) %*% v/k
  
  for (t in seq(nobs)) {
    f[t] <- - n*0.5*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[t,] %*% inv(omegav) %*% v[t,]
  }  
  lf <- - mean( f )
  return(lf)
}


#----------------------------------------------------------------------------
#  Johansen procedure 
#----------------------------------------------------------------------------
johansen <- function(y,p,r,model) {
  ty <- nrow(y)
  dy <- trimr(y,1,0)-trimr(y,0,1) 
  z0 <- trimr(dy,p-1,0)
  z1 <- trimr(y,p-1,1)  
  z2 <- c()
  
  for (j in seq(p-1)) {
    z2 <- cbind(z2, trimr(dy,p-1-j,j))
  }
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
# ----------------------  Exogeneity Testing ---------------------------------
#
coint_exogeneity <- function() {
  t <- 200 
  n <- 2
  
  # Parameters: true DGP based on Model 3   
  beta      <- c(1, -1)
  beta0     <- 0 
  beta1     <- 0
  alpha     <- c(-0.4,0) 
  alphaorth <- c(0, 1 )
  delta0    <- alphaorth*0 
  delta1    <- alphaorth*0
  psi1      <- matrix(c(-0.2, 0,
                        0, 0.4), nrow=2, byrow=T)
  omegav    <- matrix(c(1, 0.5,
                        0.5, 1), nrow=2, byrow=T)
  
  # Simulate the model   
  v <- matrix(rnorm(t*2), nrow=t) %*% chol(omegav)
  y <- zeros(t+2,2)
  
  for (j in 3:(t+2 )) {
    u      <- beta0 + beta1*j+y[j-1,]*beta
    y[j,] <- y[j-1,] + t(delta0 + delta1*j + alpha %*% cbind(u)) + (y[j-1,] - y[j-2,]) %*% psi1 + v[j-2,]
  }
  y <- trimr(y,2,0)
  
  # Estimate the vecm by maximum likelihood  
  p     <- 2      # Number of lags in VAR      
  r     <- 1      # Number of cointegrating equations   
  model <- 1
  Tp <- t-p
  # Johansen estimator of reduced rank model
  beta <- johansen(y,p,r,model)$beta
  
  cat('\nJohansen estimate of long-run parameters\n')
  print(beta)
  cat('\n')
  
  
  dy <- trimr(y,1,0)-trimr(y,0,1) 
  z0 <- trimr(dy,p-1,0)
  z1 <- trimr(y,p-1,1)
  
  z2 <- c()
  for (j in 1:(p-1)) {
    z2 <- cbind(z2, trimr(dy,p-1-j,j))
  }
  
  nobs <- nrow(z0)
  
  # Estimate the unconstrained model
  theta_0                 <- 0.1*rep(1, 7)
  estResults <- optim(theta_0, neglog1, z0=z0, z1=z1, z2=z2, method="BFGS", hessian=T)    
  theta1 <- estResults$par
  lf1 <- estResults$val
  hess <- estResults$hessian
  
  
  lf1 <- -lf1
  vc  <- (1/nobs)*inv(hess)
  
  # Wald test (y2 weakly exogenous)     
  r <- cbind(0,   1,   0,   0,   0,   0,   0)
  q <- 0
  w    <- t(r %*% theta1 - q) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta1 - q)
  
  cat('\nWald test (y2 weakly exogenous)')
  cat('\nWald statistic           = ',w)
  cat('\np-value                  = ',1-pchisq(w,1))
  cat('\n ')
  
  
  # Wald test (y2 strongly exogenous)       
  r <- cbind(0,   1,   0,   0,  0,  1,   0)
  q <- 0
  w    <- t(r %*% theta1 - q) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta1 - q)
  
  
  cat('\nWald test (y2 strongly exogenous)')
  cat('\nWald statistic           = ',w)
  cat('\np-value                  = ',1-pchisq(w,2))
  cat('\n ')
  
  
  # Estimate constrained model (y2 is weakly exogenous)   
  theta_0             <- c(theta1[1], theta1[3:7])
  estResults <- optim(theta_0, neglog0, z0=z0, z1=z1, z2=z2, method="BFGS")
  lf0 <- estResults$val
  
  
  lf0 <- -lf0
  
  # LR test of weak exogeneity
  lr <- -2*nobs*(lf0 - lf1)
  
  cat('\nLR test (y2 weakly exogenous)')
  cat('\nLR statistic           = ',lr)
  cat('\np-value                = ',1-pchisq(lr,1))
  cat('\n ')
  
  # Estimate constrained model (y2 is strongly exogenous)   
  theta_0         <- c(theta1[1], theta1[3:5], theta1[7])
  estResults <- optim(theta_0, neglog00, z0=z0, z1=z1, z2=z2, method="BFGS")
  lf00 <- estResults$val
  
  lf00 <- -lf00
  lr <- -2*(lf00 - lf1)
  
  cat('\nLR test (y2 strongly exogenous)')
  cat('\nLR statistic           = ',lr)
  cat('\np-value                = ',1-pchisq(lr,2))
  cat('\n ')
  
  # Estimate partial model just consisting of an augmented equation 
  # for y1 by that assuming y2 is weakly exogenenous     
  
  bhat <- lm(z0[,1] ~ cbind(z1, z2, z0[,2]) - 1)$coef
  
  cat('\nEstimate of long-run parameter based on partial model\n')
  print(unname(-bhat[2]/bhat[1]))
}




