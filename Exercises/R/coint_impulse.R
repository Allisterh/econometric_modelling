#============================================================================
#
#   Program to estimate impulse responses for a VECM and a VAR
#   of money demand
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr
source("EMTSUtil.R")
# Load required library - flipud, eye, ones, zeros
library(matlab)

#----------------------------------------------------------------------------
#  Johansen procedure 
#----------------------------------------------------------------------------
johansenH <- function(y,p,r,model,H) {
  ty <- nrow(y)
  dy <- trimr(y,1,0)-trimr(y,0,1) 
  z0 <- trimr(dy,p-1,0)
  z1 <- trimr(y,p-1,1) %*% H
  
  z2 <- c()
  for (j in seq(p-1)) {
    z2 <- cbind(z2, trimr(dy,p-1-j,j) )  
  }
  
  if (model == 1) {
    z1 <- z1
    z2 <- z2
  }else if(model == 2) {
    z1 <- cbind(trimr(y,p-1,1), rep(1, ty-p))
    z2 <- z2
  }else if (model == 3 ) {
    z1 <- z1
    z2 <- cbind(rep(1, ty-p),  z2)
  }else if (model == 4) {
    z1 <- cbind(z1,  seqa(1,1,ty-p))
    z2 <- cbind(rep(1, ty-p),  z2)
  } else if (model == 5) {
    z1 <- z1
    z2 <- cbind(rep(1, ty-p),  seqa(1,1,ty-p),  z2)
  }
  if (p == 1 && model <= 2) {
    r0 <- z0 
    r1 <- z1
  }else {
    r0 <- lm(z0 ~ z2 -1)$residuals
    r1 <- lm(z1 ~ z2 -1 )$residuals
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
  
  
  
  return(list(alpha=alpha,beta=beta,param=param,lam=tmp, logl=logl,maxt=maxt,tracet=tracet))
}

#----------------------------------------------------------------------------
#  C contains shocks (eg chol of variance matrix of v(t)
#  n is dimension of y
#----------------------------------------------------------------------------
impulse <- function (A,C,h,n) {
  rA <- nrow(A)
  cA <- ncol(A)
  
  rC <- nrow(C)
  cC <- ncol(C)
  
  p <- cA/rA
  
  #  VAR(1) companion form
  A <- rbind(A,
             cbind(diag((p-1)*n), array(0, c((p-1)*n,n))) )
  
  C <- rbind(cbind(C,  array(0, c(rC,(p-1)*cC))), array(0, c((p-1)*cC,p*cC)) )
  
  # Impulse responses  
  AC <- C
  
  tmp <- AC[1:n,1:n]
  imp <- as.numeric(tmp)
  
  for (j in seq(h)) {    
    AC  <- A %*% AC
    
    tmp <- AC[1:n,1:n]
    imp <- rbind(imp, c(tmp))    
  }  
  
  return(imp)
}


# 
#---------------------------------- Demand for Money ------------------------
#
coint_impulse <- function( ) {
  # Load the data
  load('moneydemand.Rdata')
  data <- as.matrix(money)
  cpi <- data[,1]
  fedfunds <- data[,2]
  gdp <- data[,3]
  m2 <- data[,4]
  tbill <- data[,5]
  
  lrm <- log(m2/cpi)
  lry <- log(gdp/cpi)
  ym  <- cbind(lrm,  lry,  tbill/100,  fedfunds/100)
  
  
  t <- 188
  y <- ym[1:t,]
  
  p     <- 4      # Lags in VAR                    
  r     <- 1      # Number of cointegrating equations   
  model <- 3
  
  # Impose unit coefficient on y
  H <- matrix(c(1, 0, 0,
                -1, 0, 0,
                0, 1, 0,
                0, 0, 1), nrow=4, byrow=T)
  
  johansenH.model <- johansenH(y,p,r,model,H)
  lam <- johansenH.model$lam
  beta <- johansenH.model$beta
  
  betahat <- H %*% beta 
  
  #Estimate remaining VECM coefficients, model 3
  dy <- trimr(y,1,0)-trimr(y,0,1) 
  x <- cbind(rep(1, t-p), trimr(y,p-1,1) %*% betahat)
  
  for (j in seq(p-1)) {
    x <- cbind(x, trimr(dy,p-1-j,j))
  }
  
  bhat <- lm(trimr(dy,p-1,0) ~ x - 1)$coef
  vhat <- trimr(dy,p-1,0)-x %*% bhat
  
  ghat <- trimr(bhat,r+1,0)
  
  c    <- ncol(y)
  gam1 <- eye(c)
  
  for (j in seq(p-1)) {
    gam1 <- gam1 - t(ghat[((j-1)*c+1):(j*c),])
  }
  
  
  # Put VECM in levels VAR form  
  alphahat <- t(bhat[2:(r+1),])
  ghat     <- cbind(t(trimr(bhat,r+1,0)), array(0, c(c,c)))
  A        <- eye(c)+t(alphahat) %*% t(betahat)+ghat[,1:c]
  for (j in 2:p) {
    tmp1 <- (j-1)*c
    tmp2 <- (j-2)*c
    A <- cbind(A, ghat[,(tmp1+1):(j*c)]-ghat[,(tmp2+1):tmp1] )
  }
  
  vecm_impulse <- impulse(A,t(chol(t(vhat) %*% vhat/nrow(vhat))), 40, c )
  
  #  Estimate levels VAR for comparison with VECM
  x <- ones(t-p,1)
  
  for (j in seq(p)) {
    x <- cbind(x,  trimr(y,p-j,j))  
  }
  
  bhat <- lm(trimr(y,p,0) ~ x - 1)$coef
  vhat <- trimr(y,p,0)-x %*% bhat
  
  var_impulse <- impulse(t(trimr(bhat,1,0)),t( chol(t(vhat) %*% vhat/nrow(vhat))),40,c)
  
  # VECM impulse responses to interpret coefficient on income
  k <- rep(0, 4)
  k <- cbind(k, zeros(4,3))
  
  k[2,2] <- 1
  k[1,2] <- 1
  
  vecm1_impulse <- impulse(A,gam1 %*% k,40,c)
  
  hh <- 0:40
  cat('\n   Illustrating the unit long run elasticity of money wrt income' )
  cat('\n----------------------------------------------------------------\n')
  time <- hh
  money <- vecm1_impulse[,5]
  income <- vecm1_impulse[,6]
  tbill <- vecm1_impulse[,7]
  fedfunds <- vecm1_impulse[,8]
  
  print(cbind(hh,  money, income, tbill, fedfunds ))
  
  
  figure()
  par(mfrow=c(3,1))
  
  
  matplot(hh,var_impulse[,5:8], type="l",
          main = 'Impulses associated with an income shock: VAR in levels',
          xlab = '',
          ylab = '',
          bty = 'l')
  matplot(hh,vecm_impulse[,5:8], type='l',
          main = 'Impulses associated with an income shock: VECM',
          xlab = '',
          ylab = '',
          bty = 'l')
  
  matplot(hh,vecm1_impulse[,5:8], type='l',
          main = 'Impulses associated with an income shock: VECM with unit restriction',
          xlab = '',
          ylab = '',
          bty = 'l')
  
}

