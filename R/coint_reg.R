#============================================================================
#
#   Program to compare single equation cointegration estimators
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123)

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr,recserar
source("EMTSUtil.R")
# Load required library - flipud, ones, zeros
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
  
  for (j in seq(p-1)) {    
    if(j > 0)      
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
  
  # Statistics
  logl   <- -0.5*tc*(log(2*pi)+1) - 0.5*log(det(s00)) - 0.5*sum( log(1-tmp[1:r]))
  tracet <- -tr*flipud(cumsum(log(1-flipud(tmp))))                     
  maxt   <- -tr*log(1-tmp) 
  return(list(alpha=alpha,beta=beta, logl=logl,maxt=maxt,tracet=tracet))
}

#----------------------------------------------------------------------------
#  Long-run variance
#----------------------------------------------------------------------------
lrvars <- function(z) {
  t <- nrow(z)
  k <- ncol(z)
  p     <- floor(8*(t/100)^(2/9))
  
  
  J0 <- t(z) %*% z 
  J1 <- 0
  
  for (j in seq(p)) {
    Gj <- t(trimr(z,j,0)) %*% trimr(z,0,j)
    J0 <- J0 + Gj + t(Gj)
    J1 <- J1 + 2*j*Gj
  }
  
  i <- ones(k,1)
  v0 <- t(i) %*% J0 %*% i 
  v1 <- t(i) %*% J1 %*% i
  p <- min( rbind(floor(1.1447*((v1/v0)^2*t)^(1/3)), p ) )
  
  Omega <- t(z) %*% z/t 
  Delta <- t(z) %*% z/t
  
  for (j in seq(p)) {
    Gj    <- t(trimr(z,j,0)) %*% trimr(z,0,j)/t
    Omega <- Omega + (1-j/(p+1))*(Gj + t(Gj))
    Delta <- Delta + (1-j/(p+1))*Gj
  }
  return(list(Delta=Delta, Omega=Omega))
}
#-------------------------------------------------------------------------
#  Johansen VECM estimator
#-------------------------------------------------------------------------
vecm <- function(y) {
  t <- nrow(y)
  k <- ncol(y)
  
  pmax  <- 12 
  
  logL <- c() 
  x    <- c()
  y0   <- trimr(y,pmax,0) 
  
  tmp  <- nrow(y0)
  logL <- -0.5*tmp*log(det(t(y0) %*% y0/tmp))
  
  
  
  for (j in seq(pmax)) {
    x    <- cbind(x, trimr(y,pmax-j,j))
    vhat <- lm(y0 ~ x - 1)$residuals    
    logL <- rbind(logL, -0.5*nrow(vhat)*log(det(t(vhat) %*% vhat/nrow(vhat))) )
  }
  nparams <- k+k^2*seqa(0,1,pmax+1)
  HQ      <- -2*logL/tmp+2*nparams*log(log(tmp))/tmp
  p  <- which.min(HQ)
  p       <- p-1
  
  beta <- johansen(y,p,1,1)$beta
  
  bhat <- -beta[2]
  
  return(bhat)
}


# 
#---------------------------------- Cointegrating regression ------------------------
#             
coint_reg <- function( ) {
  reps <- 10000
  Tv   <- c(100,200,400)
  rhov <- c(0,0.4,0.8)
  av   <- c(0,0.4,0.8)
  
  bhat  <- zeros(reps,3) 
  
  for (ac in seq(av)) {
    for (rhoc in seq(rhov)) {
      for (tc in seq(Tv) ) {
        t <- Tv[tc]
        for (rep in seq(reps)) {
          u      <- matrix(rnorm(t*2), t,2)
          u[,2] <- rhov[rhoc]*u[,1]+sqrt(1-rhov[rhoc]^2)*u[,2] 
          u[,1] <- recserar(cbind(u[,1]),cbind(u[1,1]),cbind(av[ac]))
          u[,1] <- u[,1]/sd(u[,1])
          
          x <- cumsum(u[,2])
          y <- x + u[,1]        
          
          # OLS
          bhat[rep,1] <- lm(y ~ x - 1)$coef
          e           <- y-x * bhat[rep,1]
          
          # FMOLS
          dx <- trimr(x,1,0)-trimr(x,0,1)
          lvars.fn <- lrvars( cbind(trimr(e,1,0), dx))
          Omega <- lvars.fn$Omega
          Delta <- lvars.fn$Delta
          ys            <- trimr(y,1,0) - dx/Omega[2,2]*Omega[2,1]
          Deltas        <- Delta[1,2]-Delta[2,2]/Omega[2,2]*Omega[2,1]
          
          bhat[rep,2]   <- (ys %*% trimr(x,1,0)-t*Deltas)/(trimr(x,1,0) %*% trimr(x,1,0))
          
          
          # VECM
          bhat[rep,3] <- vecm(cbind(y, x))        
          
        }
        cat('\n\nSample size       = ',t) 
        cat('\nalpha             = ',av[ac])
        cat('\nrho               = ',rhov[rhoc])
        
        cat('\n            Bias                            Std Dev           ')
        cat('\n   --------------------------------    ----------------------------------------')
        cat('\n      OLS        FMOLS         VECM         OLS       FMOLS        VECM\n')
        
        cat(colMeans(bhat - 1), apply(bhat, 2, sd) )
      }
    }
  }
  
  
}








