#============================================================================
#
#   Program to compare the performance of alternative tests of cointegration
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234)
# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr,recserar
source("EMTSUtil.R")
library(matlab)

reps <- 10000
t    <- 100                                 
n    <- 2        
r    <- seqa(0,1,n+1)                                                     
m    <- n*r + (n-r)*r+ n*(n+1)/2  # Number of free parameters                               
trcv <- c(12.285, 4.173)            # 5% cvs for common trends of n-r = 2,1  

# Parameters
phi1 <- 0.8                             
phi2 <- 0.8
rho  <- 0.0
phi  <- diag(c(phi1, phi2) )

omegav      <- matrix(c(1.0,  rho,
                        rho,  1.0), nrow=2, byrow=T)
chol_omegav <- chol(omegav)

# Initialise arrays
rtr <- rep(0, reps)
rbic <- rep(0, reps)
raic <- rep(0, reps)
rhic <- rep(0,reps)

for (i in seq(reps)) {
  u  <- matrix(rnorm(t*n), nrow=t, ncol=n) %*% chol_omegav                                        
  y  <- recserar(u,rbind(rep(0,n)), rbind(diag(phi)))
  
  # Construct variables for VECM
  dy <- trimr(y,1,0) - trimr(y,0,1)   
  y1 <- trimr(y,0,1)
  r0 <- dy 
  r1 <- y1
  
  # Perform eigen decomposition      
  tmp <- nrow(r1)
  s11 <- t(r1) %*% r1/tmp 
  s10 <- t(r1) %*% r0/tmp 
  s00 <- t(r0) %*% r0/tmp
  
  l <- t(chol(s11))
  lam  <- eigen( inv(l) %*% s10 %*% inv(s00) %*% t(s10) %*% inv(t(l)) )$values
  lam <- flipud(cbind(lam))
  
  k <- nrow(r0)
  c <- ncol(r0)
  lnl <- -c*k*0.5*(log(2*pi) + 1)  - 0.5*k*log(det(s00)) - 0.5*k*(cbind(c(0,cumsum(log(1 - lam))) )) 
  
  
  #  Compute test statistics
  tr <- -2*( trimr(lnl,0,1) - lnl[n+1] )     
  
  # LR Trace
  ind <- which.max(rbind(tr < trcv, 1))
  rtr[i]  <- r[ind]
  
  # AIC
  ind <- which.min(-2*lnl + 2*m)
  raic[i] <- r[ind]
  
  # BIC
  ind <- which.min(-2*lnl + log(t)*m)
  rbic[i] <- r[ind]
  
  # HIC
  ind <- which.min(-2*lnl + 2*log(log(t))*m)
  rhic[i] <- r[ind]
}

# Collect results
zeros <- array(0, c(reps, n+1))
trace <- zeros
aic <- zeros
bic <- zeros
hic <- zeros

for (i in 1:(n+1)) {  
  trace[,i] <- rtr == i-1
  aic[,i]   <- raic == i-1
  bic[,i]   <- rbic == i-1
  hic[,i]   <- rhic == i-1  
}
trace <- colMeans(trace)
aic   <- colMeans(aic)
bic   <- colMeans(bic)
hic   <- colMeans(hic)
print(cbind(Rank=r, "LR(trace)"=trace, AIC=aic, BIC=bic, HIC=hic))

