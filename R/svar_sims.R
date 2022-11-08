#============================================================================
#
#   Estimate the Sims 6-variate SVAR model with long-run restrictions
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, inv, diagrv, lag.matrix
source("EMTSUtil.R")

#
#-------------------------- Sims Model --------------------------------------
#

svar_sims <- function () {
  
  ytdata <- as.matrix(read.table("sims_data.dat"))
  # Define varaibles
  r    <- ytdata[,1]
  lex  <- log( ytdata[,2] )
  lcp  <- log( ytdata[,3] )
  lm   <- log( ytdata[,4])
  lp   <- log( ytdata[,5])
  lo   <- log( ytdata[,6])
  sdum <- ytdata[,7:17]
  y    <- cbind(r , lex , lcp , lm , lp , lo)
  
  # Estimate the VAR with p lags  
  p <- 14  	
  ylag <- lag.matrix(y,1:p)
  ylag <- ylag[-(1:p),]
  
  xmat <- cbind(rep(1, nrow(ylag)),  trimr(sdum,0,p), ylag)
  ymat <- trimr(y,p,0)
  
  bar <- lm(ymat ~ xmat - 1)$coef
  v   <- ymat - xmat %*% bar
  vc  <- t(v) %*% v/nrow(v)
  
  #		Reduced form approach	
  u1 <- v[,1]
  
  reg <- lm(v[,2] ~ v[,1] - 1)
  a2 <- reg$coef
  u2 <- reg$residuals
  
  reg <- lm(v[,3] ~ v[,1:2] - 1)
  a3 <- reg$coef
  u3 <- reg$residuals
  
  
  reg <- lm(v[,4] ~ v[,1:3] - 1)
  a4 <- reg$coef
  u4 <- reg$residuals
  
  
  reg <- lm(v[,5] ~ v[,1:4] - 1)
  a5 <- reg$coef
  u5 <- reg$residuals
  
  reg <- lm(v[,6] ~ v[,1:5] - 1)
  a6 <- reg$coef
  u6 <- reg$residuals
  
  b0_rf <- matrix(c(1   ,    0     ,    0    ,    0    ,    0    ,   0 ,
                -a2[1]  ,    1     ,    0    ,    0    ,    0    ,   0,
                -a3[1]  , -a3[2]   ,    1    ,    0    ,    0    ,   0,
                -a4[1]  , -a4[2]   , -a4[3]  ,    1    ,    0    ,   0,
                -a5[1]  , -a5[2]   , -a5[3]  , -a5[4]  ,    0    ,   0,
                -a6[1]  , -a6[2]   , -a6[3]  , -a6[4]  , -a6[5]  ,   1), nrow=6, byrow=T)   
  u     <- cbind(u1 , u2 , u3 , u4 , u5 , u6)
  d_rf  <- t(u) %*% u/nrow(u)
  
  # Choleski decomposition approach	
  s      <- t(chol(vc))
  b0inv  <- t(apply(s, 1, '/', diag(s)))
  b0_cd <- inv(b0inv)
  d_cd <- diagrv(diag(6),diag(s)^2) 
  
  cat('\nB0: Reduced form\n')
  print(unname(b0_rf))
  
  cat('\nB0: Choleski decomposition\n')
print(unname(b0_cd ))
}



