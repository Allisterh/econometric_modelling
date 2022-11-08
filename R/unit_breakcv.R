#============================================================================
#
#   5% cvs for unit root tests with break in trend at known point
#
#============================================================================
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, inv, recserar
source("EMTSUtil.R")

#------------------------------------------------------------------------- --
# Detrending function: 
# cbar=-7 constant; 
# cbar=-13.5 linear trend 
# cbar=-T for OLS detrending 
#----------------------------------------------------------------------------
glsdetrend <- function(y,x,cbar)  {  
  t <- length(y) 
  yc <- matrix(c(y[1],(trimr(y,1,0)-(1 + cbar/ t) * trimr(y,0,1))), nrow=t)
  xc <- matrix(c(x[1, ],(trimr(x,1,0)-(1 + cbar/ t) * trimr(x,0,1))), nrow=t )
  b <- lm(yc ~ xc - 1)$coef
  u <- y- x %*% b 
  
  return(u)
} 

#---------------------------------------------------------------------------- 
# ADF coefficient and t tests:u must be already de-trended 
#---------------------------------------------------------------------------- 
adftests <- function(u,k)  {
  
  du <- cbind(trimr(u,1,0)-trimr(u,0,1) )
  x <- trimr(u,0,1) 
  
  # Set up lags of du 
  if (k > 0)  { 
    ldu <- lag.matrix(du,seq(from = 1,to = k,by = 1)) 
    x <- trimr(cbind(x,ldu),k,0) 
  } 
  
  xxi <- inv( t(x) %*% x ) 
  b <- xxi %*% t(x) %*% trimr(du,k,0) 
  e <- trimr(du,k,0)-x %*% b 
  s2 <- t(e) %*% e/ length(e)
  
  adfstats <- rbind(length(u)*b[1],b[1]/ sqrt(s2*xxi[1,1]))
  
  return(adfstats)
}


#---------------------------------------------------------------------------- 
# M tests 
#---------------------------------------------------------------------------- 
mtests <- function(u,k)  {  
  s2 <- ar1rvar(u,k) 
  n <- length(u) 
  tmp <- sum(u[1:(n-1)]^2)   
  
  u2 <- tmp/n^2 
  mza <- (u[n]^2/n-s2)/ (2*u2)
  msb <- sqrt(u2/ (s2))
  mzt <- msb %*% mza 
  
  tests <- cbind(mza,msb,mzt) 
  
  return(tests)
}
#---------------------------------------------------------------------------- 
# Autoregressive long run variance estimator 
#---------------------------------------------------------------------------- 
ar1rvar <- function(u,k)  {
  
  du <- cbind(trimr(u,1,0)-trimr(u,0,1) )
  x <- trimr(u,0,1) 
  
  if (k > 0)  { 
    x <- cbind(x,lag.matrix(du,seqa(1,1,k))) 
    x <- cbind(x[-(1:k), ])
  }
  x <- cbind(x)  
  b <- lm(trimr(du,k,0) ~ x - 1)$coef
  
  e <- trimr(du,k,0)- x %*% b 
  s2 <- t(e) %*% e/ length(e)  
  if (k > 0)  { 
    s2 <- s2/ (1-sum(trimr(b,1,0)))^2
  }   
  return(s2)
}


unit_breakcv <- function() {
  reps <- 1000
  
  Tv   <- c(25,50,100,250,500,1000)
  tauv <- seqa(0.15,0.05,15)
  cbar <- c(-17.6,-17.8,-18.2,-18.4,-18.6,-18.4,-18.4,-18.2,-18.0,-17.6,-17.4,-17.0,-16.6,-16.0,-15.2)
  
  # Allocate arrays
  dfols    <- array(0, c(reps,length(tauv)))
  dfgls    <- dfols 
  mztols   <- dfols
  dfolscv  <- array(0, c(length(Tv),length(tauv)))
  dfglscv  <- dfolscv 
  mztolscv <- dfolscv
  
  for (Tc in 1:length(Tv)) {
    t <- Tv[Tc]
    cat('\nCurrent sample size = ',t)
    
    for (tauc in 1:length(tauv)) {
      TB <- floor(tauv[tauc]*t)
      DT <- cbind(c(rep(0, TB), seqa(1,1,t-TB)))
      x  <- cbind(rep(1, t), seqa(1,1,t), DT)
    
      set.seed(42)
      
      for (rep in seq(reps)) {
        y <- cumsum(rnorm(t))       
        uols <- glsdetrend(y,x,-t)
        ugls <- glsdetrend(y,x,cbar[Tc])
        dfols[rep,tauc] <- trimr(adftests(uols,0),1,0)
        dfgls[rep,tauc] <- trimr(adftests(ugls,0),1,0)
        mztols[rep,tauc] <- trimr(t(mtests(uols,0)),2,0)
      }  
      
    }    
    dfolscv[Tc,] <- quantile(dfols,0.05)
    dfglscv[Tc,] <- quantile(dfgls,0.05)
    mztolscv[Tc,] <- quantile(mztols,0.05)    
  }
                          
  cat('\nDF-OLS')
  cat('\n    tauv\n')
  print(cbind(Tv, dfolscv))
  cat('\n')
  cat('\nDF-GLS')
  cat('\n   tauv\n')
  print(cbind(Tv, dfglscv))
  
  cat('\n')
  
  cat('\nMZt-OLS')
  cat('\n    tauv')
  print(cbind(Tv, mztolscv))
}