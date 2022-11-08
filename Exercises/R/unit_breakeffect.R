#=========================================================================
#
#   The effect of a neglected break in trend on unit root tests 
#   DF and MZt tests with OLS and GLS de-trending
#
#=========================================================================
rm(list = ls(all=T))
graphics.off()
set.seed(1234)

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

unit_breakeffect <- function() {
  # Parameters
  t      <- 200
  reps   <- 1000
  tauv   <- c(0.25,0.5,0.75)
  breakv <- c(-0.4,-0.2,0.2,0.4)
  phiv   <- c(1,0.95,0.9,0.85)
  cv     <- c(-3.42, -2.86, -3.13, -2.86)
  
  rej  <- matrix(0, nrow=length(phiv), ncol=4)
  s   <- matrix(nrow=reps,ncol=4)
  x   <- cbind(rep(1, t), seqa(1,1,t))
  
  # Results with no break for comparison  
  for (phic in seq(phiv)) {
    pb <- txtProgressBar(min=0, max=length(phiv), style=3)
    for (rep in seq(reps)) {
      y    <- recserar(cbind(rnorm(t)),cbind(rnorm(1)),cbind(phiv[phic]))
      uols <- glsdetrend(y,x,-t)
      ugls <- glsdetrend(y,x,-13.5)
      
      s[rep,1] <- trimr(adftests(uols,0),1,0)
      s[rep,2] <- trimr(adftests(ugls,0),1,0)
      s[rep,3] <- trimr(t(mtests(uols,0)),2,0)
      s[rep,4] <- trimr(t(mtests(ugls,0)),2,0)        
    }
    rej[phic,] <- colMeans(s < cv)    
    setTxtProgressBar(pb, phic)
  }
  close(pb)
  
  cat('\nNo Break')
  cat('\n    phiv       DF-OLS    DF-GLS    MZt-OLS   MZt-GLS\n')
  print(cbind(phiv, rej))
  
  
  # Result with break in dgp
  cat('Simulating results. This may take several minutes...')
  for (tc in seq(tauv)) {
    TB <- floor(tauv[tc]*t)
    DT <- cbind(c(rep(0, TB), seqa(1,1,t-TB)))
    
    x0 <- cbind(x, DT)
    
    for (bc in seq(breakv)) {
      beta <- rbind(0, 0, breakv[bc])
      for (phic in seq(phiv)) {
        for (rep in seq(reps)) {
          u <- recserar(cbind(rnorm(t)),cbind(rnorm(1)) ,cbind(phiv[phic]))
          y <- x0 %*% beta+u
          uols <- glsdetrend(y,x,-t)
          ugls <- glsdetrend(y,x,-13.5)
          
          s[rep,1] <- trimr(adftests(uols,0),1,0)
          s[rep,2] <- trimr(adftests(ugls,0),1,0)
          s[rep,3] <- trimr(t(mtests(uols,0)),2,0)
          s[rep,4] <- trimr(t(mtests(ugls,0)),2,0)                
        }      
        rej[phic,] <- colMeans(s < cv)          
      } 
    }
    cat('\n')
    cat('\nBreak point    = ', tauv[tc])
    cat('\nBreak size     = ', breakv[bc])
    cat('\n    phi       DF-OLS    DF-GLS    MZt-OLS   MZt-GLS')
    print(cbind(phiv, rej))
  }
  
  
}


