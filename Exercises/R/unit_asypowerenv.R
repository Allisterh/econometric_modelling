#=========================================================================
#
#   Approximate asymptotic power envelope and power curves 
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# load required functions -  figure, recserar, trimr, inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Detrending function: 
#       cbar = -7 constant 
#       cbar = -13.5 linear trend  
#       cbar = -T for OLS detrending
#----------------------------------------------------------------------------

glsdetrend <- function(y,x,cbar) {
  t <- length(y)
  yc <- cbind(c(y[1], (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1))) )
  xc <- rbind(x[1,], cbind(trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)))    
  b <- lm(yc ~ xc - 1)$coef
  u <- y - x %*% b  
  return(u)
}

#----------------------------------------------------------------------------
#  ADF coefficient and t tests: u must be already de-trended
#----------------------------------------------------------------------------
adftests <- function(u,k) {
  du <- trimr(u,1,0) - trimr(u,0,1)
  x  <- trimr(u,0,1)
  
  xxi <- inv(t(x) %*% x) 
  b   <- xxi %*% t(x) %*% trimr(du,k,0) 
  e   <- trimr(du,k,0)- x %*% b   
  s2  <- (t(e) %*% e)/length(e)
  
  adfstats <- cbind(c(length(u)*b[1], b[1]/sqrt(s2 %*% xxi[1,1])) )
  
  return(adfstats)
}


#----------------------------------------------------------------------------
#  P-test
#----------------------------------------------------------------------------

Ptest <- function(y,x,cbar,k) {
   n  <- length(y)
   uc <- glsdetrend(y,x,cbar) 
   u0 <- glsdetrend(y,x,0)
   s2 <- ar1rvar(uc,k)
   uc <- cbind(c(uc[1], trimr(uc,1,0)-(1+cbar/n)*trimr(uc,0,1)) )
   u0 <- cbind(c(u0[1], trimr(u0,1,0)-trimr(u0,0,1)) )
   pt <- (t(uc) %*% uc-(1+cbar/n) %*% t(u0) %*% u0)/s2   
   return(pt)
}

#----------------------------------------------------------------------------
#  ar1rvar
#----------------------------------------------------------------------------

ar1rvar <- function(u,k) {
  du <- trimr(u,1,0)-trimr(u,0,1) 
  x  <- cbind(trimr(u,0,1))
  
  b <- lm(trimr(du,k,0) ~ x - 1)$coef
  e <- trimr(du,k,0) - x %*% b 
  s2 <- (t(e) %*% e)/length(e)
  if (k > 0) 
    s2 <- s2/(1-sum(trimr(b,1,0)))^2    
  return (s2)
}

#----------------------------------------------------------------------------
# Compute and return the rejection frequencies for 2 matrix inputs
#   compute the quantiles of M1
#   check if M2 < quatiles
# M1 and M2 must be matrices
# prob is probability value
#----------------------------------------------------------------------------
reject.freq <- function(M1, M2, prob=0.05) {
  tmp.quantiles <- apply(M1, 2,  quantile, probs=prob)
  tmp.lt <- apply(M2, 1, '<', tmp.quantiles)
  rMeans <- rowMeans(tmp.lt)
  return(rMeans)
}


#
#----------------------- Asymptotic Power Envelope --------------------------
#

unit_asypowerenv <- function() {
  # Critical values to compute the power envelope with cv = 0 -> size of 0.05           
  cv   <- seq(from=-30, by=1, length.out=30)  
  n    <- length(cv)
  t    <- 1000                       
  nreps <- 50000
  
  # Detrending parameters
  x    <- cbind(rep(1, t), seq(from=1, by=1, length.out=t))
  cbar <--13.5                     
  
  # Allocate memory
  zerosn <- array(0, c(nreps,n) )
  zeros1 <- rep(0, nreps)
  pcc    <- zerosn
  pc0    <- zerosn
  pcbarc <- zerosn 
  pcbar0 <- zeros1
  dfols  <- zerosn 
  dfols0 <- zeros1
  dfgls  <- zerosn 
  dfgls0 <- zeros1

  pb <- txtProgressBar(min=0, max=n, style=3)
  for (k in seq(n)) {
    set.seed(1234, kind="Mersenne-Twister")    
    c <- cv[k]   
    for (j in seq(nreps)) {      
      u  <- rnorm(t)
      
      yc <- recserar(cbind(u),cbind(u[1]),cbind(1+c/t))
      y0 <- recserar(cbind(u),cbind(u[1]),cbind(1))
      
      pcc[j,k]    <- Ptest(yc,x,c,0)
      
      pc0[j,k]    <- Ptest(y0,x,c,0)
      pcbarc[j,k] <- Ptest(yc,x,cbar,0)
      dfols[j,k]  <- trimr(adftests(glsdetrend(yc,x,-t),0),1,0)         
      dfgls[j,k]  <- trimr(adftests(glsdetrend(yc,x,cbar),0),1,0)  
            
      if (k == 1) {
        pcbar0[j] <- Ptest(y0,x,cbar,0) 
        dfols0[j] <- trimr(adftests(glsdetrend(y0,x,-t),0),1,0)
        dfgls0[j] <- trimr(adftests(glsdetrend(y0,x,cbar),0),1,0)        
      } 
           
    } 
    setTxtProgressBar(pb, k)    
  }
  close(pb)
  

  # Compute rejection frequencies for alternative detrending methods 
  rejpc    <- cbind(c( reject.freq(pc0, pcc), 0.05))

  pcbar0.quantiles <- quantile(pcbar0, probs=0.05)
  rejpcbar <- cbind(c( colMeans(pcbarc < pcbar0.quantiles),  0.05 ))
  
  dfols0.quantiles <- quantile(dfols0, probs=0.05)
  rejdfols <- cbind(c( colMeans(dfols < dfols0.quantiles),  0.05 ))
  
  dfgls.quantiles <- quantile(dfgls0, probs=0.05)
  rejdfgls <- cbind(c( colMeans(dfgls < dfgls.quantiles),  0.05 ))
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************

  figure()
  dvec <- cbind(c(cv, 0))

  #--------------------------------------------------------#
  matplot(dvec,cbind(rejpc,rejpcbar,rejdfols,rejdfgls),type="l",
          ylab = "Power",
          xlab = "Critical Values",
          xlim = c(-30, 0),
          ylim = c(0, 1),
          bty = "l")
}


