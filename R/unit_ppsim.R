#=============================================================================
#
#   Monte Carlo experiment on the PP coefficient test
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#
# load required functions -  trimr
source("EMTSUtil.R")

#
#--------------------------- Functions -----------------------------------
# 
#-------------------------------------------------------------------------
#  Long-run variance
#-------------------------------------------------------------------------
scQS <- function(u) {
  a <- lm(trimr(u,1,0) ~ trimr(u,0,1) - 1)$coef
  
  if (abs(a)>0.97) 
    a <- 0.97*a/abs(a)  
  
  up <- trimr(u,1,0)-trimr(u,0,1)*a
  k  <- length(up)
  g  <- rep(0, k)
  
  for (j in 0:(k-1)) {
    g[j+1] <- t( trimr(up,j,0) ) %*% trimr(up,0,j)/k
  }
 
  n  <- round(3*(k/100)^(2/25))
  s0 <- g[1]+2*sum(g[2:(n+1)])
  s2 <- 2*(seqa(1,1,n)^2) %*% g[2:(n+1)]
  
  ST  <- 1.3221*(abs(s2/s0)^0.4)*k^0.2
  x   <- seqa(1,1,k-1)/ST
  kQS <- (25/(12*pi^2))/(x^2)*(sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5))
  om <- (g[1]+2*kQS %*% g[2:k])/(1-a)^2  
  return(om)
}



#
#----------------------- Phillips-Perron Test ------------------------------
#
unit_ppsim <- function( ) {
  reps <- 10000
  Tv   <- c(100,200)
  arv  <- c(0,0.3,0.6,0.9)
  mav  <- c(-0.8,-0.4,0,0.4,0.8)
  
  rej <- array(0, c(length(arv),length(mav)))
  phitilde <- rep(0, reps)
  
  cat('\n---------------------------------------------\n')
  for (i in 1:length(Tv)) {
    for (j in 1:length(arv)) {
      for (m in 1:length(mav)) {
        set.seed(1234)
        
        for (k in seq(reps)) {
          v <- rnorm(Tv[i]+1)
          y <- cumsum(recserar(cbind(trimr(v,1,0)+mav[m]*trimr(v,0,1)),cbind(v[1]),cbind(arv[j])))
          
          phihat <- lm(trimr(y,1,0) ~ trimr(y,0,1) - 1)$coef
          vhat   <- trimr(y,1,0)-phihat * trimr(y,0,1)
          
          sig2   <- t(vhat) %*% vhat/length(vhat)
          om2    <- scQS(vhat)
          
          phitilde[k] <- phihat - 0.5*(om2-sig2)/mean(trimr(y,0,1)^2)          
        }     
        rej[j,m] <- mean(length(y)*(phitilde-1) < -8.08)              
        print(cbind("T"=Tv[i], arv=arv[j], mav=mav[m], "Rej Freq"=rej[j,m]))
      }
    }
  }  
}

