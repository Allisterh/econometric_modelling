#============================================================================
#
#   Monte Carlo experiment on the size of KPSS test
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#
# Load require functions - seqa, trimr, recserar
source("EMTSUtil.R")

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
    g[j+1] <- trimr(up,j,0) %*% trimr(up,0,j)/k
  }
  
  n  <- round(3*(k/100)^(2/25))
  s0 <- g[1] + 2*sum( g[2:(n+1)])
  s2 <- 2 * (seqa(1,1,n)^2) %*% g[2:(n+1)]
  
  ST  <- 1.3221*(abs(s2/s0)^0.4)*k^0.2
  x   <- seqa(1,1,k-1)/ST
  kQS = (25/(12*pi^2))/(x^2) * (sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5))
  
  om = (g[1]+2*kQS %*% g[2:k])/(1-a)^2
  
  return(om)
}



#
#------------- Testing the Null Hypothesis of Stationarity ------------------
#
unit_kpssmc <- function() {
  reps <- 10000
  Tv   <- c(100,200)
  arv  <- c(0,0.3,0.6,0.9)
  mav  <- c(-0.8,-0.4,0,0.4,0.8)
  
  rej <- array(0, c(length(arv),length(mav)))
  kpss <- rep(0, reps)
  
  cat('\n---------------------------------------------')
  totRows = length(Tv)*length(arv)*length(mav)
  rnames <- seq(totRows)
  cnames <- c('T', 'arv','mav', 'Rej Freq')
  tabResults <- matrix(nrow=totRows, ncol=length(arv), dimnames=list(rnames, cnames))
  row <- 1
  for (i in seq(Tv)) {
    for (j in seq(arv)) {
      for (m in seq(mav)) {
        pb <- txtProgressBar(min=0, max=reps, style=3)
        for (rep in seq(reps)) {
          v <- rnorm(Tv[i]+1,1)
          y <- recserar(cbind(trimr(v,1,0)+mav[m]*trimr(v,0,1)),cbind(v[1]),cbind(arv[j]))
          z <- y-mean(y)        
          kpss[rep] <- mean(cumsum(z)^2)/(length(z)*scQS(z))  
          setTxtProgressBar(pb, rep)
        }
        close(pb)
        rej[j,m] <- mean(kpss>0.458)
        tabResults[row, ] <- c(Tv[i], arv[j], mav[m], rej[j,m])          
        row <- row+1     
      }
      
    }
  }
  print(tabResults)
}

