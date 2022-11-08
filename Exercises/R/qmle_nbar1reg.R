#========================================================================
#
#   Program to compute sampling properties of quasi-maximum likelihood 
#   estimator where true distribution is a negative binomial and the 
#   misspecified distribution is Poisson
#
#========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345678, kind="Mersenne-Twister")

#
#------------------------- Helper Functions-----------------------------
#

#load required functions - inv, trimr
source("EMTSUtil.R")

#-----------------------------------------------------------------------
# Negative unconstrained log-likelihood
#-----------------------------------------------------------------------
neglog <- function(b,y,x) {
  b[1] <- 1/(1+exp(-b[1]))  
  b[2] <- b[2]^2
  lf   <- -mean(y*log(x %*% b)-(x %*% b))  
  return(lf)  
}

#
#---------------------- Sampling Properties of QMLE -------------------
#
  
qmle_nbar1reg <- function() {
    Tv <- c(25,50,100,200,400,800,1600)
  b0    <- 0.5
  nreps <- 100
  tstat <- array(0, c(nreps,2))
  tp    <- array(0, c(nreps,2))
  bnv   <- array(0, c(nreps,1))
  bpv   <- array(0, c(nreps,1))

  # Results
  cnames <- c("T", "meanN", "meanP", "varN", "varP",  "tOLS", "tOLSW", "tPH", "tPHJH")
  tblResults <- matrix(nrow=length(Tv), ncol=length(cnames), 
                       dimnames = list(rep("", length(Tv)), cnames))  
  pb <- txtProgressBar(min=0, max=length(Tv), style=3)
  for (k in seq(Tv) ) {
    t <- Tv[k]
    for (j in seq(nreps)) {
      y <- rep(1, (t+1))
      for (i in 2:t+1) {
        # Generate negative binomial random numbers
        y[i] <- rnbinom(1, size=1+b0*y[i-1], prob=0.5) 
      }
      x <- cbind(trimr(y,0,1), rep(1, t))
      y <- trimr(y,1,0)     
            
      xxi   <- inv(t(x) %*% x) 
      bhat  <- xxi %*% t(x) %*% y
      uhat  <- y- x %*% bhat      
      s2    <- t(uhat) %*% uhat/t 
      seols <- sqrt(s2 %*% xxi[1,1])
      
      xuhat <- x * c(uhat)
      
      # White's estimate of variance
      whitevar <- xxi %*% (t(xuhat) %*% xuhat) %*% xxi 
      seW      <- sqrt(whitevar[1,1])
  
      # Save results
      bnv[j] <- bhat[1]
      tstat[j,1] <- (bhat[1]-b0)/seols
      tstat[j,2] <- (bhat[1]-b0)/seW
      
      ba <- c(log(b0/(1-b0)), 1)      
      estResults <- optim(ba, neglog, y=y, x=x, method="BFGS")
      bp <- estResults$par
      
      bp[1] <- 1/(1+exp(-bp[1])) 
      bp[2] <- bp[2]^2      
      
      # right array divide
      xdxb <- x / c(x %*% bp)      
      
      Hi   <- -inv(t(x) %*% xdxb)       
      
      seHi <- sqrt(abs(Hi[1,1]))
      g    <- xdxb * c(y-x %*% bp)
      
      J <- t(g) %*% g
      
      HiJHi <- Hi %*% J %*% Hi 
      seHJH <- sqrt(HiJHi[1,1])
      
      bpv[j]  <- bp[1]
      tp[j,1] <- (bp[1]-b0)/seHi
      tp[j,2] <- (bp[1]-b0)/seHJH  
    }    
    tblResults[k,] <- c(t, mean(bnv), mean(bpv), var(bnv), var(bpv), 
                       colMeans(abs(tstat)>1.96), colMeans(abs(tp) > 1.96))  
    setTxtProgressBar(pb, k)    
  }
  close(pb)
  print(tblResults)  
}



