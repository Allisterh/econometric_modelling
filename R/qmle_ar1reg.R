#=========================================================================
#
#   Sampling properties of the White estimator 
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")

#
#------------------------- Helper Functions------------------------------
#

#load required functions - inv, trimr
source("EMTSUtil.R")

#
#--------------------  White estimation --------------------------------
#
qmle_ar1reg <- function () {
  b0    <- 0.5
  Tv    <- c(25,50,100,200,400,800,1600)
  nreps <- 10000
  tstat <- array(0, c(nreps,2))
  bhatv <- array(0, c(nreps,1))
  
  # Results
  cnames <- c("T", "Mean(bhat)", "Var(bhat)", "OLS", "White")
  tblResults <- matrix(nrow=length(Tv), ncol=length(cnames), 
                       dimnames = list(rep("", length(Tv)), cnames))
  pb <- txtProgressBar(min=0, max=length(Tv), style=3)
  for (j in seq(Tv)) {
    t <- Tv[j]
    for (k in seq(nreps)) {
      v <- rnorm(t+1)
    	y <- v
      for (i in 2:t+1) {
        y[i] <- b0*y[i-1] + sqrt(1+0.5*y[i-1]^2)*v[i] 
      }
      x     <- trimr(y,0,1) 
      y     <- trimr(y,1,0)
      xxi   <- inv(t(x) %*% x) 
  		bhat  <- xxi %*% t(x) %*% y
  		uhat  <- y-x %*% bhat
  		s2    <- t(uhat) %*% uhat/t 
      seols <- sqrt(s2 %*% xxi[1,1])
  		xuhat <- x * uhat
  
      # White's estimate of variance
  		whitevar <- xxi %*% (t(xuhat) %*% xuhat) %*% xxi 
      seW      <- sqrt(whitevar[1,1])
  
      # Save results
  		bhatv[k]   <- bhat[1]
  		tstat[k,1] <- (bhat[1]-b0)/seols
  		tstat[k,2] <- (bhat[1]-b0)/seW    
    }        
    sd_bhatv <- apply(bhatv, 2, sd)
    tblResults[j,] <- c(t, mean(bhatv), sd_bhatv^2, colMeans(abs(tstat) > 1.96))
    setTxtProgressBar(pb, j)
  }
  close(pb)
  print(tblResults)
}
