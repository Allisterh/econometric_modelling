#============================================================================
#
#   Approximate asymptotic critical values for PP coefficient test
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
#----------------------- Phillips-Perron Test ------------------------------
#
unit_ppcv <- function() {
  reps <- 100000
  t    <- 1000
  s    <- rep(0, reps)
  pb <- txtProgressBar(min=0, max=reps, style=3)
  for (j in seq(reps)) {
    b    <- cumsum(rnorm(t)/sqrt(t))
    s[j] <- 0.5*(b[t]^2-1)/mean(trimr(b,0,1)^2)    
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  cat('\n Critical Values')
  cat('\n      1%        5%        10% ')
  cat('\n', quantile(s,probs=c(0.01, 0.05, 0.1)))
  
}

